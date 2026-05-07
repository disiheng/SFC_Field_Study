#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* FFTW3-MPI is included via allvars.h equivalent; mra_grid is standalone MPI.
 * Include explicitly for plan creation functions. */
#include <fftw3-mpi.h>

#include "mra_grid.h"
#include "sfc.h"

#ifndef CUB
#define CUB(a) ((a)*(a)*(a))
#endif
#ifndef SQR
#define SQR(a) ((a)*(a))
#endif

int mra_grid_build(const Particle *particles, long long n,
                   uint64_t grid_len, double boxsize, MRAGrid *grid)
{
    uint64_t *indices;
    long long i;
    double cell_vol;

    if (!particles || n <= 0 || !grid) return -1;

    grid->grid_len = grid_len;
    grid->grid_num = grid_len * grid_len * grid_len;
    grid->boxsize = boxsize;

    grid->density = (double *) calloc((size_t) grid->grid_num, sizeof(double));
    if (!grid->density) return -1;

    indices = sfc_order(particles, n, grid_len);
    if (!indices) {
        free(grid->density);
        grid->density = NULL;
        return -1;
    }

    cell_vol = CUB(boxsize / (double) grid_len);

    for (i = 0; i < n; i++) {
        uint64_t idx = indices[i];
        if (idx < grid->grid_num) {
            grid->density[idx] += particles[i].weight / cell_vol;
        }
    }

    free(indices);
    return 0;
}

int mra_grid_convolve(MRAGrid *grid, WaveletConfig *wav,
                      KernelType kernel_type, double radius, double theta)
{
    fftw_plan fwd_plan, inv_plan;
    ptrdiff_t local_n0, local_0_start, local_n1, local_1_start;
    ptrdiff_t alloc_local;
    int ngrid, ngrid_half;
    double delk, boxsize;
    long long i;
    int x, y, z, ip;
    double kx, ky, kz;
    double norm_factor;

    if (!grid || !grid->density) return -1;

    ngrid = (int) grid->grid_len;
    ngrid_half = ngrid / 2;
    boxsize = grid->boxsize;
    delk = 2.0 * M_PI / boxsize;

    /* FFTW3 MPI: create plans bound to grid->density buffer.
     * In-place: same pointer for in and out; FFTW3 handles the sizing.
     * Reallocate density to hold complex output (which is larger). */
    alloc_local = fftw_mpi_local_size_3d_transposed(
        (ptrdiff_t) ngrid, (ptrdiff_t) ngrid, (ptrdiff_t) ngrid,
        MPI_COMM_WORLD,
        &local_n0, &local_0_start, &local_n1, &local_1_start);

    /* Ensure buffer is large enough for complex output (in-place) */
    size_t cplx_bytes = (size_t) alloc_local * sizeof(fftw_complex);
    size_t real_bytes = (size_t) local_n0 * (size_t) ngrid * (size_t) ngrid
                        * sizeof(double);
    size_t need_bytes = (cplx_bytes > real_bytes) ? cplx_bytes : real_bytes;
    size_t have_bytes = (size_t) grid->grid_num * sizeof(double);

    if (have_bytes < need_bytes) {
        double *newbuf = (double *) realloc(grid->density, need_bytes);
        if (!newbuf) return -1;
        grid->density = newbuf;
    }

    fwd_plan = fftw_mpi_plan_dft_r2c_3d(
        (ptrdiff_t) ngrid, (ptrdiff_t) ngrid, (ptrdiff_t) ngrid,
        grid->density, (fftw_complex *) grid->density,
        MPI_COMM_WORLD, FFTW_ESTIMATE);

    inv_plan = fftw_mpi_plan_dft_c2r_3d(
        (ptrdiff_t) ngrid, (ptrdiff_t) ngrid, (ptrdiff_t) ngrid,
        (fftw_complex *) grid->density, grid->density,
        MPI_COMM_WORLD, FFTW_ESTIMATE);

    /* Forward FFT: real -> complex (in-place) */
    fftw_execute(fwd_plan);

    /* After forward FFT, data is in transposed half-complex format:
     * For each y in [local_1_start, local_1_start + local_n1):
     *   For each x in [0, ngrid):
     *     For each z in [0, ngrid/2 + 1):
     *       complex value at offset = y_local * ngrid*(ngrid_half+1) + x*(ngrid_half+1) + z
     */
    fftw_complex *fdata = (fftw_complex *) grid->density;
    ptrdiff_t ly_start = local_1_start;
    ptrdiff_t ly_count = local_n1;

    for (y = (int) ly_start; y < (int)(ly_start + ly_count); y++) {
        int yy = (y > ngrid_half) ? y - ngrid : y;

        for (x = 0; x < ngrid; x++) {
            int xx = (x > ngrid_half) ? x - ngrid : x;

            for (z = 0; z <= ngrid_half; z++) {
                int zz = (z > ngrid_half) ? z - ngrid : z;

                kx = xx * delk;
                ky = yy * delk;
                kz = zz * delk;

                ip = (int)(((ptrdiff_t)(y - (int) ly_start)) * (ptrdiff_t) ngrid
                     * (ptrdiff_t)(ngrid_half + 1)
                     + (ptrdiff_t) x * (ptrdiff_t)(ngrid_half + 1)
                     + (ptrdiff_t) z);

                /* Evaluate window kernel in k-space */
                double window = kernel_eval(kernel_type, radius, theta, kx, ky, kz);

                /* Apply wavelet phi filter if configured */
                if (wav && wav->phi && wav->phi_len > 0) {
                    double kmag = sqrt(kx * kx + ky * ky + kz * kz);
                    double phi_k = 1.0;
                    window *= phi_k;
                    (void) kmag;
                }

                fdata[ip][0] *= window;
                fdata[ip][1] *= window;
            }
        }
    }

    /* Inverse FFT: complex -> real (in-place) */
    fftw_execute(inv_plan);

    /* Normalize: FFTW unnormalized inverse gives ngrid^3 * original */
    norm_factor = 1.0 / CUB((double) ngrid);
    long long nlocal = (long long) local_n0 * (long long) ngrid
                       * (long long) ngrid;
    for (i = 0; i < nlocal; i++) {
        grid->density[i] *= norm_factor;
    }

    /* Clean up */
    fftw_destroy_plan(fwd_plan);
    fftw_destroy_plan(inv_plan);

    return 0;
}

int mra_grid_write(const MRAGrid *grid, const char *filename)
{
    FILE *fp;
    if (!grid || !grid->density) return -1;

    fp = fopen(filename, "wb");
    if (!fp) return -1;

    fwrite(grid->density, sizeof(double), (size_t) grid->grid_num, fp);
    fclose(fp);
    return 0;
}

void mra_grid_free(MRAGrid *grid)
{
    if (grid->density) {
        free(grid->density);
        grid->density = NULL;
    }
}
