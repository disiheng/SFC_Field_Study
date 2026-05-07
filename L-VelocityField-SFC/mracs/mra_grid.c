#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef NOTYPEPREFIX_FFTW
#include <rfftw_mpi.h>
#else
#ifdef DOUBLEPRECISION_FFTW
#include <drfftw_mpi.h>
#else
#include <srfftw_mpi.h>
#endif
#endif

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
    rfftwnd_mpi_plan fwd_plan, inv_plan;
    int local_nx, local_x_start, local_ny, local_y_start, total_local_size;
    fftw_real *workspace;
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

    /* Create FFTW plans for MRA grid dimensions */
    fwd_plan = rfftw3d_mpi_create_plan(MPI_COMM_WORLD,
                         ngrid, ngrid, ngrid,
                         FFTW_REAL_TO_COMPLEX,
                         FFTW_ESTIMATE | FFTW_IN_PLACE);

    inv_plan = rfftw3d_mpi_create_plan(MPI_COMM_WORLD,
                         ngrid, ngrid, ngrid,
                         FFTW_COMPLEX_TO_REAL,
                         FFTW_ESTIMATE | FFTW_IN_PLACE);

    rfftwnd_mpi_local_sizes(fwd_plan, &local_nx, &local_x_start,
                            &local_ny, &local_y_start, &total_local_size);

    /* Allocate workspace for FFTW transpose */
    workspace = (fftw_real *) malloc(total_local_size * sizeof(fftw_real));
    if (!workspace) {
        rfftwnd_mpi_destroy_plan(fwd_plan);
        rfftwnd_mpi_destroy_plan(inv_plan);
        return -1;
    }

    /* Forward FFT: real -> complex (in-place) */
    rfftwnd_mpi(fwd_plan, 1, grid->density, workspace, FFTW_TRANSPOSED_ORDER);

    /* After forward FFT, data is in half-complex format:
     * For each y in [local_y_start, local_y_start + local_ny):
     *   For each x in [0, ngrid):
     *     For each z in [0, ngrid/2 + 1):
     *       complex value at ip = ngrid*(ngrid_half+1)*y_local + (ngrid_half+1)*x + z
     */
    fftw_complex *fdata = (fftw_complex *) grid->density;

    for (y = local_y_start; y < local_y_start + local_ny; y++) {
        int yy = (y > ngrid_half) ? y - ngrid : y;

        for (x = 0; x < ngrid; x++) {
            int xx = (x > ngrid_half) ? x - ngrid : x;

            for (z = 0; z <= ngrid_half; z++) {
                int zz = (z > ngrid_half) ? z - ngrid : z;

                kx = xx * delk;
                ky = yy * delk;
                kz = zz * delk;

                ip = (y - local_y_start) * ngrid * (ngrid_half + 1)
                   + x * (ngrid_half + 1) + z;

                /* Evaluate window kernel in k-space */
                double window = kernel_eval(kernel_type, radius, theta, kx, ky, kz);

                /* Apply wavelet phi filter if configured */
                if (wav && wav->phi && wav->phi_len > 0) {
                    /* Evaluate phi in k-space: FFT of scaling function.
                     * For simplicity, use the scaling function's Fourier transform
                     * approximated at the k-magnitude.
                     * Full implementation needs phi's k-space representation. */
                    double kmag = sqrt(kx * kx + ky * ky + kz * kz);
                    double phi_k = 1.0;  /* placeholder — phi k-space eval deferred */
                    window *= phi_k;
                    (void) kmag;  /* suppress unused warning */
                }

                fdata[ip].re *= window;
                fdata[ip].im *= window;
            }
        }
    }

    /* Inverse FFT: complex -> real (in-place) */
    rfftwnd_mpi(inv_plan, 1, grid->density, workspace, FFTW_TRANSPOSED_ORDER);

    /* Normalize: FFTW unnormalized inverse gives ngrid^3 * original */
    norm_factor = 1.0 / CUB((double) ngrid);
    long long nlocal = (long long) local_nx * ngrid * ngrid;
    for (i = 0; i < nlocal; i++) {
        grid->density[i] *= norm_factor;
    }

    /* Clean up */
    free(workspace);
    rfftwnd_mpi_destroy_plan(fwd_plan);
    rfftwnd_mpi_destroy_plan(inv_plan);

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
