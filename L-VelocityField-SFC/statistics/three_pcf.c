#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "three_pcf.h"
#include "../mracs/sfc.h"
#include "../mracs/mra_grid.h"

int three_pcf(const Particle *particles, long long n,
              WaveletConfig *wav, KernelType kernel,
              double radius, double boxsize,
              int nbins_r, Result *three_pcf_out)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (!particles || n <= 0 || !three_pcf_out) return -1;
    if (nbins_r <= 0) return -1;

    /* Clear output bins */
    int i;
    for (i = 0; i < nbins_r; i++) {
        three_pcf_out[i].n = 0;
        three_pcf_out[i].data = 0.0;
        three_pcf_out[i].data2 = 0.0;
        three_pcf_out[i].x = 0.0;
        three_pcf_out[i].y = 0.0;
    }

    /* Resolution: grid_len = 2^J, using wav->resolution or default J=8 */
    uint64_t grid_len = (wav && wav->resolution > 0)
                        ? (1ULL << wav->resolution) : (1ULL << 8);

    /* Step 1: Build MRA density grid */
    MRAGrid mgrid;
    int ret = mra_grid_build(particles, n, grid_len, boxsize, &mgrid);
    if (ret != 0) {
        if (rank == 0) fprintf(stderr, "three_pcf: grid build failed\n");
        return -1;
    }
    if (rank == 0) printf("3PCF: Built MRA grid (%llu cells)\n", mgrid.grid_num);

    /* Step 2: Multi-scale smoothing via convolution */
    int nscales = 3;
    int scale;
    double dr = radius / (double) nbins_r;

    for (scale = 0; scale < nscales; scale++) {
        double R_scale = radius * pow(1.5, (double) scale);

        /* Copy density for this scale */
        MRAGrid work_grid;
        work_grid.grid_len = mgrid.grid_len;
        work_grid.grid_num = mgrid.grid_num;
        work_grid.boxsize = mgrid.boxsize;
        work_grid.density = (double *) malloc(
            (size_t) mgrid.grid_num * sizeof(double));
        if (!work_grid.density) continue;

        memcpy(work_grid.density, mgrid.density,
               (size_t) mgrid.grid_num * sizeof(double));

        ret = mra_grid_convolve(&work_grid, wav, kernel, R_scale, 0.0);
        if (ret != 0) {
            mra_grid_free(&work_grid);
            continue;
        }

        if (rank == 0) {
            printf("3PCF: Scale %d (R=%.2f Mpc/h) smoothed\n", scale, R_scale);
        }

        /* Step 3: Sample the smoothed field to estimate 3PCF.
         * Q(r) = <delta(x)delta(x+r1)delta(x+r2)> / <delta^2>^2
         *
         * Approximate by sampling triples from the grid.
         * Use a subset of cells to keep computation tractable.
         */
        double cell_size = boxsize / (double) grid_len;
        long long total_cells = (long long) mgrid.grid_num;
        long long sample_step = total_cells / 10000;
        if (sample_step < 1) sample_step = 1;

        long long j;
        double sum_corr = 0.0, count = 0.0;

        for (j = 0; j < total_cells; j += sample_step) {
            /* Get the density at this cell */
            double d1 = work_grid.density[j];

            /* Compute cell indices */
            uint64_t iz = j % grid_len;
            uint64_t iy = (j / grid_len) % grid_len;
            uint64_t ix = j / (grid_len * grid_len);

            /* Find a neighbor at distance ~r by stepping cells */
            int r_cells = (int)(radius / cell_size);
            if (r_cells < 1) r_cells = 1;

            uint64_t jx = (ix + r_cells) % grid_len;
            uint64_t jy = (iy + r_cells / 2) % grid_len;
            uint64_t jz = iz;  /* same z-plane */

            long long idx2 = (long long)(jx * grid_len * grid_len
                                        + jy * grid_len + jz);

            uint64_t kx = (ix + r_cells / 3) % grid_len;
            uint64_t ky = (iy + r_cells) % grid_len;
            uint64_t kz = (iz + r_cells / 2) % grid_len;

            long long idx3 = (long long)(kx * grid_len * grid_len
                                        + ky * grid_len + kz);

            if (idx2 >= 0 && idx2 < total_cells
                && idx3 >= 0 && idx3 < total_cells) {
                double d2 = work_grid.density[idx2];
                double d3 = work_grid.density[idx3];
                sum_corr += d1 * d2 * d3;
                count += 1.0;
            }
        }

        double global_sum_corr, global_count;
        MPI_Allreduce(&sum_corr, &global_sum_corr, 1, MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);
        MPI_Allreduce(&count, &global_count, 1, MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);

        double q3 = (global_count > 0.0) ? global_sum_corr / global_count : 0.0;

        /* Bin by scale radius */
        int bin = (int)(R_scale / dr);
        if (bin >= 0 && bin < nbins_r) {
            three_pcf_out[bin].n += 1;
            three_pcf_out[bin].data += q3;
            three_pcf_out[bin].data2 += q3 * q3;
            three_pcf_out[bin].x = R_scale;
        }

        mra_grid_free(&work_grid);
    }

    /* Set bin centers */
    for (i = 0; i < nbins_r; i++) {
        if (three_pcf_out[i].x == 0.0) {
            three_pcf_out[i].x = dr * ((double)i + 0.5);
        }
    }

    mra_grid_free(&mgrid);
    return 0;
}
