#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "mra_grid.h"
#include "sfc.h"

#ifndef CUB
#define CUB(a) ((a)*(a)*(a))
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
    /* TODO: FFT-based 3D convolution using FFTW2.
     * This requires the MPI FFTW context which is initialized in velocity_field.c.
     * The full implementation will:
     * 1. Forward FFT of grid->density
     * 2. Multiply by kernel in k-space
     * 3. Multiply by wavelet phi in k-space
     * 4. Inverse FFT
     * Deferred to the integration step when FFTW context is available.
     */
    (void) grid;
    (void) wav;
    (void) kernel_type;
    (void) radius;
    (void) theta;
    fprintf(stderr, "mra_grid_convolve: not yet implemented (needs FFTW context)\n");
    return -1;
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
