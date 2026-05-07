#ifndef MRA_GRID_H
#define MRA_GRID_H

#include <stdint.h>
#include "particle_io.h"
#include "wavelet.h"
#include "kernel.h"

/* MRA density grid built from particles via SFC ordering */
typedef struct {
    uint64_t grid_len;    /* side = 2^J */
    uint64_t grid_num;    /* total cells = grid_len^3 */
    double *density;      /* density field on grid */
    double boxsize;
} MRAGrid;

/* Build density grid from particles using SFC ordering */
int mra_grid_build(const Particle *particles, long long n,
                   uint64_t grid_len, double boxsize, MRAGrid *grid);

/* 3D convolution of grid with window kernel in Fourier space.
 * Uses user-provided FFTW plan context (caller must init FFTW).
 * Output overwrites grid->density with convolved field. */
int mra_grid_convolve(MRAGrid *grid, WaveletConfig *wav,
                      KernelType kernel_type, double radius, double theta);

/* Write an MRAGrid to disk as raw doubles */
int mra_grid_write(const MRAGrid *grid, const char *filename);

/* Free MRAGrid */
void mra_grid_free(MRAGrid *grid);

#endif
