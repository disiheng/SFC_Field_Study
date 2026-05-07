#ifndef THREE_PCF_H
#define THREE_PCF_H

#include "../mracs/particle_io.h"
#include "../mracs/mra_grid.h"
#include "../mracs/kernel.h"
#include "bispectrum_fft.h"  /* for Result */

/* 3-point correlation function using MRA wavelet decomposition.
 * Operates directly on particle catalog. */
int three_pcf(const Particle *particles, long long n,
              WaveletConfig *wav, KernelType kernel,
              double radius, double boxsize,
              int nbins_r, Result *three_pcf_out);

#endif
