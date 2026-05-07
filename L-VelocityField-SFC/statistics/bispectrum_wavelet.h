#ifndef BISPECTRUM_WAVELET_H
#define BISPECTRUM_WAVELET_H

#include "../mracs/mra_grid.h"
#include "../mracs/kernel.h"
#include "bispectrum_fft.h"  /* for Result */

/* Wavelet-based bispectrum using pre-built MRA grid.
 * Decomposes the density field via wavelet transform, then
 * computes bispectrum of wavelet coefficients at each scale. */
int bispectrum_wavelet(const MRAGrid *grid, WaveletConfig *wav,
                       KernelType kernel, double radius, int nbins,
                       Result *bispec_out);

#endif
