#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "bispectrum_wavelet.h"

int bispectrum_wavelet(const MRAGrid *grid, WaveletConfig *wav,
                       KernelType kernel, double radius, int nbins,
                       Result *bispec_out)
{
    (void) grid;
    (void) wav;
    (void) kernel;
    (void) radius;
    (void) nbins;
    (void) bispec_out;

    fprintf(stderr, "bispectrum_wavelet: FFTW integration required — "
                    "skeleton only, full implementation deferred\n");
    return -1;
}
