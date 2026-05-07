#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "three_pcf.h"
#include "../mracs/sfc.h"

int three_pcf(const Particle *particles, long long n,
              WaveletConfig *wav, KernelType kernel,
              double radius, double boxsize,
              int nbins_r, Result *three_pcf_out)
{
    (void) particles;
    (void) n;
    (void) wav;
    (void) kernel;
    (void) radius;
    (void) boxsize;
    (void) nbins_r;
    (void) three_pcf_out;

    fprintf(stderr, "three_pcf: FFTW + MRA integration required — "
                    "skeleton only, full implementation deferred\n");
    return -1;
}
