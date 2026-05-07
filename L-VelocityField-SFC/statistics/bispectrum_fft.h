#ifndef BISPECTRUM_FFT_H
#define BISPECTRUM_FFT_H

#include "../allvars.h"

/* Result structure for binned statistics */
typedef struct {
    long long n;
    double x;
    double y;
    double data;
    double data2;
} Result;

/* FFT-based bispectrum computation on an already-FFT'd grid.
 * grid: complex FFT data (in-place layout)
 * ngrid: grid dimension
 * boxsize: simulation box size in Mpc/h
 * nbins: number of k bins
 * output: pre-allocated Result arrays */
void bispectrum_fft(fftw_complex *grid, int ngrid, double boxsize,
                    int nbins, Result *bispec_out, Result *p1_out,
                    Result *p2_out, Result *p3_out);

#endif
