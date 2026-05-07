#ifndef WAVELET_H
#define WAVELET_H

typedef struct {
    int resolution;       /* J, grid side = 2^J */
    int base_type;        /* 0 = B-Spline, 1 = Daubechies */
    int genus;            /* Daubechies order n, or B-Spline order */
    int samp_rate;        /* samples per unit support interval */
    int phi_len;          /* total samples in phi array */
    double *phi;          /* scaling function samples */
} WaveletConfig;

/* Load Daubechies phi from pre-computed .bin file */
int daubechies_load(const char *datadir, int genus, WaveletConfig *cfg);

/* Generate B-Spline scaling function in memory (does not require .bin file) */
int bspline_generate(int order, int samp_rate, WaveletConfig *cfg);

/* Free wavelet config */
void wavelet_free(WaveletConfig *cfg);

#endif
