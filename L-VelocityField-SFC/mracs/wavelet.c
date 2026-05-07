#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "wavelet.h"

int daubechies_load(const char *datadir, int genus, WaveletConfig *cfg)
{
    char fname[1024];
    int phi_start = 0;
    int phi_end = 2 * genus - 1;
    int phi_support = phi_end - phi_start;
    int nsamp, i;

    snprintf(fname, sizeof(fname), "%s/DaubechiesG%dPhi.bin", datadir, genus);

    FILE *fp = fopen(fname, "rb");
    if (!fp) {
        fprintf(stderr, "Error: cannot open wavelet file %s\n", fname);
        return -1;
    }

    nsamp = cfg->samp_rate * phi_support;
    cfg->phi = (double *) malloc(nsamp * sizeof(double));
    if (!cfg->phi) {
        fprintf(stderr, "Error: malloc failed for phi (%d doubles)\n", nsamp);
        fclose(fp);
        return -1;
    }

    size_t nread = fread(cfg->phi, sizeof(double), nsamp, fp);
    fclose(fp);

    if (nread != (size_t) nsamp) {
        fprintf(stderr, "Error: read %zu/%d values from %s\n", nread, nsamp, fname);
        free(cfg->phi);
        cfg->phi = NULL;
        return -1;
    }

    cfg->genus = genus;
    cfg->base_type = 1;  /* Daubechies */
    cfg->phi_len = nsamp;

    return 0;
}

int bspline_generate(int order, int samp_rate, WaveletConfig *cfg)
{
    int support = order + 1;
    int nsamp = support * samp_rate;
    int i, k;
    double dx = 1.0 / (double) samp_rate;

    cfg->phi = (double *) calloc(nsamp, sizeof(double));
    if (!cfg->phi) return -1;

    /* Generate order-0 B-spline (box function) */
    for (i = 0; i < samp_rate; i++) {
        cfg->phi[i] = 1.0;
    }

    /* Recursive convolution to build order-n B-spline */
    for (k = 1; k <= order; k++) {
        double *tmp = (double *) calloc(nsamp, sizeof(double));
        if (!tmp) { free(cfg->phi); return -1; }

        for (i = 0; i < nsamp; i++) {
            double sum = 0.0;
            int j;
            for (j = 0; j < samp_rate; j++) {
                int idx = i - j;
                if (idx >= 0 && idx < nsamp) {
                    sum += cfg->phi[idx];
                }
            }
            tmp[i] = sum / (double) samp_rate;
        }

        /* Normalize */
        double norm = 0.0;
        for (i = 0; i < nsamp; i++) norm += tmp[i];
        norm *= dx;
        if (norm > 0.0) {
            for (i = 0; i < nsamp; i++) tmp[i] /= norm;
        }

        free(cfg->phi);
        cfg->phi = tmp;
    }

    cfg->genus = order;
    cfg->base_type = 0;  /* B-Spline */
    cfg->samp_rate = samp_rate;
    cfg->phi_len = nsamp;

    return 0;
}

void wavelet_free(WaveletConfig *cfg)
{
    if (cfg->phi) {
        free(cfg->phi);
        cfg->phi = NULL;
    }
    cfg->phi_len = 0;
}
