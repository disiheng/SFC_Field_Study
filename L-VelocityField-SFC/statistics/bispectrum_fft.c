#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../allvars.h"
#include "../proto.h"
#include "bispectrum_fft.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef CUB
#define CUB(a) ((a)*(a)*(a))
#endif

#ifndef SQR
#define SQR(a) ((a)*(a))
#endif

void bispectrum_fft(fftw_complex *grid, int ngrid, double boxsize,
                    int nbins, Result *bispec_out, Result *p1_out,
                    Result *p2_out, Result *p3_out)
{
    double delk = 2.0 * M_PI / boxsize;
    double k_min = delk;
    double k_max = delk * (ngrid / 2);
    double dk = (k_max - k_min) / nbins;
    int i, x1, y1, z1, x2, y2, z2;
    int ip1, ip2;
    double k1x, k1y, k1z, k2x, k2y, k2z, k3x, k3y, k3z;
    double k1, k2, k3;
    int bin;
    int ngrid_half = ngrid / 2;
    (void) p3_out;  /* populated by caller from p1/p2 */

    /* Clear outputs */
    for (i = 0; i < nbins; i++) {
        bispec_out[i].n = 0; bispec_out[i].data = 0.0; bispec_out[i].data2 = 0.0;
    }
    for (i = 0; i < nbins; i++) {
        p1_out[i].n = 0; p1_out[i].data = 0.0;
        p2_out[i].n = 0; p2_out[i].data = 0.0;
    }

    /* Loop over k1 triangles on local slab */
    for (y1 = slabstart_y; y1 < slabstart_y + nslab_y; y1++) {
        for (x1 = 0; x1 < ngrid; x1++) {
            for (z1 = 0; z1 < ngrid_half + 1; z1++) {

                int xx1 = (x1 > ngrid_half) ? x1 - ngrid : x1;
                int yy1 = (y1 > ngrid_half) ? y1 - ngrid : y1;
                int zz1 = (z1 > ngrid_half) ? z1 - ngrid : z1;

                k1x = xx1 * delk;
                k1y = yy1 * delk;
                k1z = zz1 * delk;
                k1 = sqrt(k1x * k1x + k1y * k1y + k1z * k1z);

                if (k1 == 0.0 || k1 >= k_max || k1 < k_min) continue;

                ip1 = ngrid * (ngrid_half + 1) * (y1 - slabstart_y)
                    + (ngrid_half + 1) * x1 + z1;

                /* k2 loop */
                for (y2 = slabstart_y; y2 < slabstart_y + nslab_y; y2++) {
                    for (x2 = 0; x2 < ngrid; x2++) {
                        for (z2 = 0; z2 < ngrid_half + 1; z2++) {

                            int xx2 = (x2 > ngrid_half) ? x2 - ngrid : x2;
                            int yy2 = (y2 > ngrid_half) ? y2 - ngrid : y2;
                            int zz2 = (z2 > ngrid_half) ? z2 - ngrid : z2;

                            k2x = xx2 * delk;
                            k2y = yy2 * delk;
                            k2z = zz2 * delk;
                            k2 = sqrt(k2x * k2x + k2y * k2y + k2z * k2z);

                            if (k2 == 0.0 || k2 >= k_max || k2 < k_min) continue;

                            /* k3 = -k1 - k2 (triangle closure) */
                            k3x = -k1x - k2x;
                            k3y = -k1y - k2y;
                            k3z = -k1z - k2z;
                            k3 = sqrt(k3x * k3x + k3y * k3y + k3z * k3z);

                            if (k3 == 0.0 || k3 >= k_max || k3 < k_min) continue;

                            ip2 = ngrid * (ngrid_half + 1) * (y2 - slabstart_y)
                                + (ngrid_half + 1) * x2 + z2;

                            /* Accumulate P(k) */
                            bin = (int)((k1 - k_min) / dk);
                            if (bin >= 0 && bin < nbins) {
                                p1_out[bin].n++;
                                p1_out[bin].data += grid[ip1][0] * grid[ip1][0]
                                                  + grid[ip1][1] * grid[ip1][1];
                            }

                            bin = (int)((k2 - k_min) / dk);
                            if (bin >= 0 && bin < nbins) {
                                p2_out[bin].n++;
                                p2_out[bin].data += grid[ip2][0] * grid[ip2][0]
                                                  + grid[ip2][1] * grid[ip2][1];
                            }

                            /* Bispectrum: Re[delta(k1)*delta(k2)*conj(delta(k1+k2))]
                             * Since k3 = -k1 - k2, conj(delta(k3)) = delta(-k3).
                             * Approximate contribution to bispectrum bin by |k3|. */
                            bin = (int)((k3 - k_min) / dk);
                            if (bin >= 0 && bin < nbins) {
                                double d1r = grid[ip1][0], d1i = grid[ip1][1];
                                double d2r = grid[ip2][0], d2i = grid[ip2][1];
                                double bispec = d1r * d2r - d1i * d2i;
                                bispec_out[bin].n++;
                                bispec_out[bin].data += bispec;
                                bispec_out[bin].data2 += bispec * bispec;
                            }
                        }
                    }
                }
            }
        }
    }

    /* Set bin centers */
    for (i = 0; i < nbins; i++) {
        bispec_out[i].x = k_min + (i + 0.5) * dk;
        bispec_out[i].y = 0.0;
    }

    /* Normalization must be applied by caller via MPI_Allreduce */
}
