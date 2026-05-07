#ifndef SFC_H
#define SFC_H

#include <stdint.h>
#include "particle_io.h"

/* Space-Filling Curve ordering: maps particles to MRA grid cells.
 * grid_len = 2^J (side length of MRA grid)
 * Returns array of grid indices [0, grid_len^3) for each particle.
 * Caller must free the returned array.
 */
uint64_t *sfc_order(const Particle *particles, long long n, uint64_t grid_len);

/* Same as sfc_order but with an offset applied before SFC lookup */
uint64_t *sfc_order_offset(const Particle *particles, long long n,
                            uint64_t grid_len, double dx, double dy, double dz,
                            double boxsize);

#endif
