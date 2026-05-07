#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#include "sfc.h"

/* Morton Z-order: interleave bits of x, y, z into a 64-bit key */
static uint64_t morton_key(uint64_t x, uint64_t y, uint64_t z)
{
    uint64_t key = 0;
    int i;
    /* 21 bits per dimension fits in 64-bit (63 bits total) */
    for (i = 0; i < 21; i++) {
        key |= ((x >> i) & 1ULL) << (3 * i);
        key |= ((y >> i) & 1ULL) << (3 * i + 1);
        key |= ((z >> i) & 1ULL) << (3 * i + 2);
    }
    return key;
}

uint64_t *sfc_order(const Particle *particles, long long n, uint64_t grid_len)
{
    return sfc_order_offset(particles, n, grid_len, 0.0, 0.0, 0.0, 1.0);
}

uint64_t *sfc_order_offset(const Particle *particles, long long n,
                            uint64_t grid_len, double dx, double dy, double dz,
                            double boxsize)
{
    uint64_t *indices;
    long long i;
    uint64_t ix, iy, iz;
    double inv_cell;

    if (n <= 0 || !particles) return NULL;

    indices = (uint64_t *) malloc(n * sizeof(uint64_t));
    if (!indices) return NULL;

    inv_cell = (double) grid_len / boxsize;

    for (i = 0; i < n; i++) {
        double x = (particles[i].x + dx) * inv_cell;
        double y = (particles[i].y + dy) * inv_cell;
        double z = (particles[i].z + dz) * inv_cell;

        /* Periodic wrap */
        while (x < 0.0) x += (double) grid_len;
        while (y < 0.0) y += (double) grid_len;
        while (z < 0.0) z += (double) grid_len;

        ix = (uint64_t) x % grid_len;
        iy = (uint64_t) y % grid_len;
        iz = (uint64_t) z % grid_len;

        indices[i] = morton_key(ix, iy, iz);
    }

    return indices;
}
