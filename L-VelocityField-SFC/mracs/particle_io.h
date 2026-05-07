#ifndef PARTICLE_IO_H
#define PARTICLE_IO_H

#include <stdint.h>

/* Particle with position and weight */
typedef struct {
    double x, y, z;
    double weight;
} Particle;

/* Millennium Run galaxy (23 float fields in binary) */
typedef struct {
    float x, y, z;
    float vx, vy, vz;
    float mag_u, mag_g, mag_r, mag_i, mag_z;
    float bulge_mag_u, bulge_mag_g, bulge_mag_r, bulge_mag_i, bulge_mag_z;
    float stellar_mass, bulge_mass, cold_gas, hot_gas, ejected_mass, black_hole_mass, sfr;
} Galaxy;

long long read_millennium_galaxies(const char *path, Particle **particles, int use_mass_weight);
long long read_tng_particles(const char *path, Particle **particles);
void free_particles(Particle *p);

#endif
