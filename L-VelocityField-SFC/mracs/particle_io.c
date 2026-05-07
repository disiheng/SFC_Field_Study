#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "particle_io.h"

long long read_millennium_galaxies(const char *path, Particle **particles, int use_mass_weight)
{
    FILE *fp;
    long long ngal, i;
    float *buf;

    fp = fopen(path, "rb");
    if (!fp) {
        fprintf(stderr, "Error: cannot open %s\n", path);
        return -1;
    }

    /* Determine file size */
    fseek(fp, 0, SEEK_END);
    long long file_size = ftell(fp);
    rewind(fp);

    /* Each galaxy = 23 floats = 92 bytes */
    ngal = file_size / (23 * sizeof(float));

    *particles = (Particle *) malloc(ngal * sizeof(Particle));
    if (!*particles) {
        fprintf(stderr, "Error: malloc failed for %lld particles\n", ngal);
        fclose(fp);
        return -1;
    }

    buf = (float *) malloc(23 * sizeof(float));

    for (i = 0; i < ngal; i++) {
        if (fread(buf, sizeof(float), 23, fp) != 23) break;

        (*particles)[i].x = (double) buf[0];
        (*particles)[i].y = (double) buf[1];
        (*particles)[i].z = (double) buf[2];

        if (use_mass_weight) {
            (*particles)[i].weight = (double) buf[16];  /* stellar_mass */
        } else {
            (*particles)[i].weight = 1.0;
        }
    }

    free(buf);
    fclose(fp);
    return i;
}

long long read_tng_particles(const char *path, Particle **particles)
{
    FILE *fp;
    long long npart, i;
    double pos[3];

    fp = fopen(path, "rb");
    if (!fp) {
        fprintf(stderr, "Error: cannot open %s\n", path);
        return -1;
    }

    fseek(fp, 0, SEEK_END);
    long long file_size = ftell(fp);
    rewind(fp);

    /* Each particle = 3 doubles (x,y,z) */
    npart = file_size / (3 * sizeof(double));

    *particles = (Particle *) malloc(npart * sizeof(Particle));
    if (!*particles) {
        fprintf(stderr, "Error: malloc failed for %lld particles\n", npart);
        fclose(fp);
        return -1;
    }

    for (i = 0; i < npart; i++) {
        if (fread(pos, sizeof(double), 3, fp) != 3) break;

        (*particles)[i].x = pos[0];
        (*particles)[i].y = pos[1];
        (*particles)[i].z = pos[2];
        (*particles)[i].weight = 1.0;
    }

    fclose(fp);
    return i;
}

void free_particles(Particle *p)
{
    free(p);
}
