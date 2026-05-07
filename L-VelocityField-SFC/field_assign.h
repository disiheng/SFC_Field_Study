#ifndef FIELD_ASSIGN_H
#define FIELD_ASSIGN_H

#include "allvars.h"

void density(Particle *P, fftw_real *grid, long long NumPart, int mode, int value);
void mapping(Particle *P, fftw_real *grid, long long FirstPart, long long NumPart, int mode, int value);

#endif
