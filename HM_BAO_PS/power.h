
#ifndef POWER_H
#define POWER_H

#define POWERSPEC_FOLDFAC 32 

#include "density.h"

void power( char fname[], Grid grid, Grid grid_x, int flag, float box); 
void divergence(Grid grid, Grid grid_vx, Grid grid_vy, Grid grid_vz, Grid * theta_grid, float boxsize);

double deconvolution_factor(double kx, double ky, double kz, int ngrid);
double modulus(int x, int ngrid);
int inv_modulus(double kx, int ngrid);
int fftw_index(double kx, double ky, double kz, int ngrid, int local_y_start, int *sign);

void gather_results(result *in, int size, double norm);
void write_results(result *in, int sizex, int sizey, int type, char buf[]);

#endif
