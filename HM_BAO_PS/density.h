
#ifndef DENSITY_H
#define DENSITY_H

#include <stdlib.h>
#include "read_data.h"

#define TAG_N          10
#define TAG_PDATA      11
#define TAG_MDATA      12
#define TAG_POSDATA    13

#ifdef DOUBLEPRECISION
typedef double  GridFloat;
typedef double  MyFloat;
#define SEND_TYPE MPI_DOUBLE
#else
typedef float  GridFloat;
typedef float  MyFloat;
#define SEND_TYPE MPI_FLOAT
#endif

typedef struct Grid_ {

    GridFloat *data;
#ifdef CENTER_OF_MASS
    GridFloat *pos;
#endif
    long long Ngrid;
    int local_nx;
    int local_x_start;
	double m1;
	double m2;
	double m3;
	double dm2;
	double dm3;
	double to_slab_fac;
	double inv_to_slab_fac;
	int total_local_size;
	long total;
    int *slab_to_task;
    int *slabs_per_task;
    int *first_slab_of_task;
    rfftwnd_mpi_plan fft_forward_plan;
    rfftwnd_mpi_plan fft_inverse_plan;
    int local_ny;
    int local_y_start; 
} Grid;

void smooth(Grid *orig, Grid *smooth, float SmoothingScale);
Grid init_density (int Ngrid, int foldnum);
void density(Particle *P, long long NumPart, Grid *grid, int nfold, int value);
void grid_clean(Grid grid);
void grid_init(Grid *grid, int foldnum);
void mapping( Particle *P, long FirstPart, long NumPart, Grid *grid, int nfold, int mode, int value);
double density_statistics (Grid *grid, const char *label, int foldnum);
void density_save (char fname[], Grid grid);
Particle *create_density (Grid N_grid, float *bias, float N2, Grid DM_grid, long *N_NPart);

#endif
