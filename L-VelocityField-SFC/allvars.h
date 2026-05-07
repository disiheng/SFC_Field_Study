#ifndef ALLVARS_H
#define ALLVARS_H

#include <stdio.h>
#include <stdlib.h>
#include <fftw3-mpi.h>

#define MAX(a,b) ((a>b)?a:b)
#define MIN(a,b) ((a<b)?a:b)
#define SQR(a) (a*a)
#define CUB(a) (a*a*a)

/* FFTW3 uses double for real data directly */
typedef double fftw_real;
#define SEND_TYPE MPI_DOUBLE

typedef long long large_array_offset;
typedef unsigned long long peanokey;

#define TAG_N          10
#define TAG_PDATA      11
#define TAG_MDATA      12
#define TAG_POSDATA    13

#ifndef LONGIDS
typedef unsigned int MyIDType;
#else
typedef unsigned long long MyIDType;
#endif


extern time_t t0;

/* FFTW3 plans */
extern fftw_plan fft_forward_plan, fft_inverse_plan;
extern int *first_slab_of_task, *slabs_per_task, *slab_to_task;
extern int *first_slab_file, *nslabs_file;
extern fftw_complex *fft_Grid, *fft_sGrid, *fft_Vel;
extern fftw_real *Grid, *sGrid, *Vel;
extern fftw_real *workspace;
extern peanokey *Keys;
extern int slabstart_x, nslab_x, slabstart_y, nslab_y, smallest_slab;

/* FFTW3 transposed local sizes */
extern ptrdiff_t local_n0, local_0_start, local_n1, local_1_start;

extern size_t HighMark_run;

extern int NumFilesWrittenInParallel;
extern int CYCLES;

#ifndef PARTICLE_TYPEDEF_DONE
typedef struct Particle_ {
  float pos[3];
  float vel[3];
} Particle;
#define PARTICLE_TYPEDEF_DONE
#endif

typedef struct
{
  int npart[6];
  double mass[6];
  double time;
  double redshift;
  int flag_sfr;
  int flag_feedback;
  unsigned int npartTotal[6];
  int flag_cooling;
  int num_files;
  double BoxSize;
  double Omega0;
  double OmegaLambda;
  double HubbleParam;
  int flag_stellarage;
  int flag_metals;
  unsigned int npartTotalHighWord[6];
  int flag_entropy_instead_u;
  int flag_doubleprecision;
  int flag_ic_info;
  float lpt_scalingfactor;
  char fill[48];
}
io_header;


extern Particle *PP;
extern double BoxSize, Smoothing;
extern double Time, Omega0;
extern float to_slab_fac;
extern int *first_slab_of_task, *slabs_per_task;
extern long long Ndim, NFiles;
extern int Dim;
extern int NTask, ThisTask, PTask, rep;

extern ptrdiff_t fftsize, maxfftsize;

extern int SnapshotNum;
extern long long TotNumHalo;
extern int LocNgroups;

extern char OutputDir[512];
extern char SimulationDir[512];
extern char FileNameBase[512];
extern char Mode[64];

extern int  LastSnapShotNr;
extern int  MaxMemSize;


#define  terminate(x) {char __buf[1000]; sprintf(__buf, "code termination on task=%d, function '%s()', file '%s', line %d: '%s'\n", ThisTask, __FUNCTION__, __FILE__, __LINE__, x); printf(__buf); fflush(stdout);                           MPI_Abort(MPI_COMM_WORLD, 1); exit(0);}

#include "mymalloc.h"


#endif
