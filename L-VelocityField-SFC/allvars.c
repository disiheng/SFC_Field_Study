
#include "allvars.h"

time_t t0;

char OutputDir[512];
char SimulationDir[512];
char FileNameBase[512];

int LastSnapShotNr, Dim;

long long NFiles, Ndim;
int  bits_per_dimension;

float to_slab_fac;
double BoxSize, Smoothing;
double Time, Omega0;

int NTask, ThisTask, PTask, rep;

int fftsize, maxfftsize;

int NumFilesWrittenInParallel;
int CYCLES;

int SnapshotNum;
long long TotNumHalo;
int LocNgroups;

int LastSnapShotNr;
int MaxMemSize;

size_t HighMark_run;

size_t AllocatedBytes;
size_t HighMarkBytes;
float HighMarkBytes_Step;
size_t FreeBytes;

rfftwnd_mpi_plan fft_forward_plan, fft_inverse_plan;
int *first_slab_of_task, *slabs_per_task, *slab_to_task;
int *first_slab_file, *nslabs_file;
fftw_complex *fft_Grid, *fft_sGrid, *fft_Vel;
fftw_real *Grid, *sGrid, *Vel;
fftw_real *workspace;

int slabstart_x, nslab_x, slabstart_y, nslab_y, smallest_slab;

Particle *PP;

