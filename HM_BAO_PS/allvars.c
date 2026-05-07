#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <assert.h>
#include <mpi.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "allvars.h"
#include "read_data.h"

int ThisTask, PTask, NTasks, NProcRead, NCell;
float MinMag, MaxMag;
int  NRegrid, NFilesDM, NFilesGr, NGens, NNear,NBins, NTimesHal;
int  Geometry, AssgnScheme, npow, SimFormat;
int  flag_write_df, flag_redshift,flag_hdf5, flag_log, flag_periodic, flag_fits, flag_mpi;
int iseed, Ngrid2, Num, NumFilesWrittenInParallel;
double Hubble, ExpFactor,ToMpc,Center, kTMin, ML1,ML2;
double CellSize;
char  OutFile[NMAXCHAR], SimFile[NMAXCHAR], SnapBase[NMAXCHAR], OutFile[NMAXCHAR];
double MinConc, MaxConc,Redshift, Omega_m, DiluteFactor, minScale, maxScale;
double *MyBoxSize; //, BoxSize;
int *Send_offset, *Send_count, *Recv_count, *Recv_offset, *Sendcount;

int local_x_start, local_nx, local_y_start, local_ny, total_local_size, max_size;
int SnapNum, FileFormat, NFiles;
int *slab_to_task, *slab_to_task_local, *slabs_per_task, *first_slab_of_task;
float N1, N2;
double BoxSize, invCellSize, to_slab_fac;
double rms, mean, sd, skew, kurtosis;
char DMFile[255], HalFile[255];
char OutputDir[255], SimuTag[255], SnapshotFileBase[255], SimulationDir[255], CrossOutputDir[255];
unsigned int LocNgroups;
time_t t0;
Particle *PP;

int *CataFileMap, *IterstoLoad, *FileNumTab;
int SubVolDim, CataVolDim;
int NumIterstoLoad, iDiv;
double DomainCenter[3], DomainDepth[3];


