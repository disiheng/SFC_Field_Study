#ifndef ALLVARS_H
#define ALLVARS_H

#include <stdio.h>
#include <stdlib.h>

#ifdef NOTYPEPREFIX_FFTW
#include      <rfftw_mpi.h>
#else
#ifdef DOUBLEPRECISION_FFTW
#include     <drfftw_mpi.h> /* double precision FFTW */
#define SEND_TYPE MPI_DOUBLE
#else
#include     <srfftw_mpi.h>
#define SEND_TYPE MPI_FLOAT
#endif
#endif

#define MAX(a,b) ((a>b)?a:b)
#define MIN(a,b) ((a<b)?a:b)
#define SQR(a) (a*a)
#define CUB(a) (a*a*a)


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

extern rfftwnd_mpi_plan fft_forward_plan, fft_inverse_plan;
extern int *first_slab_of_task, *slabs_per_task, *slab_to_task;
extern int *first_slab_file, *nslabs_file;
extern fftw_complex *fft_Grid, *fft_sGrid, *fft_Vel;
extern fftw_real *Grid, *sGrid, *Vel;
extern fftw_real *workspace;
extern peanokey *Keys;
extern int slabstart_x, nslab_x, slabstart_y, nslab_y, smallest_slab;

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
  int npart[6];         /*!< number of particles of each type in this file */
  double mass[6];       /*!< mass of particles of each type. If 0, then the masses are explicitly
                   stored in the mass-block of the snapshot file, otherwise they are omitted */
  double time;          /*!< time of snapshot file */
  double redshift;      /*!< redshift of snapshot file */
  int flag_sfr;         /*!< flags whether the simulation was including star formation */
  int flag_feedback;        /*!< flags whether feedback was included (obsolete) */
  unsigned int npartTotal[6];   /*!< total number of particles of each type in this snapshot. This can be
                   different from npart if one is dealing with a multi-file snapshot. */
  int flag_cooling;     /*!< flags whether cooling was included  */
  int num_files;        /*!< number of files in multi-file snapshot */
  double BoxSize;       /*!< box-size of simulation in case periodic boundaries were used */
  double Omega0;        /*!< matter density in units of critical density */
  double OmegaLambda;       /*!< cosmological constant parameter */
  double HubbleParam;       /*!< Hubble parameter in units of 100 km/sec/Mpc */
  int flag_stellarage;      /*!< flags whether the file contains formation times of star particles */
  int flag_metals;      /*!< flags whether the file contains metallicity values for gas and star
                   particles */
  unsigned int npartTotalHighWord[6];   /*!< High word of the total number of particles of each type */
  int flag_entropy_instead_u;   /*!< flags that IC-file contains entropy instead of u */
  int flag_doubleprecision; /*!< flags that snapshot contains double-precision instead of single precision */
  
  int flag_ic_info;             /*!< flag to inform whether IC files are generated with ordinary Zeldovich approximation,
                                     or whether they ocontains 2nd order lagrangian perturbation theory initial conditions.
                                     For snapshots files, the value informs whether the simulation was evolved from
                                     Zeldoch or 2lpt ICs. Encoding is as follows:
                                        FLAG_ZELDOVICH_ICS     (1)   - IC file based on Zeldovich
                                        FLAG_SECOND_ORDER_ICS  (2)   - Special IC-file containing 2lpt masses
                                        FLAG_EVOLVED_ZELDOVICH (3)   - snapshot evolved from Zeldovich ICs
                                        FLAG_EVOLVED_2LPT      (4)   - snapshot evolved from 2lpt ICs
                                        FLAG_NORMALICS_2LPT    (5)   - standard gadget file format with 2lpt ICs
                                     All other values, including 0 are interpreted as "don't know" for backwards compatability.
                                 */
  float lpt_scalingfactor;      /*!< scaling factor for 2lpt initial conditions */
  char fill[48];        /*!< fills to 256 Bytes */
}
io_header;             /*!< holds header for snapshot files */


extern Particle *PP;
extern double BoxSize, Smoothing; 
extern double Time, Omega0; 
extern float to_slab_fac;
extern int *first_slab_of_task, *slabs_per_task;
extern long long Ndim, NFiles;
extern int Dim;
extern int NTask, ThisTask, PTask, rep;

extern int fftsize, maxfftsize;

extern int SnapshotNum;
extern long long TotNumHalo;
extern int LocNgroups;

extern char OutputDir[512];
extern char SimulationDir[512];
extern char FileNameBase[512];

extern int  LastSnapShotNr;
extern int  MaxMemSize;


#define  terminate(x) {char __buf[1000]; sprintf(__buf, "code termination on task=%d, function '%s()', file '%s', line %d: '%s'\n", ThisTask, __FUNCTION__, __FILE__, __LINE__, x); printf(__buf); fflush(stdout);                           MPI_Abort(MPI_COMM_WORLD, 1); exit(0);}

#include "mymalloc.h"


#endif



