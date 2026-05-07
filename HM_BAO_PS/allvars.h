
#ifdef DOUBLEPRECISION_FFTW
#include     <drfftw_mpi.h>
#else
#include     <srfftw_mpi.h>
#endif


#include "io.h"
#include "read_data.h"

#define TAG_N          10
#define TAG_PDATA      11
#define CYCLES 4
#define FOLD_FACTOR 8 //16
#define FOLDS 0

/*
#ifdef DM
#define NBINS 1  
#define SEL_DM
#else
#define NBINS 6  
//#define SEL_MASS
//#define SEL_SFR
#define CENTRAL_HALO
#endif
*/

#define NMAX     6000000
#define NMAXRDM  100000
#define NCELLMAX 10000
#define NMAXCHAR 255
#define NBINMAX  100
#define Pi       3.14159
#define D2R (acos(-1.0)/180.0)
#define MIN(a,b) ((a)<(b)?(a):(b))
#define MAX(a,b) ((a)<(b)?(b):(a))
#define SQR(a) ((a)*(a))
#define CUB(a) ((a)*(a)*(a))

#ifndef DEBUG
  #define DEBUG 0
#endif

typedef struct result_data
{
  long long n; 
  float x; 
  float y;
  float r;
  float data;
  float data2;
} result;

extern int ThisTask, PTask, NTasks, NProcRead, NCell;
extern float MinMag, MaxMag;
extern int  NRegrid, NFilesDM, NFilesGr, NGens, NNear,NBins, NTimesHal;
extern int  Geometry, AssgnScheme, npow, SimFormat;
extern int  flag_write_df, flag_redshift,flag_hdf5, flag_log, flag_periodic, flag_fits, flag_mpi;
extern int iseed, Ngrid2, Num, NumFilesWrittenInParallel;
extern double Hubble, ExpFactor,ToMpc,Center, kTMin, ML1,ML2;
extern double CellSize;
extern char  OutFile[NMAXCHAR], SimFile[NMAXCHAR], SnapBase[NMAXCHAR], OutFile[NMAXCHAR];
extern double MinConc, MaxConc,Redshift, Omega_m, DiluteFactor, minScale, maxScale;
typedef int int4byte;
extern double *MyBoxSize; //, BoxSize;
extern int *Send_offset, *Send_count, *Recv_count, *Recv_offset, *Sendcount;

extern int local_x_start, local_nx, local_y_start, local_ny, total_local_size, max_size;
extern int SnapNum, FileFormat, NFiles;
extern int *slab_to_task, *slab_to_task_local, *slabs_per_task, *first_slab_of_task;
extern float N1, N2;
extern double BoxSize, invCellSize, to_slab_fac;
extern double rms, mean, sd, skew, kurtosis;
extern char DMFile[255], HalFile[255];
extern time_t t0;
extern char OutputDir[255], SimuTag[255], SnapshotFileBase[255], SimulationDir[255], CrossOutputDir[255];
extern unsigned int LocNgroups;
extern Particle *PP;

extern int *CataFileMap, *IterstoLoad, *FileNumTab;
extern int SubVolDim, CataVolDim;
extern int NumIterstoLoad, iDiv;
extern double DomainCenter[3], DomainDepth[3];

typedef unsigned long long peanokey;
#define  ITERVOL_BITS_PER_DIMENSION 1 //  2
#define  CATAVOL_BITS_PER_DIMENSION 3 //  8


/*
struct data_struct
{
  int   Type; // Galaxy type: 0 for central galaxies of a main halo, 1 for central galaxies in sub-halos, 2 for satellites without halo.
  int   SnapNum; // The snapshot number where this galaxy was identified.
  float CentralMvir; // 10^10/h Msun virial mass of background (FOF) halo containing this galaxy
  float Pos[3]; // 1/h Mpc - Galaxy Positions
  float BulgeMass; // 10^10/h Msun - Mass in the bulge
  float DiskMass;
  float Distance;
  float Time;
  float Redshift;
}*halo_in;
*/

//   float tau[256][5];
