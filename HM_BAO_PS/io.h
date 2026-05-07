

#ifndef IO_H
#define IO_H

#include "density.h"
struct data_struct
{
   int Type; 
   int SnapNum; 
   float CentralMvir;
   float Pos[3];
   float Vel[3];
   float Mvir; 
   float Rvir;
   float Vvir; 
   float DistanceToCentralGal[3];
   float ColdGas;
   float BulgeMass;
   float DiskMass;
   float HotGas;
   float BlackHoleMass; 
   float SFR; 
}
 *halo_in;

struct gal_data_struct
{
   int Type; 
   int SnapNum; 
   float CentralMvir;
   float Pos[3];
   float Vel[3];
   float Mvir; 
   float Rvir;
   float Vvir; 
   float DistanceToCentralGal[3];
   float ColdGas;
   float BulgeMass;
   float DiskMass;
   float HotGas;
   float BlackHoleMass; 
   float SFR; 
}
 *gal_in;

//struct simple_data_struct
//{
//   float Pos[3];
//   float Vel[3];
//   int Len; 
//   float Vmax; 
//   float VmaxRad; 
//}
// *data_in;

struct simple_data_struct
{
  float Pos[3];
  float Vel;
  float Zred; 
  float Mvir; 
}
 *data_in;

struct simple_data_struct2
{
  int FileIndex;
  int TreeIndex;
  int HaloIndex;
  int SnapNum;
  int SnapNumatInFall;
  int Len;
  int Type;
  float M_Crit200_atInFall;
  float M_Mean200_atInFall;
  float M_TopHat_atInFall;
  float VmaxatInFall;
  float Pos[3];
}
 *data_in2;


void load_catalogue();
void load_bruno_single( int nr);
void load_dm_single( int nr, float nsel);
void load_catalogue_single( int nr);
void dump_particles(int np);

double shuffle_periodic(double x);

void get_slabs_file();
void load_field(Grid *grid);

int *first_slab_file, *nslabs_file;

void load_dm_single_hdf5(int nr, float nsel);

void load_file_tab(void);
void load_halo_catalogue_single_photons(int nr);
int ShouldWeLoadIter(float x, float y, float z, float IterHalfSize);
int ShouldWeLoadHalo(float x, float y, float z);
void mark_domain_center(int iDiv);
void mark_Iter_nums(void);

void load_halo_catalogue_single(int nr);

void load_simple_catalogue_single( int nr);
void load_simple_catalogue_single2( int nr);

double periodic_wrap(double x);
#endif

