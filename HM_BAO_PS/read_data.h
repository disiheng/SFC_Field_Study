
#ifndef READ_H
#define READ_H

#include <stdlib.h>

#ifdef GADGET2 
typedef struct
{
  int npart[6];
  double mass[6];
  double time;
  double redshift;
  int flag_sfr;
  int flag_feedback;
  int npartTotal[6];
  int flag_cooling;
  int num_files;
  double BoxSize;
  double Omega0;
  double OmegaLambda;
  double HubbleParam;
  int flag_stellarage;
  int flag_metals;
  int hashtabsize;
  char fill[84];                /* fills to 256 Bytes */
}
io_header;

#else 
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
  char fill[48];        /*!< fills to 256 Bytes */
}
io_header;
#endif

typedef struct Particle_ {
  float pos[3];
  float mass;

#if defined(DIVERGENCE) || defined(VELOCITY)
  float vel[3];
#endif

#ifdef ID
  unsigned long long ID;
#endif

#ifdef GAL_SHUFFLE 
  int Type;
  float HaloMass;
#endif

#ifdef ONEHALOTERM
  long long fofid;
#endif
} Particle;

typedef struct Halo_ {
  float pos[3];
  float vel[3];
  float mass;
  unsigned long long ID;
} Halo;


Particle *read_particles (char *File, int NFiles, int Num, long *NPart, 
		void (*readfile)(char *filename, Particle *P, long *NPart, float DiluteFactor),
		long (*get_number_particles)(char *, int Num)
		 );

void read_fof_file(char *filename, Particle * P, long *NPart, float DiluteFactor);
long get_number_of_fof(char *filename);
long get_number_of_haloes(char *file, int Num);
void read_gadget_file(char *filename, Particle *P, long *NPart, float DiluteFactor);
void read_halo_file( char *filename, Halo *P, long *NPart, float DiluteFactor);
void read_halo_file_lg3( char *filename, Particle *P, long *NPart, float DiluteFactor);

int find_block(FILE *fd,char *label);
unsigned long long *read_id_particles_in_groups(char *OutputDir, int Num, int iFile, int *Nid, int *GroupOff, int *GroupLen);
int get_number_of_dm_files(char *file, int Num);
int get_number_of_halo_files(char *file, int Num);

void read_scaled_halo_file(char *filename, Particle * P, long *NPart, float diluteFactor);
long get_number_of_scaled_halos(char *filename);
int blksize;

#define SKIP  {my_fread(&blksize,sizeof(int),1,fd);}
//#define FORMAT2
#define NTASKREAD 8
#define DILUTE 1.0

int PHOTOZFILE;
int flag_photoz;
int flag_redshift;
int flag_ids;
int flag_distribute;
int flag_vel;
float photoz;
char photoz_file[255];

int MinMass, MaxMass;

#endif

