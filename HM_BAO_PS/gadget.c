
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "gadget.h"
#include "byte_swap.h"
#include "utils.h"


#ifndef PARTICLE_TYPE
#define PARTICLE_TYPE 1+2+4+8+16+32 
#endif

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
  char fill[84];		/* fills to 256 Bytes */
}
io_header;

typedef struct
{
  int len;
  int file;
  int offset;
}
cat_data;

static double BoxSize;

int millenniumdir(int Num, int readFile)
{
  FILE *fd;
  char buf[1000];
  int i, DirNum = -1;

  while(DirNum == -1)
    {
      for(i = 0; i < 14; i++)
	{
	  sprintf(buf, "/data/milli%02d/d%02d/snapshot/snap_millennium_%03d.%d", i, Num, Num, readFile);
	  if((fd = fopen(buf, "r")))
	    {
	      DirNum = i;
	      fclose(fd);
	    }
	}
      if(DirNum == -1)
	printf("couldn't find any directory %s\n", buf);
    }
  return DirNum;
}

int hrdir(int Num, int readFile)
{
  FILE *fd;
  char buf[1000];
  int i, DirNum = -1;

  while(DirNum == -1)
    {
      for(i = 0; i < 14; i++)
	{
	  sprintf(buf, "/data/milli%02d/hr1%02d/snapshot/hr1_snap_%03d.%d", i, Num, Num, readFile);
	  if((fd = fopen(buf, "r")))
	    {
	      DirNum = i;
	      fclose(fd);
	    }
	}
      if(DirNum == -1)
	{
	  printf("hrdir: couldn't find any directory for %d in Num %d  \n", readFile, Num);
	  printf("hrdir: %s  \n", buf);
	  exit(0);
	}
    }
  return DirNum;
}

int get_subhalo_file(long *FileMap, int Num, int NFiles, int type)
{
  FILE *fd;
  int i, dummy, DirNum;
  int Ngroups, Nids, TotNgroups, Nfiles, Nsubs;
  long TotNsubs;
  char buf[1000];

  // sprintf(buf, "%s/snapdir_%03d/%s_%03d.0", OutputDir,Num,SnapBase, Num);
  TotNsubs = 0;
  for(i = 0; i < NFiles; i++)
    {
      if(type <= 2)
	{
	  DirNum = millenniumdir(Num, i);
	  sprintf(buf, "/data/milli%02d/d%02d/sub_tab/sub_tab_%03d.%d", DirNum, Num, Num, i);
	}
      if(type >= 3 && type <= 5)
	{
	  DirNum = hrdir(Num, i);
	  sprintf(buf, "/data/milli%02d/hr1%02d/sub_tab/sub_tab_%03d.%d", DirNum, Num, Num, i);
	}

      if(type >= 6 && type <= 6)
	sprintf(buf, "/data/rw15/reangulo/PS-Sim/ps1/Output/postproc_%03d/sub_tab_%03d.%d", Num, Num, i);

      if(type == 7)
	sprintf(buf, "/data/rw15/reangulo/hr1_time_test/postproc_%03d/sub_tab_%03d.%d", Num, Num, i);

      if(type == 8)
	sprintf(buf, "/data/rw15/reangulo/mill_time_test/postproc_%03d/sub_tab_%03d.%d", Num, Num, i);

      fd = my_fopen(buf, "r");

      my_fread(&Ngroups, sizeof(int), 1, fd);
      my_fread(&Nids, sizeof(int), 1, fd);
      my_fread(&TotNgroups, sizeof(int), 1, fd);
      my_fread(&Nfiles, sizeof(int), 1, fd);
      my_fread(&Nsubs, sizeof(int), 1, fd);

      fclose(fd);

#ifdef SWAPENDIAN
      byte_swap_int(1, &Ngroups);
      byte_swap_int(1, &Nids);
      byte_swap_int(1, &TotNgroups);
      byte_swap_int(1, &Nfiles);
      byte_swap_int(1, &Nsubs);
#endif

      //printf("%d %d %d %d %d\n",Ngroups, Nids, TotNgroups, Nfiles, Nsubs);
      //TotNsubs += Nsubs;
      TotNsubs += Ngroups;
      FileMap[i] = TotNsubs;
    }
  return 0;
}


float get_hashcell_size(char *OutputDir, char *SnapBase, int Num, int HashTabSize, int *base, int *hashbits)
{
  FILE *fd;
  int hashbit, dummy, bas, DirNum;
  char buf[1000];
  io_header header;

  sprintf(buf, "%s/%s_%03d.%d", OutputDir, SnapBase, Num, 0);

#ifdef MILLENNIUM
  DirNum = millenniumdir(Num, 0);
  sprintf(buf, "/data/milli%02d/d%02d/snapshot/snap_millennium_%03d.%d", DirNum, Num, Num, 0);
#endif

#ifdef HR1
  DirNum = hrdir(Num, 0);
  sprintf(buf, "/data/milli%02d/hr1%02d/snapshot/hr1_snap_%03d.%d", DirNum, Num, Num, 0);
#endif

  fd = my_fopen(buf, "r");

  my_fread(&dummy, sizeof(int), 1, fd);
  my_fread(&header, sizeof(io_header), 1, fd);
  my_fread(&dummy, sizeof(int), 1, fd);
  fclose(fd);

  bas = 1;
  while(bas * bas * bas < HashTabSize)
    bas *= 2;

  hashbit = 1;
  while((1 << hashbit) < bas)
    hashbit++;

  *base = bas;
  *hashbits = hashbit;

#ifdef SWAPENDIAN
  byte_swap_double(1, &(header.BoxSize));
#endif

  return header.BoxSize / bas;

}

int get_nfiles(char *filename)
{
  FILE *fd;
  int dummy, nfiles;
  char buf[1000];
  io_header header;

  fd = fopen(filename, "r");
  if(!fd)
    {
      sprintf(buf, "%s.0", filename);
      fd = my_fopen(buf, "r");
    }

  my_fread(&dummy, sizeof(int), 1, fd);
  my_fread(&header, sizeof(io_header), 1, fd);
  my_fread(&dummy, sizeof(int), 1, fd);
  fclose(fd);

#ifdef SWAPENDIAN
  byte_swap_int(1, &(header.num_files));
#endif
  nfiles = header.num_files;

  return nfiles;
}


float get_redshift(char *filename, int Num)
{
  FILE *fd;
  int dummy, DirNum;
  float redshift;
  char buf[1000];
  io_header header;

#ifdef MILLENNIUM
  DirNum = millenniumdir(Num, 0);
  sprintf(filename, "/data/milli%02d/d%02d/snapshot/snap_millennium_%03d.%d", DirNum, Num, Num, 0);
#endif

#ifdef HR1
  DirNum = hrdir(Num, 0);
  sprintf(filename, "/data/milli%02d/hr1%02d/snapshot/hr1_snap_%03d.%d", DirNum, Num, Num, 0);
#endif


  fd = fopen(filename, "r");
  if(!fd)
    {
      sprintf(buf, "%s.0", filename);
      fd = my_fopen(buf, "r");
    }

  my_fread(&dummy, sizeof(int), 1, fd);
  my_fread(&header, sizeof(io_header), 1, fd);
  my_fread(&dummy, sizeof(int), 1, fd);
  fclose(fd);

#ifdef SWAPENDIAN
  byte_swap_double(1, &(header.redshift));
#endif
  redshift = header.redshift;

  return redshift;

}

float get_mass(char *filename)
{
  FILE *fd;
  int dummy;
  float mass;
  char buf[1000];
  int DirNum, Num = 20;
  io_header header;

  sprintf(buf, "%s.0", filename);

#ifdef MILLENNIUM
  DirNum = millenniumdir(Num, 0);
  sprintf(buf, "/data/milli%02d/d%02d/snapshot/snap_millennium_%03d.%d", DirNum, Num, Num, 0);
#endif

#ifdef HR1
  DirNum = hrdir(Num, 0);
  sprintf(buf, "/data/milli%02d/hr1%02d/snapshot/hr1_snap_%03d.%d", DirNum, Num, Num, 0);
#endif

  fd = my_fopen(buf, "r");

  my_fread(&dummy, sizeof(int), 1, fd);
  my_fread(&header, sizeof(io_header), 1, fd);
  my_fread(&dummy, sizeof(int), 1, fd);
  fclose(fd);

#ifdef SWAPENDIAN
  byte_swap_double(1, &(header.mass[1]));
#endif
  mass = header.mass[1];

  return mass;

}

//Reetunrm the boxsize in a gadget file
double get_boxsize(char *filename, int Num)
{
  FILE *fd;
  int DirNum, dummy;
  long TotNumPart;
  io_header header;
  char buf[255];

#ifdef MILLENNIUM
  DirNum = millenniumdir(Num, 0);
  sprintf(filename, "/data/milli%02d/d%02d/snapshot/snap_millennium_%03d.%d", DirNum, Num, Num, 0);
#endif

#ifdef HR1
  DirNum = hrdir(Num, 0);
  sprintf(filename, "/data/milli%02d/hr1%02d/snapshot/hr1_snap_%03d.%d", DirNum, Num, Num, 0);
#endif


  fd = fopen(filename, "r");
  if(!fd)
    {
      sprintf(buf, "%s.0", filename);
      fd = my_fopen(buf, "r");
    }

#ifdef FORMAT2
  find_block(fd, "HEAD");
#endif

  my_fread(&dummy, sizeof(int), 1, fd);
  my_fread(&header, sizeof(io_header), 1, fd);
  my_fread(&dummy, sizeof(int), 1, fd);

#ifdef SWAPENDIAN
  byte_swap_double(1, &(header.BoxSize));
#endif

  BoxSize = header.BoxSize;
  fclose(fd);

  return BoxSize;
}


//Reetunrm the number of particles in a gadget file
long get_number_of_particles(char *filename, int Num)
{
  FILE *fd;
  int i, DirNum, dummy;
  long TotNumPart;
  io_header header;
  char buf[255];
/*
	
#ifdef MILLENNIUM
    DirNum = millenniumdir(Num,0);
    sprintf(filename, "/data/milli%02d/d%02d/snapshot/snap_millennium_%03d.%d",DirNum, Num, Num, 0);
#endif

#ifdef HR1
    DirNum = hrdir(Num,0);
    sprintf(filename, "/data/milli%02d/hr1%02d/snapshot/hr1_snap_%03d.%d",DirNum, Num, Num, 0);
#endif
*/

  fd = fopen(filename, "r");
  if(!fd)
    {
      sprintf(buf, "%s.0", filename);
      fd = my_fopen(buf, "r");
    }

#ifdef FORMAT2
  find_block(fd, "HEAD");
#endif
  my_fread(&dummy, sizeof(int), 1, fd);
  fread(&header, sizeof(io_header), 1, fd);
  fread(&dummy, sizeof(int), 1, fd);


#ifdef SWAPENDIAN
  for (i=0;i<6;i++)
    byte_swap_int(1, &(header.npart[i]));
#endif

  TotNumPart = 0 ;
  for (i=0;i<6;i++)
     if( ((1 << i) & (PARTICLE_TYPE)))
       TotNumPart += header.npart[i];
  fclose(fd);

  return TotNumPart;
}

long long get_total_number_of_particles(char *filename, int Num)
{
  FILE *fd;
  int i, DirNum, dummy;
  long long TotNumPart;
  char buf[255];
  io_header header;

#ifdef MILLENNIUM
  DirNum = millenniumdir(Num, 0);
  sprintf(filename, "/data/milli%02d/d%02d/snapshot/snap_millennium_%03d.%d", DirNum, Num, Num, 0);
#endif

#ifdef HR1
  DirNum = hrdir(Num, 0);
  sprintf(filename, "/data/milli%02d/hr1%02d/snapshot/hr1_snap_%03d.%d", DirNum, Num, Num, 0);
#endif

  fd = fopen(filename, "r");
  if(!fd)
    {
      sprintf(buf, "%s.0", filename);
      fd = my_fopen(buf, "r");
    }

  my_fread(&dummy, sizeof(int), 1, fd);
  my_fread(&header, sizeof(io_header), 1, fd);
  my_fread(&dummy, sizeof(int), 1, fd);

#ifdef SWAPENDIAN
  for (i=0;i<6;i++)
  byte_swap_int(1, &(header.npartTotal[i]));
#endif

  TotNumPart = 0 ;
  for (i=0;i<6;i++)
   if( ((1 << i) & (PARTICLE_TYPE)))
     TotNumPart += header.npartTotal[i];
  fclose(fd);

  return TotNumPart;
}

long get_number_of_groups(char *filename, int Num, int iFile)
{
  int DirNum, Ngroups, TotNgroups;
  char buf[1000];
  FILE *fd;

#ifdef MILLENNIUM
  DirNum = millenniumdir(Num, iFile);
  sprintf(buf, "/data/milli%02d/d%02d/snapshot/snap_millennium_%03d.%d", DirNum, Num, Num, iFile);
#endif

#ifdef HR1
  DirNum = hrdir(Num, iFile);
  sprintf(buf, "/data/milli%02d/hr1%02d/group_tab/group_tab_%03d.%d", DirNum, Num, Num, iFile);
#endif

  fd = my_fopen(buf, (char *) "r");

  my_fread(&Ngroups, sizeof(int), 1, fd);
  my_fread(&TotNgroups, sizeof(int), 1, fd);

#ifdef SWAPENDIAN
  byte_swap_long(1, &Ngroups);
  byte_swap_long(1, &TotNgroups);
#endif

  fclose(fd);
  return Ngroups;
}

long get_number_of_groups_total(char *filename, int Num)
{
  int Ngroups, TotNgroups, DirNum;
  char buf[1000];
  FILE *fd;

#ifdef MILLENNIUM
  DirNum = millenniumdir(Num, 0);
  sprintf(buf, "/data/milli%02d/d%02d/snapshot/snap_millennium_%03d.%d", DirNum, Num, Num, 0);
#endif

#ifdef HR1
  DirNum = hrdir(Num, 0);
  sprintf(buf, "/data/milli%02d/hr1%02d/group_tab/group_tab_%03d.%d", DirNum, Num, Num, 0);
#endif


  fd = my_fopen(filename, (char *) "r");

  my_fread(&Ngroups, sizeof(int), 1, fd);
  my_fread(&TotNgroups, sizeof(int), 1, fd);

#ifdef SWAPENDIAN
  byte_swap_long(1, &Ngroups);
  byte_swap_long(1, &TotNgroups);
#endif

  fclose(fd);
  return TotNgroups;
}

int get_total_number_of_groups(char *OutputDir, int Num)
{
  int DirNum, Ngroups, Nids, TotNgroups, NTask;
  char buf[1000];
  FILE *fd;

  sprintf(buf, "%s/group_tab_%03d.0", OutputDir, Num);

#ifdef MILLENNIUM
  DirNum = millenniumdir(Num, 0);
  sprintf(buf, "/data/milli%02d/d%02d/snapshot/snap_millennium_%03d.%d", DirNum, Num, Num, 0);
#endif

#ifdef HR1
  DirNum = hrdir(Num, 0);
  sprintf(buf, "/data/milli%02d/hr1%02d/group_tab/group_tab_%03d.%d", DirNum, Num, Num, 0);
#endif

  fd = my_fopen(buf, "r");

  my_fread(&Ngroups, sizeof(int), 1, fd);
  my_fread(&Nids, sizeof(int), 1, fd);
  my_fread(&TotNgroups, sizeof(int), 1, fd);
  my_fread(&NTask, sizeof(int), 1, fd);
  fclose(fd);

#ifdef SWAPENDIAN
  byte_swap_int(1, &TotNgroups);
#endif

  return TotNgroups;
}

int get_group_catalogue(char *OutputDir, int Num, int *GroupLen, int *GroupFileNr, int *GroupNr)
{
  int Ngroups, Nids, TotNgroups, NTask, count, i, filenr;
  cat_data *temp;
  char buf[1000];
  FILE *fd;

  sprintf(buf, "%s/group_tab_%03d.0", OutputDir, Num);
  if(!(fd = fopen(buf, "r")))
    {
      printf("can't open file `%s'\n", buf);
      return -1;
    }

  my_fread(&Ngroups, sizeof(int), 1, fd);
  my_fread(&Nids, sizeof(int), 1, fd);
  my_fread(&TotNgroups, sizeof(int), 1, fd);
  my_fread(&NTask, sizeof(int), 1, fd);
  fclose(fd);

  for(filenr = 0, count = 0; filenr < NTask; filenr++)
    {
      sprintf(buf, "%s/group_tab_%03d.%d", OutputDir, Num, filenr);
      if(!(fd = fopen(buf, "r")))
	{
	  printf("can't open file `%s'\n", buf);
	  return -1;
	}

      my_fread(&Ngroups, sizeof(int), 1, fd);
      my_fread(&Nids, sizeof(int), 1, fd);
      my_fread(&TotNgroups, sizeof(int), 1, fd);
      my_fread(&NTask, sizeof(int), 1, fd);

      my_fread(&GroupLen[count], sizeof(int), Ngroups, fd);

      for(i = 0; i < Ngroups; i++)
	{
	  GroupFileNr[i + count] = filenr;
	  GroupNr[i + count] = i;
	}

      count += Ngroups;

      fclose(fd);
    }
/*
  temp = my_malloc(sizeof(cat_data) * TotNgroups);

  for(i = 0; i < TotNgroups; i++)
    {
      temp[i].len = GroupLen[i];
      temp[i].file = GroupFileNr[i];
      temp[i].offset = GroupNr[i];
    }

  qsort(temp, TotNgroups, sizeof(cat_data), id_sort_groups);

  for(i = 0; i < TotNgroups; i++)
    {
      GroupLen[i] = temp[i].len;
      GroupFileNr[i] = temp[i].file;
      GroupNr[i] = temp[i].offset;
    }

  free(temp);
*/
  return TotNgroups;
}

int id_sort_groups(const void *a, const void *b)
{
  if(((cat_data *) a)->len > ((cat_data *) b)->len)
    return -1;

  if(((cat_data *) a)->len < ((cat_data *) b)->len)
    return +1;

  return 0;
}

int get_hash_table_size(char *OutputDir, int Num, char *SnapBase, int *NFiles)
{
  int dummy, DirNum;
  io_header header;
  char buf[1000];
  FILE *fd;

  sprintf(buf, "%s/%s_%03d.0", OutputDir, SnapBase, Num);

#ifdef MILLENNIUM
  DirNum = millenniumdir(Num, 0);
  sprintf(buf, "/data/milli%02d/d%02d/snapshot/snap_millennium_%03d.%d", DirNum, Num, Num, 0);
#endif

#ifdef HR1
  DirNum = hrdir(Num, 0);
  sprintf(buf, "/data/milli%02d/hr1%02d/snapshot/hr1_snap_%03d.%d", DirNum, Num, Num, 0);
#endif

  fd = my_fopen(buf, "r");

  my_fread(&dummy, sizeof(int), 1, fd);
  my_fread(&header, sizeof(io_header), 1, fd);
  my_fread(&dummy, sizeof(int), 1, fd);
  fclose(fd);

#ifdef SWAPENDIAN
  byte_swap_int(1, &(header.num_files));
  byte_swap_int(1, &(header.hashtabsize));
#endif

  *NFiles = header.num_files;

  return header.hashtabsize;
}

int get_hash_table(char *OutputDir, int Num, char *SnapBase, int *Hashtable, int *FileTable, int *LastHashCell, int *NInFiles)
{
  int DirNum, dummy, filenr, NTask;
  int i, firstcell, lastcell;
  io_header header;
  char buf[1000];
  FILE *fd;

  NTask = 1;			/* will be taken from header */

  for(filenr = 0; filenr < NTask; filenr++)
    {
      sprintf(buf, "%s/%s_%03d.%d", OutputDir, SnapBase, Num, filenr);

#ifdef MILLENNIUM
      DirNum = millenniumdir(Num, filenr);
      sprintf(buf, "/data/milli%02d/d%02d/snapshot/snap_millennium_%03d.%d", DirNum, Num, Num, filenr);
#endif

#ifdef HR1
      DirNum = hrdir(Num, filenr);
      sprintf(buf, "/data/milli%02d/hr1%02d/snapshot/hr1_snap_%03d.%d", DirNum, Num, Num, filenr);
#endif


      fd = my_fopen(buf, "r");

      fread(&dummy, sizeof(int), 1, fd);
      fread(&header, sizeof(io_header), 1, fd);
      fread(&dummy, sizeof(int), 1, fd);
#ifdef SWAPENDIAN
      byte_swap_int(1, &(header.num_files));
#endif
      NTask = header.num_files;

      /* skip positions */
      fread(&dummy, sizeof(int), 1, fd);
#ifdef SWAPENDIAN
      byte_swap_int(1, &(dummy));
#endif
      fseek(fd, dummy, SEEK_CUR);
      fread(&dummy, sizeof(int), 1, fd);

      /* skip velocities */
      fread(&dummy, sizeof(int), 1, fd);
#ifdef SWAPENDIAN
      byte_swap_int(1, &(dummy));
#endif
      fseek(fd, dummy, SEEK_CUR);
      fread(&dummy, sizeof(int), 1, fd);

      /* skip IDs */
      fread(&dummy, sizeof(int), 1, fd);
#ifdef SWAPENDIAN
      byte_swap_int(1, &(dummy));
#endif
      fseek(fd, dummy, SEEK_CUR);
      fread(&dummy, sizeof(int), 1, fd);

      fread(&dummy, sizeof(int), 1, fd);
      fread(&firstcell, sizeof(int), 1, fd);
      fread(&lastcell, sizeof(int), 1, fd);
      fread(&dummy, sizeof(int), 1, fd);

#ifdef SWAPENDIAN
      byte_swap_int(1, &(firstcell));
      byte_swap_int(1, &(lastcell));
      byte_swap_int(1, &(header.npart[1]));
      byte_swap_int(1, &(header.hashtabsize));
#endif

      NTask = header.num_files;
      LastHashCell[filenr] = lastcell;
      NInFiles[filenr] = header.npart[1];

      if(firstcell < 0 || firstcell >= header.hashtabsize || lastcell < 0 || lastcell >= header.hashtabsize)
	{
	  printf("bummer!\n");
	  return -1;
	}

      my_fread(&dummy, sizeof(int), 1, fd);
      my_fread(&Hashtable[firstcell], sizeof(int), lastcell - firstcell + 1, fd);
      my_fread(&dummy, sizeof(int), 1, fd);

#ifdef SWAPENDIAN
      byte_swap_int(lastcell - firstcell + 1, &Hashtable[firstcell]);
#endif


      for(i = 0; i < lastcell - firstcell + 1; i++)
	FileTable[firstcell + i] = filenr;

      fclose(fd);
    }

  return 0;
}
/*
int get_group_coordinates(char *OutputDir, int Num, char *Snapbase,
			  int *Hashtable, int *Filetable, int hashtabsize, int *LastHashCell,
			  int *NInFiles, int GroupNr, int FileNr, int GroupLen, float *x, float *y, float *z, float *vx, float *vy, float *vz)
{

  int dummy, NTask, ind;
  int i, j, n;
  io_header header;
  char buf[1000];
  FILE *fd;
  int Ngroups, TotNgroups, count, Nids;
  int len, offset;
  float *Pos, *Vel;
  double Sx, Sy, Sz;
  double sx, sy, sz;
  double svx, svy, svz;
  double refx, refy, refz;
  float *pos_read;
  float *vel_read;
  int hash_key;
  unsigned long long id;
  unsigned long long *ids, *id_read;


  sprintf(buf, "%s/group_tab_%03d.%d", OutputDir, Num, FileNr);
  if(!(fd = fopen(buf, "r")))
    {
      printf("can't open file `%s'\n", buf);
      return -1;
    }

  my_fread(&Ngroups, sizeof(int), 1, fd);
  my_fread(&Nids, sizeof(int), 1, fd);
  my_fread(&TotNgroups, sizeof(int), 1, fd);
  my_fread(&NTask, sizeof(int), 1, fd);
  // skip group-len table 
  fseek(fd, sizeof(int) * Ngroups, SEEK_CUR);
  fseek(fd, sizeof(int) * GroupNr, SEEK_CUR);
  my_fread(&offset, sizeof(int), 1, fd);
  fclose(fd);

  ids = my_malloc(GroupLen * sizeof(long long));

  sprintf(buf, "%s/group_ids_%03d.%d", OutputDir, Num, FileNr);
  if(!(fd = fopen(buf, "r")))
    {
      printf("can't open file `%s'\n", buf);
      return -1;
    }
  my_fread(&Ngroups, sizeof(int), 1, fd);
  my_fread(&Nids, sizeof(int), 1, fd);
  my_fread(&TotNgroups, sizeof(int), 1, fd);
  my_fread(&NTask, sizeof(int), 1, fd);

  fseek(fd, sizeof(long long) * offset, SEEK_CUR);
  my_fread(ids, sizeof(long long), GroupLen, fd);
  fclose(fd);

  //HashMap    = my_malloc( hashtabsize * sizeof(int) );
  // ids contains the id of the group's particles


  for(i = 0; i < hashtabsize; i++)
    Filetable[i] <<= 1;		// Filetable[i] = 2^1*Filetable[i];

  for(i = 0, count = 0; i < GroupLen; i++)
    {
      hash_key = (ids[i] >> 34);

      if(hash_key >= hashtabsize)
	{
	  printf("bummer\n");
	  return -1;
	}

      if((Filetable[hash_key] & 1) == 0)
	{
	  Filetable[hash_key] |= 1;

	  if(hash_key != LastHashCell[Filetable[hash_key] >> 1])
	    count += Hashtable[hash_key + 1] - Hashtable[hash_key];
	  else
	    count += NInFiles[Filetable[hash_key] >> 1] - Hashtable[hash_key];
	}
    }

  pos_read = my_malloc(3 * count * sizeof(float));
  vel_read = my_malloc(3 * count * sizeof(float));
  id_read = my_malloc(count * sizeof(long long));

  for(i = 0, count = 0; i < hashtabsize; i++)
    {
      if((Filetable[i] & 1))
	{
	  if(i != LastHashCell[Filetable[i] >> 1])
	    len = Hashtable[i + 1] - Hashtable[i];
	  else
	    len = NInFiles[Filetable[i] >> 1] - Hashtable[i];

	  sprintf(buf, "%s/%s_%03d.%d", OutputDir, Snapbase, Num, Filetable[i] >> 1);

	  fd = my_fopen(buf, "r");
	  my_fread(&dummy, sizeof(int), 1, fd);
	  my_fread(&header, sizeof(io_header), 1, fd);
	  my_fread(&dummy, sizeof(int), 1, fd);
	  my_fread(&dummy, sizeof(int), 1, fd);
	  fseek(fd, 3 * Hashtable[i] * sizeof(float), SEEK_CUR);
	  my_fread(&pos_read[3 * count], sizeof(float), 3 * len, fd);
	  fclose(fd);

	  fd = my_fopen(buf, "r");
	  my_fread(&dummy, sizeof(int), 1, fd);
	  my_fread(&header, sizeof(io_header), 1, fd);
	  my_fread(&dummy, sizeof(int), 1, fd);
	  // skip positions 
	  my_fread(&dummy, sizeof(int), 1, fd);
	  fseek(fd, dummy, SEEK_CUR);
	  my_fread(&dummy, sizeof(int), 1, fd);
	  my_fread(&dummy, sizeof(int), 1, fd);
	  fseek(fd, 3 * Hashtable[i] * sizeof(float), SEEK_CUR);
	  my_fread(&vel_read[3 * count], sizeof(float), 3 * len, fd);
	  fclose(fd);

	  fd = my_fopen(buf, "r");
	  my_fread(&dummy, sizeof(int), 1, fd);
	  my_fread(&header, sizeof(io_header), 1, fd);
	  my_fread(&dummy, sizeof(int), 1, fd);
	  // skip positions 
	  my_fread(&dummy, sizeof(int), 1, fd);
	  fseek(fd, dummy, SEEK_CUR);
	  my_fread(&dummy, sizeof(int), 1, fd);
	  // skip velocities 
	  my_fread(&dummy, sizeof(int), 1, fd);
	  fseek(fd, dummy, SEEK_CUR);
	  my_fread(&dummy, sizeof(int), 1, fd);
	  my_fread(&dummy, sizeof(int), 1, fd);
	  fseek(fd, Hashtable[i] * sizeof(long long), SEEK_CUR);
	  my_fread(&id_read[count], sizeof(long long), len, fd);
	  fclose(fd);

	  count += len;
	}
    }

  for(i = 0; i < GroupLen; i++)
    ids[i] &= ((((long long) 1) << 34) - 1);

  for(i = 0; i < count; i++)
    id_read[i] = ((id_read[i] & ((((long long) 1) << 34) - 1)) << 29) + i;

  qsort(ids, GroupLen, sizeof(long long), id_sort_compare_key);
  qsort(id_read, count, sizeof(long long), id_sort_compare_key);

  Pos = my_malloc(3 * GroupLen * sizeof(float));
  Vel = my_malloc(3 * GroupLen * sizeof(float));

  if(count < GroupLen)
    {
      printf("PROBLEM! count=%d GroupLen=%d\n", count, GroupLen);
      printf("Something seems to be wrong here=!\n");
    }

  for(i = 0, n = 0; i < GroupLen && n < count; i++)
    {
      id = (id_read[n] >> 29);

      while(n < count && id < ids[i])
	{
	  n++;
	  id = (id_read[n] >> 29);
	}

      if(id != ids[i])
	printf("We have a mismatch! (i=%d) something is wrong here.\n", i);

      ind = (id_read[n] & ((1 << 29) - 1));

      for(j = 0; j < 3; j++)
	{
	  Vel[i * 3 + j] = vel_read[ind * 3 + j];
	  Pos[i * 3 + j] = pos_read[ind * 3 + j];
	}
    }

  sx = sy = sz = 0;
  svx = svy = svz = 0;
  refx = Pos[0];
  refy = Pos[1];
  refz = Pos[2];

  for(i = 0; i < GroupLen; i++)
    {
      Pos[i * 3 + 0] = fof_periodic(Pos[i * 3 + 0] - refx, BoxSize);
      Pos[i * 3 + 1] = fof_periodic(Pos[i * 3 + 1] - refy, BoxSize);
      Pos[i * 3 + 2] = fof_periodic(Pos[i * 3 + 2] - refz, BoxSize);
    }

  for(i = 0; i < GroupLen; i++)
    {
      sx += Pos[i * 3 + 0];
      sy += Pos[i * 3 + 1];
      sz += Pos[i * 3 + 2];
      svx += Vel[i * 3 + 0];
      svy += Vel[i * 3 + 1];
      svz += Vel[i * 3 + 2];
    }

  sx /= GroupLen;
  sy /= GroupLen;
  sz /= GroupLen;
  svx /= GroupLen;
  svy /= GroupLen;
  svz /= GroupLen;

  for(i = 0; i < GroupLen; i++)
    {
      Pos[i * 3 + 0] -= sx;
      Pos[i * 3 + 1] -= sy;
      Pos[i * 3 + 2] -= sz;
    }

  for(i = 0; i < hashtabsize; i++)
    Filetable[i] >>= 1;		// Filetable[i] = Filetable[i] / 2^1 ;

  free(id_read);
  free(pos_read);
  free(vel_read);
  free(ids);
  free(Pos);
  free(Vel);

  *x = fof_periodic_wrap(sx + refx, BoxSize);
  *y = fof_periodic_wrap(sy + refy, BoxSize);
  *z = fof_periodic_wrap(sz + refz, BoxSize);

  *vx = svx;
  *vy = svy;
  *vz = svz;

  return 0;
}
*/

double fof_periodic_wrap(double x, double BoxSize)
{

  while(x >= BoxSize)
    x -= BoxSize;

  while(x < 0)
    x += BoxSize;

  if(x >= BoxSize)
    printf(" %f Somthing wrong here\n", x);
  if(x < 0)
    printf(" %f Somthing wrong here\n", x);

  return x;
}


double fof_periodic(double x, double BoxSize)
{
  if(x >= 0.5 * BoxSize)
    x -= BoxSize;

  if(x < -0.5 * BoxSize)
    x += BoxSize;

  return x;
}


int id_sort_compare_key(const void *a, const void *b)
{
  if(*((long long *) a) < *((long long *) b))
    return -1;

  if(*((long long *) a) > *((long long *) b))
    return +1;

  return 0;
}


// ----- peano-hilbert key related ----- //

static int quadrants[24][2][2][2] = {
  /* rotx=0, roty=0-3 */
  {{{0, 7}, {1, 6}}, {{3, 4}, {2, 5}}},
  {{{7, 4}, {6, 5}}, {{0, 3}, {1, 2}}},
  {{{4, 3}, {5, 2}}, {{7, 0}, {6, 1}}},
  {{{3, 0}, {2, 1}}, {{4, 7}, {5, 6}}},
  /* rotx=1, roty=0-3 */
  {{{1, 0}, {6, 7}}, {{2, 3}, {5, 4}}},
  {{{0, 3}, {7, 4}}, {{1, 2}, {6, 5}}},
  {{{3, 2}, {4, 5}}, {{0, 1}, {7, 6}}},
  {{{2, 1}, {5, 6}}, {{3, 0}, {4, 7}}},
  /* rotx=2, roty=0-3 */
  {{{6, 1}, {7, 0}}, {{5, 2}, {4, 3}}},
  {{{1, 2}, {0, 3}}, {{6, 5}, {7, 4}}},
  {{{2, 5}, {3, 4}}, {{1, 6}, {0, 7}}},
  {{{5, 6}, {4, 7}}, {{2, 1}, {3, 0}}},
  /* rotx=3, roty=0-3 */
  {{{7, 6}, {0, 1}}, {{4, 5}, {3, 2}}},
  {{{6, 5}, {1, 2}}, {{7, 4}, {0, 3}}},
  {{{5, 4}, {2, 3}}, {{6, 7}, {1, 0}}},
  {{{4, 7}, {3, 0}}, {{5, 6}, {2, 1}}},
  /* rotx=4, roty=0-3 */
  {{{6, 7}, {5, 4}}, {{1, 0}, {2, 3}}},
  {{{7, 0}, {4, 3}}, {{6, 1}, {5, 2}}},
  {{{0, 1}, {3, 2}}, {{7, 6}, {4, 5}}},
  {{{1, 6}, {2, 5}}, {{0, 7}, {3, 4}}},
  /* rotx=5, roty=0-3 */
  {{{2, 3}, {1, 0}}, {{5, 4}, {6, 7}}},
  {{{3, 4}, {0, 7}}, {{2, 5}, {1, 6}}},
  {{{4, 5}, {7, 6}}, {{3, 2}, {0, 1}}},
  {{{5, 2}, {6, 1}}, {{4, 3}, {7, 0}}}
};


static int rotxmap_table[24] = { 4, 5, 6, 7, 8, 9, 10, 11,
  12, 13, 14, 15, 0, 1, 2, 3, 17, 18, 19, 16, 23, 20, 21, 22
};

static int rotymap_table[24] = { 1, 2, 3, 0, 16, 17, 18, 19,
  11, 8, 9, 10, 22, 23, 20, 21, 14, 15, 12, 13, 4, 5, 6, 7
};

static int rotx_table[8] = { 3, 0, 0, 2, 2, 0, 0, 1 };
static int roty_table[8] = { 0, 1, 1, 2, 2, 3, 3, 0 };

static int sense_table[8] = { -1, -1, -1, +1, +1, -1, -1, -1 };

static int flag_quadrants_inverse = 1;
static char quadrants_inverse_x[24][8];
static char quadrants_inverse_y[24][8];
static char quadrants_inverse_z[24][8];




/*! This function computes a Peano-Hilbert key for an integer triplet (x,y,z),
 *  with x,y,z in the range between 0 and 2^bits-1.
 */
peanokey peano_hilbert_key(int x, int y, int z, int bits)
{
  int i, quad, bitx, bity, bitz;
  int mask, rotation, rotx, roty, sense;
  peanokey key;


  mask = 1 << (bits - 1);
  key = 0;
  rotation = 0;
  sense = 1;


  for(i = 0; i < bits; i++, mask >>= 1)
    {
      bitx = (x & mask) ? 1 : 0;
      bity = (y & mask) ? 1 : 0;
      bitz = (z & mask) ? 1 : 0;

      quad = quadrants[rotation][bitx][bity][bitz];

      key <<= 3;
      key += (sense == 1) ? (quad) : (7 - quad);

      rotx = rotx_table[quad];
      roty = roty_table[quad];
      sense *= sense_table[quad];
      while(rotx > 0)
	{
	  rotation = rotxmap_table[rotation];
	  rotx--;
	}

      while(roty > 0)
	{
	  rotation = rotymap_table[rotation];
	  roty--;
	}
    }

  return key;
}


void peano_hilbert_key_inverse(peanokey key, int bits, int *x, int *y, int *z)
{
  int i, keypart, bitx, bity, bitz, mask, quad, rotation, shift;
  char sense, rotx, roty;

  if(flag_quadrants_inverse)
    {
      flag_quadrants_inverse = 0;
      for(rotation = 0; rotation < 24; rotation++)
	for(bitx = 0; bitx < 2; bitx++)
	  for(bity = 0; bity < 2; bity++)
	    for(bitz = 0; bitz < 2; bitz++)
	      {
		quad = quadrants[rotation][bitx][bity][bitz];
		quadrants_inverse_x[rotation][quad] = bitx;
		quadrants_inverse_y[rotation][quad] = bity;
		quadrants_inverse_z[rotation][quad] = bitz;
	      }
    }

  shift = 3 * (bits - 1);
  mask = 7 << shift;

  rotation = 0;
  sense = 1;

  *x = *y = *z = 0;

  for(i = 0; i < bits; i++, mask >>= 3, shift -= 3)
    {
      keypart = (key & mask) >> shift;

      quad = (sense == 1) ? (keypart) : (7 - keypart);

      *x = (*x << 1) + quadrants_inverse_x[rotation][quad];
      *y = (*y << 1) + quadrants_inverse_y[rotation][quad];
      *z = (*z << 1) + quadrants_inverse_z[rotation][quad];

      rotx = rotx_table[quad];
      roty = roty_table[quad];
      sense *= sense_table[quad];

      while(rotx > 0)
	{
	  rotation = rotxmap_table[rotation];
	  rotx--;
	}

      while(roty > 0)
	{
	  rotation = rotymap_table[rotation];
	  roty--;
	}
    }
}


