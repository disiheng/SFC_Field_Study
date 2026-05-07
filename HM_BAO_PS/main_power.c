/*
This code computes the counts-in-cells moments for dark matter and haloes.
It also computes the cross-correlations of the fields.
*/

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


#include <gsl/gsl_rng.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_histogram2d.h>
#include <srfftw_mpi.h>

#include "utils.h"
#include "mymalloc.h"

//#include "gadget.h"
//#include "read_data.h"
//#include "L-Galaxies_bruno3.h"

#include "density.h"
#include "allvars.h"
#include "power.h"
#include "correl.h"
#include "io.h"

#include <hdf5.h>

void read_parameter_file(char *fname);
int mstar_compare(const void *a, const void *b);

int main(int argc, char *argv[])
{
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &NTasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);

  //FILE *fd, *fd_star;
  //int H_NPart, DM_NPart;
  //char DM_File[255], H_File[255];
 
  int  i, j, k, ii, nbins;
  char sel_type[100], sel_space[100];
  char buf[512], buf2[512], fname[256];
  int  nr;
  long long thisNP;
  long long TotNumPart;
  int foldnum;

  //float ndens[] = {1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e2};
  //float ndens[] = {1.16e-3, 6.5e-3};
 
  float ndens[] = {100.0};
 
  long long ntot, cum;
  //struct stat buffd;

  Grid DM_grid, Vx_grid, Vy_grid, Vz_grid, Theta_grid, grid_x;

  t0 = time(NULL);

  if (ThisTask == 0)
    printf("\nCode was compiled with settings:\n %s\n\n", COMPILETIMESETTINGS);

  read_parameter_file(argv[1]);

  HighMark = 0;
  mymalloc_init();

  FileFormat = atoi(argv[2]);
  nbins = 1;

  if(FileFormat == 0)
    {
      sprintf(sel_type,"DM");
      nbins = 1;

#ifdef SUBSAMPLE
      nbins = 1;
#endif
    }

  if (FileFormat == 1) sprintf(sel_type,"Gals_StellarMass");         // stellar mass
  if (FileFormat == 2) sprintf(sel_type,"Gals_StarFormationRate");   // star formation rate

  if (FileFormat == 3) sprintf(sel_type,"Halo_Mvir");                // host halo virial mass

  //if (FileFormat == 7) sprintf(sel_type,"tree_halo");
  //if (FileFormat == 7) sprintf(sel_type,"tree_halo_sub");
  if (FileFormat == 7) sprintf(sel_type,"tree_halo_sub2");
  //if (FileFormat == 7) sprintf(sel_type,"halo_vmax");

#ifdef ITERVOL
  SubVolDim = (((peanokey) 1) << (ITERVOL_BITS_PER_DIMENSION));

  if(argc == 4)
    iDiv = atoi(argv[3]);
  else
    iDiv = 0;

  sprintf(sel_type,"%s_SubVol%03d",sel_type, iDiv);
  mark_domain_center(iDiv);
#else
  SubVolDim = 1;
#endif


#ifdef GAL_SHUFFLE
  sprintf(sel_type,"%s_shuffled",sel_type);  
  if(argc == 4)
    sprintf(sel_type,"%s%d",sel_type, atoi(argv[3]));  
#endif


#ifdef Z_SPACE
  sprintf(sel_space,"zspace");
#ifdef PHOTO_SPACE
  sprintf(sel_space,"pzspace");
#endif
#else
  sprintf(sel_space,"rspace");
#endif


  MyBoxSize = (double *) mymalloc("MyBoxSize", (FOLDS + 1) * sizeof(double));
  for(i = 0; i <= FOLDS; ++i)
    MyBoxSize[i] = BoxSize / SubVolDim / pow(FOLD_FACTOR, i);

#ifdef PRELOAD
  load_catalogue(); 

  //-----------------------------------------------
  float minsep[3], maxsep[3], glob_minsep[3], glob_maxsep[3];

  for(j=0;j<3;j++)
  {
    minsep[j] = 1e20;
    maxsep[j] = -1e20;
  }

  for(i = 0; i < LocNgroups; i++)
  {
    for(k = 0; k < 3; k++)
    {
      if(PP[i].pos[k] > maxsep[k])
        maxsep[k] = PP[i].pos[k];
      if(PP[i].pos[k] < minsep[k])
        minsep[k] = PP[i].pos[k];
    }
  }

  MPI_Reduce(minsep, glob_minsep, 3, MPI_FLOAT, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce(maxsep, glob_maxsep, 3, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
  {
    printf("\n");
    for(j = 0; j < 3; j++)
      printf("  minsep[%d]=%f maxsep[%d]=%f \n", j, glob_minsep[j], j, glob_maxsep[j]);
    printf("\n");
  }

  //if (ThisTask == 0) printf("Parallel sorting\n"); 
  //parallel_sort(PP, LocNgroups, sizeof(Particle), mstar_compare);

  int *gals_per_task = mymalloc("gals_per_task",sizeof(int)*NTasks);
  MPI_Allgather(&LocNgroups, 1, MPI_INT, gals_per_task, 1, MPI_INT, MPI_COMM_WORLD);
#endif

#ifdef CROSS
  grid_x = init_density(NCell,0);

  sprintf(buf,"%s/DM_field_%s_%d/DM_field.%d", CrossOutputDir, sel_space, NCell, ThisTask);
  density_load(buf, &grid_x);
#endif

  for (i = 0; i < nbins; i++)
  {
    for (foldnum=0; foldnum <= FOLDS; foldnum++)
    {
      DM_grid = init_density(NCell,foldnum);
      ntot  = ndens[i]*CUB(BoxSize);

#ifdef DIVERGENCE
      Theta_grid = init_density(NCell,foldnum);
      Vx_grid = init_density(NCell,foldnum);
      Vy_grid = init_density(NCell,foldnum);
      Vz_grid = init_density(NCell,foldnum);
#endif

#ifdef PRELOAD
      for(cum = 0, ii = 0; ii < ThisTask; ii++)
        cum += gals_per_task[ii];

      thisNP = ntot - cum;
      if (thisNP > LocNgroups) 
         thisNP = LocNgroups;
      if (thisNP < 0) 
         thisNP = 0;

      if (thisNP > 0 && thisNP < LocNgroups)
        printf("Minimum value in task %d, equal to %f\n", ThisTask, PP[thisNP-1].mass);

      //thisNP =  LocNgroups;
      //printf("Selecting %d of %lld gals \n", thisNP, ntot);   
      //dump_particles(thisNP);

      { 
        density (PP, thisNP, &DM_grid, 1, 0);
#ifdef DIVERGENCE
        density (PP, LocNgroups, &Vx_grid, 1, 1);
        density (PP, LocNgroups, &Vy_grid, 1, 2);
        density (PP, LocNgroups, &Vz_grid, 1, 3);
#endif
        density_statistics(&DM_grid,  "Density grid",0);
        overdensity(&DM_grid);
      }
#endif

#ifdef DM
#ifdef DENSITY_LOAD
      // loading the presave density field outputed by LG3 
      //load_field(&DM_grid);
      //density_statistics(&DM_grid,  "Density grid",0);
      //overdensity(&DM_grid);

      // loading the density field outputed with OPT-DENSITY_DUMP of this code 
      sprintf(buf,"%s/DM_field_%03d_%s_%d/DM_field.%d", CrossOutputDir, SnapNum, sel_space, NCell, ThisTask);
      density_load(buf, &DM_grid);
#else //DENSITY_LOAD

      int ngroups = ceil(NFiles*1./NTasks);
      int ig, readgroups = floor(NTasks/NumFilesWrittenInParallel);
      float nsel  = ndens[i] * CUB(BoxSize);

      for(ii = 0; ii < ngroups; ii ++)
      {
        nr = ii*NTasks + ThisTask;
        for(ig = 0; ig < readgroups; ig++)
        {
          if(ThisTask == 0) printf("reading group = %d of %d in %d cycles\n",ig, readgroups, ngroups);
          if(ThisTask % readgroups == ig)
          { 
            LocNgroups = 0;
            if(nr < NFiles )
            {
#ifdef HDF5_G4
              load_dm_single_hdf5(nr, nsel);
#else
              load_dm_single(nr, nsel);
#endif //HDF5_G4
              //printf("."); fflush(stdout);
            } 
          }
          MPI_Barrier(MPI_COMM_WORLD);
        }
        //if(ThisTask == 0)
        //printf("Task:%d Have Reading = %d particles in thisTask\n", ThisTask, LocNgroups);

        density (PP, LocNgroups, &DM_grid, 1, 0);
#ifdef DIVERGENCE
        density (PP, LocNgroups, &Vx_grid, 1, 1);
        density (PP, LocNgroups, &Vy_grid, 1, 2);
        density (PP, LocNgroups, &Vz_grid, 1, 3);
#endif //DIVERGENCE
        if(LocNgroups > 0)
           myfree(PP);
        PP = NULL;
      }
      density_statistics(&DM_grid,  "Density grid",0);
      overdensity(&DM_grid);
#endif //DENSITY_LOAD
#endif //DM

#ifdef DIVERGENCE
      density_statistics(&Vx_grid,  "Vx grid",0);
      density_statistics(&Vy_grid,  "Vy grid",0);
      density_statistics(&Vz_grid,  "Vz grid",0);

      divergence(DM_grid, Vx_grid, Vy_grid, Vz_grid, &Theta_grid,  MyBoxSize[foldnum]);
      grid_clean(Vz_grid);
      grid_clean(Vy_grid);
      grid_clean(Vx_grid);

      sprintf(buf2,"%s/gals_%s_divergence_bin%d",OutputDir, SnapshotFileBase, i);
      power(buf2, Theta_grid, foldnum, MyBoxSize[foldnum]);
      grid_clean(Theta_grid);
#endif //DIVERGENCE

#if defined(DENSITY_DUMP) && defined(DM)
      sprintf(buf2,"%s/DM_field_%03d_%s_%d/", CrossOutputDir, SnapNum, sel_space, NCell);
      mkdir(buf2,02755);
      sprintf(buf2,"%s/DM_field_%03d_%s_%d/DM_field.%d", CrossOutputDir, SnapNum, sel_space, NCell, ThisTask);
      density_save(buf2, DM_grid);
#endif //DENSITY_DUMP

      sprintf(buf2,"%s/%s_%s_%s_bin%d",OutputDir, SimuTag, sel_space, sel_type, i);
      if (SnapNum != -1)
        sprintf(buf2,"%s/%s_%03d_%s_%s_bin%d",OutputDir, SimuTag, SnapNum, sel_space, sel_type, i);
      power(buf2, DM_grid, grid_x, foldnum, MyBoxSize[foldnum]);

#ifdef CORREL
      correl(buf2, DM_grid, foldnum, MyBoxSize[foldnum]);
#endif

      grid_clean(DM_grid);
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);

#ifdef CROSS
  grid_clean(grid_x);
#endif

#ifdef PRELOAD
  myfree(gals_per_task);
  //if(LocNgroups > 0)
    myfree(PP);

#if defined(ITERVOL) && defined(HDF5_Phs)
  myfree(IterstoLoad);
  myfree(CataFileMap);
#endif
#endif

  myfree(MyBoxSize);
  MPI_Finalize();

  return 0;
}


int mstar_compare(const void *a, const void *b)
{

  if(((Particle *) a)->mass < ((Particle *) b)->mass)
    return +1;

  if(((Particle *) a)->mass > ((Particle *) b)->mass)
    return -1;

 return 0;
}


void read_parameter_file(char *fname)
{
#define DOUBLE 1
#define STRING 2
#define INT 3
#define MAXTAGS 300

  FILE *fd;
  char dumm[NMAXCHAR];
  char tag[MAXTAGS][50];
  char buf[200], buf1[200], buf2[200], buf3[400];
  int id[MAXTAGS];
  int nt, i, j;
  void *addr[MAXTAGS];

  nt = 0;

  strcpy(tag[nt], "BoxSize");
  addr[nt] = &BoxSize;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "SnapNum");
  addr[nt] = &SnapNum;
  id[nt++] = INT;

  strcpy(tag[nt], "NumFilesWrittenInParallel");
  addr[nt] = &NumFilesWrittenInParallel;
  id[nt++] = INT;

  strcpy(tag[nt], "SimuTag");
  addr[nt] = SimuTag;
  id[nt++] = STRING;

  strcpy(tag[nt], "SimulationDir");
  addr[nt] = SimulationDir;
  id[nt++] = STRING;

  strcpy(tag[nt], "CrossOutputDir");
  addr[nt] = CrossOutputDir;
  id[nt++] = STRING;

  strcpy(tag[nt], "OutputDir");
  addr[nt] = OutputDir;
  id[nt++] = STRING;

  strcpy(tag[nt], "SnapshotFileBase");
  addr[nt] = SnapshotFileBase;
  id[nt++] = STRING;

  strcpy(tag[nt], "FileFormat");
  addr[nt] = &FileFormat;
  id[nt++] = INT;

  strcpy(tag[nt], "flag_redshift");
  addr[nt] = &flag_redshift;
  id[nt++] = INT;

  strcpy(tag[nt], "NFiles");
  addr[nt] = &NFiles;
  id[nt++] = INT;

  strcpy(tag[nt], "NCell");
  addr[nt] = &NCell;
  id[nt++] = INT;

  strcpy(tag[nt], "MaxMemSize");
  addr[nt] = &MaxMemSize;
  id[nt++] = INT;

  strcpy(tag[nt], "Omega_m");
  addr[nt] = &Omega_m;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "Redshift");
  addr[nt] = &Redshift;
  id[nt++] = DOUBLE;

  if((fd = fopen(fname, "r")))
    {
      if(ThisTask == 0)
        printf("\n\n CURRENT CONFIGURATION in %s: \n", fname);
      while(!feof(fd))
        {
          *buf = 0;
          fgets(buf, 200, fd);
          if(sscanf(buf, "%s%s%s", buf1, buf2, buf3) < 2)
            continue;

          if(buf1[0] == '%')	//comments
            continue;

          for(i = 0, j = -1; i < nt; i++)
            if(strcmp(buf1, tag[i]) == 0)
              {
                j = i;
                tag[i][0] = 0;
                break;
              }

          if(j >= 0)
            {
              switch (id[j])
                {
                  case DOUBLE:
                    *((double *) addr[j]) = atof(buf2);
                    if(ThisTask == 0)
                      printf("  - %-35s%g\n", buf1, *((double *) addr[j]));
                    break;
                  case STRING:
                    strcpy(addr[j], buf2);
                    if(ThisTask == 0)
                      printf("  - %-35s%s\n", buf1, buf2);
                    break;
                  case INT:
                    *((int *) addr[j]) = atoi(buf2);
                    if(ThisTask == 0)
                      printf("  - %-35s%d\n", buf1, *((int *) addr[j]));
                    break;
                }
            }
          else
            {
              fprintf(stdout, "Error in file %s:   Tag '%s' not allowed or multiple defined.\n", fname, buf1);
              exit(0);
            }
        }
      fclose(fd);
    }
  else
    {
      printf("Parameter file %s not found.\n", fname);
      exit(0);
    }


  for(i = 0; i < nt; i++)
    {
      if(*tag[i])
        {
          printf("Error. I miss a value for tag '%s' in parameter file '%s'.\n", tag[i], fname);
          exit(0);
        }
    }

  if(FileFormat == 0)		//Means gadget
#if defined (HDF5_G4)
    {
      hid_t    hdf5_file, attr;
      herr_t   ret;

      sprintf(buf, "%s/%s_%03d.%d.hdf5", SimulationDir, SnapshotFileBase, SnapNum, 0);
      if (!(fd = fopen(buf, "r")))
        {
          sprintf(buf, "%s/snapdir_%03d/%s_%03d.%d.hdf5", SimulationDir, SnapNum, SnapshotFileBase, SnapNum, 0);
          fd = my_fopen(buf,"r");
        }

	  //printf("I'm in parameter file '%s' \n", buf);
      double BoxSizeR, RedshiftR, Omega0R;
      int num_files;

      hdf5_file = H5Fopen(buf, H5F_ACC_RDONLY, H5P_DEFAULT);
      hid_t ghead = H5Gopen(hdf5_file, "/Header", H5P_DEFAULT);

      attr = H5Aopen(ghead, "BoxSize", H5P_DEFAULT);
      ret  = H5Aread(attr, H5T_NATIVE_DOUBLE, &BoxSizeR);
      ret =  H5Aclose(attr);

      attr = H5Aopen(ghead, "Redshift", H5P_DEFAULT);
      ret  = H5Aread(attr, H5T_NATIVE_DOUBLE, &RedshiftR);
      ret =  H5Aclose(attr);

      attr = H5Aopen(ghead, "NumFilesPerSnapshot", H5P_DEFAULT);
      ret  = H5Aread(attr, H5T_NATIVE_INT, &num_files);
      ret =  H5Aclose(attr);
      H5Gclose(ghead);

      hid_t gparam = H5Gopen(hdf5_file, "/Parameters", H5P_DEFAULT);

      attr = H5Aopen(gparam, "Omega0", H5P_DEFAULT);
      ret  = H5Aread(attr, H5T_NATIVE_DOUBLE, &Omega0R);
      ret =  H5Aclose(attr);
      H5Gclose(gparam);

      //if(ThisTask == 0)
      //  {
      //     printf(" Header infos %g %g %d \n", BoxSizeR, Omega0R, num_files);
      //  }

      if(Redshift != RedshiftR)
        {
          if(ThisTask == 0)
            {
               printf("\n WARNING: Redshift %f %f(header) \n", Redshift, RedshiftR);
               printf(" Using that from the header (%s)\n", buf);
            }
          Redshift = RedshiftR;
        }

      if(BoxSize != BoxSizeR && NGens == 0)
        printf("\n WARNING: BoxSize %.8f %.8f \n", BoxSize, BoxSizeR);

      if(Omega_m != Omega0R)
        printf("\n WARNING: Omega_m %f %f \n", Omega_m, Omega0R);

      if(NFiles != num_files)
        printf("\n WARNING: NFiles %d %d \n", NFiles, num_files);
    }
#elif defined (HDF5_Phs)
    {
      double BoxSizeR, RedshiftR, Omega0R;
      int num_files, dummy;
      double dummy_dble;

      load_file_tab();

      sprintf(buf, "%s/data/snapshot_%d/cfg_%d.%d", fname, FileNumTab[SnapNum], FileNumTab[SnapNum], 0);
      if((fd = fopen(buf, "r")))
        {
          fscanf(fd, "%d"   , &dummy);
          fscanf(fd, "%d"   , &dummy);

          fscanf(fd, "%lld" , &dummy_dble);
          fscanf(fd, "%lg"  , &dummy_dble);
          fscanf(fd, "%lg"  , &BoxSizeR);
          fscanf(fd, "%lg"  , &dummy_dble);
          fscanf(fd, "%lg"  , &Omega0R);
          fscanf(fd, "%lg"  , &dummy_dble);

          fscanf(fd, "%lg"  , &RedshiftR);
          fclose(fd);
        }

      if(ThisTask == 0)
        printf(" Header infos %g %g %d \n", BoxSizeR, Omega0R, num_files);

      if(Redshift != RedshiftR)
        {
          if(ThisTask == 0)
            {
               printf("\n WARNING: Redshift %f %f(header) \n", Redshift, RedshiftR);
               printf(" Using that from the header (%s)\n", buf);
            }
          Redshift = RedshiftR;
        }

      if(BoxSize != BoxSizeR && NGens == 0)
        printf("\n WARNING: BoxSize %.8f %.8f \n", BoxSize, BoxSizeR);

      if(Omega_m != Omega0R)
        printf("\n WARNING: Omega_m %f %f \n", Omega_m, Omega0R);

      if(NFiles != num_files)
        printf("\n WARNING: NFiles %d %d \n", NFiles, num_files);
    }
#else
    {
      // We shall do some checks
      int dummy;
      io_header header;

      sprintf(buf, "%s/%s_%03d.%d", SimulationDir, SnapshotFileBase, SnapNum, 0);
      if(!(fd = fopen(buf, "r")))
        {
          sprintf(buf, "%s/snapdir_%03d/%s_%03d.%d", SimulationDir, SnapNum, SnapshotFileBase, SnapNum, 0);
          fd = my_fopen(buf, "r");
        }
      my_fread(&dummy, sizeof(int), 1, fd);
      my_fread(&header, sizeof(io_header), 1, fd);
      fclose(fd);

#ifdef SWAPENDIAN
      byte_swap_double(1, &(header.redshift));
      byte_swap_int(1, &(header.num_files));
      byte_swap_double(1, &(header.BoxSize));
      byte_swap_double(1, &(header.Omega0));
#endif

      header.redshift = 1. / header.time - 1.0;
      if(Redshift != header.redshift)
        {
          if (ThisTask == 0)
          {
            printf("\n WARNING: Redshift %f %f(header) \n", Redshift, header.redshift);
            printf(" Using that from the header (%s)\n", buf);
          }
          Redshift = header.redshift;
        }

      if(BoxSize != header.BoxSize && NGens == 0)
        printf("\n WARNING: BoxSize %.8f %.8f \n", BoxSize, header.BoxSize);

      if(Omega_m != header.Omega0)
        printf("\n WARNING: Omega_m %f %f \n", Omega_m, header.Omega0);

      if(NFiles != header.num_files)
        printf("\n WARNING: NFiles %d %d \n", NFiles, header.num_files);
    }
#endif

  ExpFactor = 1. / (1 + Redshift);
  Hubble = 100. * sqrt(Omega_m * CUB(1 + Redshift) + (1 - Omega_m));

#undef DOUBLE
#undef STRING
#undef INT
#undef MAXTAGS

}


