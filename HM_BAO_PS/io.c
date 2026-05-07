#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <mpi.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_histogram2d.h>

#ifdef DOUBLEPRECISION_FFTW
#include     <drfftw_mpi.h>
#else
#include     <srfftw_mpi.h>
#endif


#include "allvars.h"
#include "read_data.h"
#include "mymalloc.h"
#include "utils.h"
#include "density.h"
#include "io.h"

#include <hdf5.h>

void load_catalogue()       /*ok */
{
  int nr, num_files, i, ig, readgroups;
  long long totalngals=0;

#ifdef ITERVOL

#ifdef HDF5_Phs
  mark_Iter_nums();
  num_files = NumIterstoLoad;
#endif

#else
  num_files = NFiles;
#endif

  readgroups = floor(NTasks/NumFilesWrittenInParallel);

  //if(NTasks > num_files)
  //  terminate("one may at most use as many MPI processes as files in the group catalogue");

  if (ThisTask == 0)
    printf("\n loading catalogue. %d files in %d groups\n",num_files, readgroups);

  LocNgroups = 0;

  for(ig = 0; ig < readgroups; ig++)
  {
    if(ThisTask == 0)
      printf("\n reading group = %d \n",ig);

    for(nr = 0; nr < num_files; nr++)
      if((nr % NTasks) == ThisTask && (ThisTask % readgroups == ig))
      {
        if(FileFormat == 3)                     // halos
          {
#if !defined(HDF5_G4) && !defined(HDF5_Phs)
            load_halo_catalogue_single(nr);
#else

#if defined(HDF5_G4)
            load_halo_catalogue_single_hdf5(nr);
#endif

#ifdef HDF5_Phs
#ifdef ITERVOL
            load_halo_catalogue_single_photons(IterstoLoad[nr]);
#else
            load_halo_catalogue_single_photons(nr);
#endif
#endif

#endif
          } 
        else if(FileFormat == 7) 				// simple catalogue
          {
            //load_simple_catalogue_single(nr);
            //load_simple_catalogue_single2(nr + 832);
            load_simple_catalogue_single2(nr + 584);
          }
        else 							        // sam galaxy catalogue  with 1 2 6
          {
            load_catalogue_single(nr);
            //printf("."); fflush(stdout);
          }
      }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  if(PP != NULL)
    PP = (Particle *) myrealloc_movable(PP, LocNgroups * sizeof(Particle));

  MPI_Allreduce(&LocNgroups, &totalngals, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
  //sumup_large_ints(1, &LocNgroups, &totalngals);

  if(ThisTask == 0)
    printf("\n Number of objects: Local %d, Total %lld \n",  LocNgroups, totalngals);

#ifdef GAL_SHUFFLE
  float central_pos[3], *new_central_pos;
  float central_mass, min_mass = 0.1, max_mass=50000;
  int i_central, new_i_central, ibin, j, k, temp, nbins=20;
  int count_mbin[20], counter[20];

  gsl_rng *random_generator;

  int **indices;

  gsl_rng_env_setup();
  random_generator = gsl_rng_alloc(gsl_rng_rand48);
  long seed = time (NULL) * getpid();
  gsl_rng_set (random_generator, seed);

  for (i = 0; i < nbins; i++)
    count_mbin[i] = counter[i] = 0;

  for (i=0; i < LocNgroups; i++)
  {
    if (PP[i].Type == 0)
    {
     ibin = floor((log(PP[i].HaloMass)-log(min_mass))/(log(max_mass)-log(min_mass))*(nbins-1));   
     if (ibin >= nbins) ibin = nbins-1;
     count_mbin[ibin]++;
    }
  }

  indices = (int **) mymalloc("indices", nbins * sizeof(*indices));
  for (i = 0; i < nbins; i++)
    indices[i] = (int *)  mymalloc("indices", count_mbin[i] * sizeof(int));

  for (i=0; i < LocNgroups; i++)
  {
    if (PP[i].Type == 0)
    {
      ibin = floor((log(PP[i].HaloMass)-log(min_mass))/(log(max_mass)-log(min_mass))*(nbins-1));   
      if (ibin >= nbins) ibin = nbins-1;
      indices[ibin][counter[ibin]++] = i;
    }
  }

  for(ibin = 0; ibin < nbins; ibin++)
    assert(counter[ibin] ==  count_mbin[ibin]);

  new_central_pos = mymalloc("new_central_pos", 3 * LocNgroups * sizeof(float));

  //shuffling
  for (ibin = 0; ibin < nbins; ibin++)
  {
    for (i = counter[ibin]-1; i >=1; i--)
    {
      j = gsl_rng_uniform_int(random_generator,i);
      temp = indices[ibin][i];
      indices[ibin][i] = indices[ibin][j];  
      indices[ibin][j] = temp;  
    }
  }

  for (ibin = 0; ibin < nbins; ibin++)
    counter[ibin] = 0;

  for (i = 0; i < LocNgroups; i++)
  {
    if (PP[i].Type == 0)
    {
      ibin = floor((log(PP[i].HaloMass)-log(min_mass))/(log(max_mass)-log(min_mass))*(nbins-1));   
      if (ibin >= nbins) ibin = nbins-1;

//    new_i_central = gsl_rng_uniform(random_generator)*(counter[ibin]-1); 
      new_i_central = counter[ibin];
      counter[ibin]++;

      for(k = 0; k < 3; k++)
        new_central_pos[3 * i + k] =  PP[indices[ibin][new_i_central]].pos[k];
    }
  }

  for (i = 0; i < LocNgroups; i++)
  {
    if (PP[i].Type == 0)
    {
      i_central = i;
      central_mass = PP[i].HaloMass;
      for(k = 0; k < 3; k++)
        central_pos[k] =  PP[i].pos[k];
    }
 
/* 
    for(k = 0; k < 3; k++)
    {
      if(fabs(shuffle_periodic(PP[i].pos[k] - central_pos[k]) > 50))  
      {    
        printf("i=%d %d %f (%f|%f|%f)\n",i, PP[i].Type, PP[i].HaloMass, PP[i].pos[0],PP[i].pos[1],PP[i].pos[2]);
        printf("        %f (%f|%f|%f)\n",central_mass,central_pos[0],central_pos[1],central_pos[2]);
      }
    }
*/
    if (fabs(central_mass-PP[i].HaloMass) < 0.01)
    {
      for(k = 0; k < 3; k++)
      {
        PP[i].pos[k] = shuffle_periodic(PP[i].pos[k] - central_pos[k]) + new_central_pos[3*i_central + k];
        if (PP[i].pos[k] >= BoxSize) PP[i].pos[k] -= BoxSize;
        if (PP[i].pos[k] < 0)        PP[i].pos[k] += BoxSize;       
      }
    }
  }
  
  myfree(new_central_pos);

  for(i = nbins - 1; i >= 0; i--)
     myfree(indices[i]);

  myfree(indices);

  gsl_rng_free(random_generator);
#endif
}


double shuffle_periodic(double x)
{
  if(x >= 0.5 * BoxSize)
    x -= BoxSize;
  if(x < -0.5 * BoxSize)
    x += BoxSize;
  return x;
}



void load_bruno_single( int nr)
{
  int j=0, i, k, rep, dummy, thisN, num_files, ngroup, totngroup, ntrees, *galsPerTree;;
  char buf[512];
  FILE *fd;
  struct data_struct *halo_read;
  float *xx, *yy, *zz;
  gsl_rng *random_generator;

  gsl_rng_env_setup();
  random_generator = gsl_rng_alloc(gsl_rng_default);

  sprintf(buf,"%s/%s",SimulationDir,SnapshotFileBase);
  //printf("%s\n",buf);
  fd = my_fopen(buf,"r");
  my_fread(&dummy, sizeof(int),1,fd);
  my_fread(&ngroup, sizeof(int),1,fd);
  my_fread(&dummy, sizeof(int),1,fd);
 printf("%d\n",ngroup);
  xx = mymalloc_movable(&xx,"xx",sizeof(float)*ngroup);
  yy = mymalloc_movable(&yy,"yy",sizeof(float)*ngroup);
  zz = mymalloc_movable(&zz,"zz",sizeof(float)*ngroup);

  my_fread(&dummy, sizeof(int),1,fd);
  my_fread(xx, ngroup*sizeof(float),1,fd);
  my_fread(&dummy, sizeof(int),1,fd);

  my_fread(&dummy, sizeof(int),1,fd);
  my_fread(yy, ngroup*sizeof(float),1,fd);
  my_fread(&dummy, sizeof(int),1,fd);

  my_fread(&dummy, sizeof(int),1,fd);
  my_fread(zz, ngroup*sizeof(float),1,fd);
  my_fread(&dummy, sizeof(int),1,fd);
  fclose(fd);
  printf("done\n");

  if(PP == NULL)
    PP = mymalloc_movable(&PP, "PP", (LocNgroups + ngroup) * sizeof(Particle));
  else
    PP = myrealloc_movable(PP, (LocNgroups + ngroup) * sizeof(Particle));

   for(i = 0; i < ngroup; i++)
   {
       PP[LocNgroups].pos[0] = xx[i];
       PP[LocNgroups].pos[1] = yy[i];
       PP[LocNgroups].pos[2] = zz[i];
//     printf("(%f|%f|%f) Box=%f\n",PP[LocNgroups].pos[0],PP[LocNgroups].pos[1],PP[LocNgroups].pos[2], BoxSize);

       for (k = 0; k < 3; k++) 
       {
         if (PP[LocNgroups].pos[k] >= BoxSize) PP[LocNgroups].pos[k] -= BoxSize;
         if (PP[LocNgroups].pos[k] < 0)  PP[LocNgroups].pos[k] += BoxSize;       

         if (PP[LocNgroups].pos[k] < 0 || PP[LocNgroups].pos[k] >= BoxSize)
           printf("(%f|%f|%f) Box=%f\n",PP[LocNgroups].pos[0],PP[LocNgroups].pos[1],PP[LocNgroups].pos[2], BoxSize);
       }
       //PP[LocNgroups].mass = gsl_rng_uniform(random_generator); 

       //if (gsl_rng_uniform(random_generator) < 0.5) 
        LocNgroups++;
     }
   printf("fere\n");
   myfree_movable(zz);
   myfree_movable(yy);
   myfree_movable(xx);
   gsl_rng_free(random_generator);
}


void get_slabs_file()
{
  int ifile, dummy, Ndim=-1, Ndim_task, ntask;
  double box;
  char fname[1000];
  FILE *fd;

  first_slab_file = mymalloc("first_slab_file",NFiles*sizeof(int));
  nslabs_file = mymalloc("nslabs_file", NFiles*sizeof(int));

  if (ThisTask == 0)
  {
      printf("building nslabs_file \n");
      for (ifile=0; ifile<NFiles; ifile++)
      {
          printf("."); fflush(stdout);

          sprintf(fname, "%s/../../fields/fields_%03d/%s_%03d.%d", SimulationDir, SnapNum, "density", SnapNum, ifile);
          fd = my_fopen(fname, "r");

          my_fread(&dummy, sizeof(int), 1, fd);
          my_fread(&ntask, sizeof(int),1,fd);
          my_fread(&box, sizeof(double),1,fd);
          my_fread(&Ndim, sizeof(int), 1, fd);
        // fseek(fd, 3* sizeof(int)+sizeof(double), SEEK_CUR);

          my_fread(&Ndim_task, sizeof(int), 1, fd);
        //  my_fread(&dummy, sizeof(int), 1, fd);

          nslabs_file[ifile] = Ndim_task;

          fclose(fd);
      }
      printf("done \n");

      first_slab_file[0] = 0;
      for (ifile=1; ifile<NFiles; ifile++)
          first_slab_file[ifile] = first_slab_file[ifile-1] + nslabs_file[ifile-1];
  }

  MPI_Bcast(&Ndim, sizeof(int), MPI_BYTE, 0, MPI_COMM_WORLD);

  MPI_Bcast(first_slab_file, sizeof(int)*NFiles, MPI_BYTE, 0, MPI_COMM_WORLD);
  MPI_Bcast(nslabs_file, sizeof(int)*NFiles, MPI_BYTE, 0, MPI_COMM_WORLD);

  if (Ndim != (first_slab_file[NFiles-1]+nslabs_file[NFiles-1]) )
  {
    char buf[500];
    sprintf(buf, "something strange here\n Task=%d, Ndim=%d \n", ThisTask, Ndim);
        terminate(buf);
  }
}


void load_field(Grid *grid)
{
  int i,j,k, igrid, nprev;
  int ifile, nread, end, start, dummy, Ndim, Ndim_task, Ndim2;
  float *tmp, *tmp2;
  double box;
  int ntask;
  char fname[1024];
  FILE *fd;

  get_slabs_file();

  printf("Task=%d nslabs=%d start=%d\n",ThisTask, grid->slabs_per_task[ThisTask],grid->first_slab_of_task[ThisTask]);
  fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);

  for (nprev= 0, ifile=0; ifile < NFiles; ifile++)
  {
    nread = -1;
    start = -1;
    end = -1;

    if (first_slab_file[ifile] >= grid->first_slab_of_task[ThisTask] ||
        first_slab_file[ifile] < (grid->first_slab_of_task[ThisTask] + grid->slabs_per_task[ThisTask]))
    {
      start = MAX(first_slab_file[ifile],grid->first_slab_of_task[ThisTask]) - first_slab_file[ifile];
      end   = MIN(first_slab_file[ifile]+nslabs_file[ifile], grid->first_slab_of_task[ThisTask] + grid->slabs_per_task[ThisTask]) - first_slab_file[ifile];
      nread = end-start;
    }

    if (nread > 0)
    {
      //printf("Task %d reading %d slabs (start=%d) from file %d\n",ThisTask, nread, start, ifile);
      printf("."); fflush(stdout);

      sprintf(fname, "%s/../../fields/fields_%03d/%s_%03d.%d", SimulationDir, SnapNum, "density", SnapNum, ifile);
      fd = my_fopen(fname, "r");

      my_fread(&dummy, sizeof(int), 1, fd);
      my_fread(&ntask, sizeof(int), 1, fd);
      my_fread(&box, sizeof(double), 1, fd);
      my_fread(&Ndim, sizeof(int), 1, fd);
      my_fread(&Ndim_task, sizeof(int), 1, fd);
      my_fread(&dummy, sizeof(int), 1, fd);

      tmp = (float *) mymalloc("tmp", Ndim * Ndim * Ndim_task * sizeof(float));
      my_fread(&dummy, sizeof(int), 1, fd);

      if (start != 0)
         fseek(fd, Ndim * Ndim * (start -1)* sizeof(float), SEEK_CUR);

      my_fread(tmp, sizeof(float), Ndim * Ndim * nread, fd);
      my_fread(&dummy, sizeof(int), 1, fd);
      fclose(fd);

      Ndim2 = (2*(Ndim/2 + 1));
      for(i = 0; i < nread; i++)
        for(j = 0; j < Ndim; j++)
          for(k = 0; k < Ndim; k++)
            {
              igrid = i + first_slab_file[ifile] -  grid->first_slab_of_task[ThisTask];
              igrid = i + nprev;
              grid->data[Ndim * (Ndim * igrid + j) + k] = tmp[Ndim * (Ndim * i + j) + k];
            }

      nprev += nread;

      myfree(tmp);
    }
  }

  myfree(nslabs_file);
  myfree(first_slab_file);
}


void load_dm_single_hdf5(int nr, float nsel)
{
  int j=0, i, k, inside, rep;//, nthis, thisN;
  unsigned long long cum, ngroup, ntot;
  unsigned int nthis, thisN;
  char buf[512];
  FILE *fd;
  float *halo_read;

  hid_t    hdf5_file, attr, hdf5_dataset, hdf5_dataspace_in_file, hdf5_dataspace_in_memory, hdf5_datatype;
  herr_t   ret;
  hsize_t  dimsm[1], dimsm2[2];

  gsl_rng *random_generator;

  gsl_rng_env_setup();
  random_generator = gsl_rng_alloc(gsl_rng_default);

  sprintf(buf, "%s/%s_%03d.%d.hdf5", SimulationDir, SnapshotFileBase, SnapNum, nr);
  if (!(fd = fopen(buf, "r")))
    {
      sprintf(buf, "%s/snapdir_%03d/%s_%03d.%d.hdf5", SimulationDir, SnapNum, SnapshotFileBase, SnapNum, nr);
      fd = my_fopen(buf,"r");
    }
  fclose(fd);

  unsigned long long NumPart[2], NumPartTot[2];
  hdf5_file = H5Fopen(buf, H5F_ACC_RDONLY, H5P_DEFAULT);
  hid_t ghead = H5Gopen(hdf5_file, "/Header", H5P_DEFAULT);

  attr = H5Aopen(ghead, "NumPart_ThisFile", H5P_DEFAULT);
  ret  = H5Aread(attr, H5T_STD_U64LE, &NumPart);
  ret =  H5Aclose(attr);

  attr = H5Aopen(ghead, "NumPart_Total", H5P_DEFAULT);
  ret  = H5Aread(attr, H5T_STD_U64LE, &NumPartTot);
  ret =  H5Aclose(attr);

  H5Gclose(ghead);

  ngroup = NumPart[1];
  ntot   = NumPartTot[1];

  int *flag = (int *) mymalloc_movable(&flag, "flag", ngroup * sizeof(int));
  nthis = ngroup;
#ifdef SUBSAMPLE
  float sample;
  sample = nsel*1./ntot;

  for(nthis = 0, i = 0; i < ngroup; i++)
  {
    if ( gsl_rng_uniform(random_generator) < sample  )
    {
      flag[i] = 1;
      nthis++;
    }
    else
      flag[i] = 0;
  }
#endif

  thisN = ngroup;
  hid_t get_part = H5Gopen(hdf5_file, "/PartType1", H5P_DEFAULT);
  dimsm2[0] = thisN;
  dimsm2[1] = 3;

  if(PP == NULL)
    PP = mymalloc_movable(&PP, "PP", (LocNgroups + nthis) * sizeof(Particle));
  else
    PP = myrealloc_movable(PP, (LocNgroups + nthis) * sizeof(Particle));

  halo_read = mymalloc_movable(&halo_read,"halo_read", 3 * thisN* sizeof(float));

  hdf5_dataset = H5Dopen(get_part, "Coordinates", H5P_DEFAULT);

  hdf5_dataspace_in_memory = H5Screate_simple(2,dimsm2,NULL);
  hdf5_dataspace_in_file   = H5Screate_simple(2,dimsm2,NULL);
  
  H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, hdf5_dataspace_in_memory, hdf5_dataspace_in_file, H5P_DEFAULT, halo_read);
  
  H5Sclose(hdf5_dataspace_in_file);
  H5Sclose(hdf5_dataspace_in_memory);
  H5Dclose(hdf5_dataset);

  for(cum=0, i = 0; i < thisN; i++)
  {
    inside = 1;
    for(k=0; k<3; k++)
    {
      PP[LocNgroups].pos[k] = halo_read[3 * i + k];

      if (PP[LocNgroups].pos[k] >= BoxSize) PP[LocNgroups].pos[k] -= BoxSize;
      if (PP[LocNgroups].pos[k] < 0)  PP[LocNgroups].pos[k] += BoxSize;

      if(PP[LocNgroups].pos[k] < 0 || PP[LocNgroups].pos[k] > BoxSize)
         inside *= 0;
    }

#ifdef SUBSAMPLE
        if( flag[cum] == 1)
#endif
        if( inside == 1)
            LocNgroups++;

       cum++;
  }
  myfree_movable(halo_read);
  assert(cum == ngroup);

#if defined(Z_SPACE) || defined(DIVERGENCE)
  LocNgroups -= nthis;

  halo_read = mymalloc_movable(&halo_read,"halo_read", 3 * thisN* sizeof(float));
  hdf5_dataset = H5Dopen(get_part, "Coordinates", H5P_DEFAULT);

  hdf5_dataspace_in_memory = H5Screate_simple(2,dimsm2,NULL);
  hdf5_dataspace_in_file   = H5Screate_simple(2,dimsm2,NULL);
  
  H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, hdf5_dataspace_in_memory, hdf5_dataspace_in_file, H5P_DEFAULT, halo_read);
  
  H5Sclose(hdf5_dataspace_in_file);
  H5Sclose(hdf5_dataspace_in_memory);
  H5Dclose(hdf5_dataset);

  for(cum=0, i = 0; i < thisN; i++)
  {
    for(k=0; k<3; k++)
    {
#ifdef Z_SPACE
      if(k == 2)
        PP[LocNgroups].pos[k] +=  halo_read[3*i + k] / Hubble / sqrt(ExpFactor);
#endif
#ifdef DIVERGENCE
        PP[LocNgroups].vel[k] =  halo_read[3*i + k] / Hubble / sqrt(ExpFactor);
#endif
        if (PP[LocNgroups].pos[k] >= BoxSize) PP[LocNgroups].pos[k] -= BoxSize;
        if (PP[LocNgroups].pos[k] < 0)  PP[LocNgroups].pos[k] += BoxSize;
    }

#ifdef SUBSAMPLE
    if ( flag[cum] == 1)
#endif
      LocNgroups++;

    cum++;
  }
  myfree_movable(halo_read);
  assert(cum == ngroup);
#endif
  H5Gclose(get_part);

  myfree_movable(flag);
  gsl_rng_free(random_generator);
}


void load_dm_single( int nr, float nsel)
{
  int j=0, i, k, inside, cum, rep, nthis, thisN, num_files, ngroup, totngroup, ntrees, *galsPerTree, dummy;
  char buf[512];
  FILE *fd;
  float *halo_read;
  io_header header;
  gsl_rng *random_generator;

  gsl_rng_env_setup();
  random_generator = gsl_rng_alloc(gsl_rng_default);

  sprintf(buf, "%s/%s_%03d.%d", SimulationDir, SnapshotFileBase, SnapNum, nr);
  if (!(fd = fopen(buf, "r")))
    {
      sprintf(buf, "%s/snapdir_%03d/%s_%03d.%d", SimulationDir, SnapNum, SnapshotFileBase, SnapNum, nr);
      fd = my_fopen(buf,"r");
    }

  my_fread(&dummy, sizeof(int),1,fd);
  my_fread(&header, sizeof(io_header),1,fd);
  my_fread(&dummy, sizeof(int),1,fd);

  ngroup = header.npart[1];
  nthis = ngroup;
  long long ntot;
  float sample;

  ntot = header.npartTotal[1];
  ntot += (((long long) header.npartTotalHighWord[1]) << 32);

  sample = nsel*1./ntot;

  int *flag = (int *) mymalloc_movable(&flag, "flag", ngroup * sizeof(int));
  nthis = ngroup;  

#ifdef SUBSAMPLE
  for(nthis = 0, i = 0; i < ngroup; i++)
  {
    if ( gsl_rng_uniform(random_generator) < sample )
    {
      flag[i] = 1;
      nthis++;
    }
    else
      flag[i] = 0;
  }
#endif

  if(PP == NULL)
    PP = mymalloc_movable(&PP, "PP", (LocNgroups + nthis) * sizeof(Particle));
  else
    PP = myrealloc_movable(PP, (LocNgroups + nthis) * sizeof(Particle));

  my_fread(&dummy, sizeof(int),1,fd);
   for (cum = 0, rep=0; rep < CYCLES; rep++)
   {
      if(rep == (CYCLES - 1))
        thisN = ngroup - (CYCLES - 1) * (ngroup / CYCLES);
      else
        thisN = (ngroup / CYCLES);

     halo_read = mymalloc_movable(&halo_read,"halo_read", 3 * thisN* sizeof(float));
     my_fread(halo_read, thisN * 3 * sizeof(float),1,fd);

     for(i = 0; i < thisN; i++)
     {
       inside = 1;
       for(k=0; k<3; k++)
       {
         PP[LocNgroups].pos[k] = halo_read[3 * i + k];

         if (PP[LocNgroups].pos[k] >= BoxSize) PP[LocNgroups].pos[k] -= BoxSize;
         if (PP[LocNgroups].pos[k] < 0)  PP[LocNgroups].pos[k] += BoxSize;       
  //       if (PP[LocNgroups].pos[k] < 0 || PP[LocNgroups].pos[k] >= BoxSize)
  //         printf("(%f|%f|%f) Box=%f\n",PP[LocNgroups].pos[0],PP[LocNgroups].pos[1],PP[LocNgroups].pos[2], BoxSize);

         if(PP[LocNgroups].pos[k] < 0 || PP[LocNgroups].pos[k] > BoxSize)
            inside *= 0;
       }

#ifdef SUBSAMPLE
        if ( flag[cum] == 1)
#endif
        if ( inside == 1)
        LocNgroups++;

       cum++;
     }
     myfree_movable(halo_read);
   }
  my_fread(&dummy, sizeof(int),1,fd);
  assert(cum == ngroup);
  
#if defined(Z_SPACE) || defined(DIVERGENCE)
   my_fread(&dummy, sizeof(int),1,fd);
   LocNgroups -= nthis;
   for (cum = 0, rep=0; rep < CYCLES; rep++)
   {
      if(rep == (CYCLES - 1))
        thisN = ngroup - (CYCLES - 1) * (ngroup / CYCLES);
      else
        thisN = (ngroup / CYCLES);

     halo_read = mymalloc_movable(&halo_read,"halo_read", 3 * thisN* sizeof(float));
     my_fread(halo_read, thisN * 3 * sizeof(float),1,fd);

     for(i = 0; i < thisN; i++)
     {
        for(k=0; k<3; k++)
        {
#ifdef Z_SPACE
          if(k == 2) 
            PP[LocNgroups].pos[k] +=  halo_read[3*i + k] / Hubble / sqrt(ExpFactor); 
#endif
#ifdef DIVERGENCE
          PP[LocNgroups].vel[k] =  halo_read[3*i + k] / Hubble / sqrt(ExpFactor); 
#endif

          if (PP[LocNgroups].pos[k] >= BoxSize) PP[LocNgroups].pos[k] -= BoxSize;
          if (PP[LocNgroups].pos[k] < 0)  PP[LocNgroups].pos[k] += BoxSize;       
          //if (PP[LocNgroups].pos[k] < 0 || PP[LocNgroups].pos[k] >= BoxSize)
          //  printf("2: (%f|%f|%f) Box=%f\n",PP[LocNgroups].pos[0],PP[LocNgroups].pos[1],PP[LocNgroups].pos[2], BoxSize);
        } 

#ifdef SUBSAMPLE
       if( flag[cum] == 1)
#endif
         LocNgroups++;

       cum++;
     }
     myfree_movable(halo_read);
   }
  my_fread(&dummy, sizeof(int),1,fd);
  assert(cum == ngroup);
#endif
  fclose(fd);
   
  myfree_movable(flag);

  gsl_rng_free(random_generator);

}


void load_halo_catalogue_single_photons(int nr)
{
  hid_t    hdf5_file, attr, hdf5_dataset, hdf5_dataspace_in_file, hdf5_dataspace_in_memory, hdf5_datatype;
  herr_t   ret;
  hsize_t  dimsm[1], dimsm2[2];

  char buf[500], fname[500];
  long long i;
  int j, k;

  float *pos, *vel;

  unsigned long long Ngroups, Nsubgroups, TotNgroups, TotNsubgroups;
  unsigned long long *firstsub, subgroup_offset;
  unsigned int *nsubhalo;
  float        *m_crit200;

  int *ProcKey;

  sprintf(fname, "%s/halo/halocat_%03d/%s_%03d_%d.hdf5", SimulationDir, SnapNum, SnapshotFileBase, SnapNum, nr);

  hdf5_file = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);


  hid_t ghead = H5Gopen(hdf5_file, "/Header", H5P_DEFAULT);

  attr = H5Aopen(ghead, "NumHalos_ThisFile", H5P_DEFAULT);
  ret  = H5Aread(attr, H5T_STD_U64LE, &Ngroups);
  ret =  H5Aclose(attr);

  attr = H5Aopen(ghead, "NumHalos_Total", H5P_DEFAULT);
  ret  = H5Aread(attr, H5T_STD_U64LE, &TotNgroups);
  ret =  H5Aclose(attr);

  attr = H5Aopen(ghead, "NumSubhalos_ThisFile", H5P_DEFAULT);
  ret  = H5Aread(attr, H5T_STD_U64LE, &Nsubgroups);
  ret =  H5Aclose(attr);

  attr = H5Aopen(ghead, "NumSubhalos_Total", H5P_DEFAULT);
  ret  = H5Aread(attr, H5T_STD_U64LE, &TotNsubgroups);
  ret =  H5Aclose(attr);

  attr = H5Aopen(ghead, "SubhalosOffset", H5P_DEFAULT);
  ret  = H5Aread(attr, H5T_STD_U64LE, &subgroup_offset);
  ret =  H5Aclose(attr);

  H5Gclose(ghead);


  if(PP == NULL)
#ifdef SUB_HALOS
    PP = mymalloc_movable(&PP, "PP", (LocNgroups + Nsubgroups) * sizeof(Particle));
#else
    PP = mymalloc_movable(&PP, "PP", (LocNgroups + Ngroups) * sizeof(Particle));
#endif
  else
#ifdef SUB_HALOS
    PP = myrealloc_movable(PP, (LocNgroups + Nsubgroups) * sizeof(Particle));
#else
    PP = myrealloc_movable(PP, (LocNgroups + Ngroups) * sizeof(Particle));
#endif

  /* ------------------------------ */ /* ------------------------------ */
  hid_t get_fof_grp = H5Gopen(hdf5_file, "/FOFHalo", H5P_DEFAULT);
  dimsm[0] = Ngroups;

  /* ------------------------------ */
  hdf5_dataset = H5Dopen(get_fof_grp, "GroupNsubs", H5P_DEFAULT);
  /* number of substructures in FOF group  */
  nsubhalo = (unsigned int *) mymalloc("nsubhalo", Ngroups * sizeof(unsigned int));

  hdf5_dataspace_in_memory = H5Screate_simple(1,dimsm,NULL);
  hdf5_dataspace_in_file   = H5Screate_simple(1,dimsm,NULL);

  H5Dread(hdf5_dataset, H5T_STD_U32LE, hdf5_dataspace_in_memory, hdf5_dataspace_in_file, H5P_DEFAULT, nsubhalo);

  H5Sclose(hdf5_dataspace_in_file);
  H5Sclose(hdf5_dataspace_in_memory);
  H5Dclose(hdf5_dataset);

  /* ------------------------------ */
  /* M_Crit200 */
  hdf5_dataset = H5Dopen(get_fof_grp, "GroupCrit200", H5P_DEFAULT);
  m_crit200 = (float *) mymalloc("m_crit200", Ngroups * sizeof(float));

  hdf5_dataspace_in_memory = H5Screate_simple(1,dimsm,NULL);
  hdf5_dataspace_in_file   = H5Screate_simple(1,dimsm,NULL);

  H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, hdf5_dataspace_in_memory, hdf5_dataspace_in_file, H5P_DEFAULT, m_crit200);

  H5Sclose(hdf5_dataspace_in_file);
  H5Sclose(hdf5_dataspace_in_memory);
  H5Dclose(hdf5_dataset);

  /* ------------------------------ */
  hdf5_dataset = H5Dopen(get_fof_grp, "GroupFirstSub", H5P_DEFAULT);
  /* first substructure in FOF group  */
  firstsub = (unsigned long long *) mymalloc("firstsub", Ngroups * sizeof(unsigned long long));

  hdf5_dataspace_in_memory = H5Screate_simple(1,dimsm,NULL);
  hdf5_dataspace_in_file   = H5Screate_simple(1,dimsm,NULL);

  H5Dread(hdf5_dataset, H5T_STD_U64LE, hdf5_dataspace_in_memory, hdf5_dataspace_in_file, H5P_DEFAULT, firstsub);

  H5Sclose(hdf5_dataspace_in_file);
  H5Sclose(hdf5_dataspace_in_memory);
  H5Dclose(hdf5_dataset);

  ///* ------------------------------ */
  //dimsm2[0] = dimsm[0];
  //dimsm2[1] = 3;

  //hdf5_dataset = H5Dopen(get_fof_grp, "GroupPos", H5P_DEFAULT);
  ///* Pos of substructure  */
  //grpos = (float *) mymalloc("grpos", Ngroups * 3 * sizeof(float));

  //hdf5_dataspace_in_memory = H5Screate_simple(2,dimsm2,NULL);
  //hdf5_dataspace_in_file   = H5Screate_simple(2,dimsm2,NULL);

  //H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, hdf5_dataspace_in_memory, hdf5_dataspace_in_file, H5P_DEFAULT, grpos);

  //H5Sclose(hdf5_dataspace_in_file);
  //H5Sclose(hdf5_dataspace_in_memory);
  //H5Dclose(hdf5_dataset);

  H5Gclose(get_fof_grp);

  /* ------------------------------ */ /* ------------------------------ */
  hid_t get_fof_sub = H5Gopen(hdf5_file, "/Subhalo", H5P_DEFAULT);
  dimsm[0] = Nsubgroups;

  /* ------------------------------ */
  dimsm2[0] = Nsubgroups;
  dimsm2[1] = 3;

  hdf5_dataset = H5Dopen(get_fof_sub, "SubhaloPos", H5P_DEFAULT);
  /* Pos of substructure  */
  pos = (float *) mymalloc("pos", Nsubgroups * 3 * sizeof(float));

  hdf5_dataspace_in_memory = H5Screate_simple(2,dimsm2,NULL);
  hdf5_dataspace_in_file   = H5Screate_simple(2,dimsm2,NULL);

  H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, hdf5_dataspace_in_memory, hdf5_dataspace_in_file, H5P_DEFAULT, pos);

  H5Sclose(hdf5_dataspace_in_file);
  H5Sclose(hdf5_dataspace_in_memory);
  H5Dclose(hdf5_dataset);

  /* ------------------------------ */
  hdf5_dataset = H5Dopen(get_fof_sub, "SubhaloVel", H5P_DEFAULT);
  /* Vel of substructure  */
  vel = (float *) mymalloc("vel", Nsubgroups * 3 * sizeof(float));

  hdf5_dataspace_in_memory = H5Screate_simple(2,dimsm2,NULL);
  hdf5_dataspace_in_file   = H5Screate_simple(2,dimsm2,NULL);

  H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, hdf5_dataspace_in_memory, hdf5_dataspace_in_file, H5P_DEFAULT, vel);

  H5Sclose(hdf5_dataspace_in_file);
  H5Sclose(hdf5_dataspace_in_memory);
  H5Dclose(hdf5_dataset);

  H5Fclose(hdf5_file);

  ProcKey = (int *) mymalloc("ProcKey", Nsubgroups * sizeof(int));

  //-----------------------------------------------
  for(i = 0; i < Nsubgroups; i++)
    ProcKey[i] = -1;

  int proc_flag;
  int cnt_proc, cnt_subs;
  float pchk[3];

  unsigned long long LocNsubs = 0;
  long long LocSubOffset = 0;

  /* ------------------------------ */ /* ------------------------------ */
/*
  LocNsubs = Nsubgroups / NTasks;
  for(i = 0, LocSubOffset = 0; i < ThisTask; i++)
    LocSubOffset += LocNsubs;

  if(ThisTask == NTasks-1)
    LocNsubs = Nsubgroups - LocSubOffset;
*/

  LocNsubs = Nsubgroups;
  LocSubOffset = 0;

  for(i = 0, cnt_proc = 0; i < Ngroups; i++)
    {
      proc_flag = 0;
      if(nsubhalo[i] > 0)
        {
          for(k = 0; k < 3; k++)
          {
            pchk[k] = pos[firstsub[i] * 3 + k] / 1000.;
#if defined(Z_SPACE) && defined(ITERVOL)
            if(k == 2)
              pchk[k] += vel[3 * i + k] / Hubble / ExpFactor; 
#endif
          }

#ifdef ITERVOL
          if(((firstsub[i] >= LocSubOffset) && (firstsub[i] < (LocSubOffset + LocNsubs))))
            {
              if(ShouldWeLoadHalo(pchk[0], pchk[1], pchk[2]))
                proc_flag = 1;
            }
#else
          proc_flag = 1;
#endif
          if(proc_flag > 0)
            {
              ProcKey[firstsub[i]] = 1;
              cnt_proc++;

#ifdef SUB_HALOS
              for(j = 1; j < nsubhalo[i]; j++)
                {
                  ProcKey[firstsub[i]+j] = 1;
                  cnt_proc++;
                }
#endif
            }
        }
    }

  //-----------------------------------------------
  for(i = 0, cnt_subs=0; i < Nsubgroups; i++)
    {
      if(ProcKey[i] > 0)
      {
        for(k = 0; k < 3; k++)
        {
          PP[LocNgroups].pos[k] = pos[3 * i + k] / 1000.; 
#ifdef Z_SPACE
          if(k == 2)
            PP[LocNgroups].pos[k] +=  vel[3 * i + k] / Hubble / ExpFactor; 
#endif
#ifdef DIVERGENCE
          PP[LocNgroups].vel[k] =  vel[3 * i +k] / Hubble / ExpFactor; 
#endif

          if (PP[LocNgroups].pos[k] >= BoxSize) PP[LocNgroups].pos[k] -= BoxSize;
          if (PP[LocNgroups].pos[k] < 0)  PP[LocNgroups].pos[k] += BoxSize;       

          PP[LocNgroups].pos[k] += DomainDepth[k] - DomainCenter[k]; 
        }

        LocNgroups++;
        cnt_subs++;
      }
    }

  assert(cnt_proc == cnt_subs);

  unsigned int LocNgroupsOld;
  LocNgroupsOld = LocNgroups;
  LocNgroups -= cnt_subs;

  for(i = 0; i < Ngroups; i++)
    {
      if((nsubhalo[i] > 0) && (ProcKey[firstsub[i]] > 0))
      {
        PP[LocNgroups].mass = m_crit200[i];
        LocNgroups++;

#ifdef SUB_HALOS
        for(j = 1; j < nsubhalo[i]; j++)
        {
          assert(ProcKey[firstsub[i]+j] > 0);

          PP[LocNgroups].mass = m_crit200[i];
          LocNgroups++;
        }
#endif 
      }
    }

  if(LocNgroups != LocNgroupsOld)
  {
    printf("Task:%d -- loading catalogue. %u vs %u groups with %d \n", ThisTask, LocNgroupsOld, LocNgroups, cnt_subs);
    exit(12121);
  }

  //assert(LocNgroups == LocNgroupsOld);

  myfree(ProcKey);
  myfree(vel);
  myfree(pos);

  myfree(firstsub);
  myfree(m_crit200);
  myfree(nsubhalo);
}


int ShouldWeLoadIter(float x, float y, float z, float IterHalfSize)
{
  float dx, dy, dz;

  dx = periodic_wrap(x - DomainCenter[0]);
  dy = periodic_wrap(y - DomainCenter[1]);
  dz = periodic_wrap(z - DomainCenter[2]);

  if((fabs(dx) < DomainDepth[0] + IterHalfSize) &&
     (fabs(dy) < DomainDepth[1] + IterHalfSize) &&
     (fabs(dz) < DomainDepth[2] + IterHalfSize))
    return 1;

  return 0;
}


int ShouldWeLoadHalo(float x, float y, float z)
{
  float dx, dy, dz;

  dx = periodic_wrap(x - DomainCenter[0]);
  dy = periodic_wrap(y - DomainCenter[1]);
  dz = periodic_wrap(z - DomainCenter[2]);

  if((dx >= -DomainDepth[0]) && (dx < DomainDepth[0]) &&
     (dy >= -DomainDepth[1]) && (dy < DomainDepth[1]) &&
     (dz >= -DomainDepth[2]) && (dz < DomainDepth[2]))
    return 1;

  return 0;
}


void mark_domain_center(int iDiv)
{
  int i_sub, j_sub, k_sub;

  if(ThisTask == 0)
    {
      i_sub =   iDiv / (SubVolDim * SubVolDim);
      j_sub = ((iDiv % (SubVolDim * SubVolDim)) / SubVolDim);
      k_sub = ((iDiv % (SubVolDim * SubVolDim)) % SubVolDim);
      printf("Testing: ijk=(%d %d %d) of SubVolDim %d \n", i_sub, j_sub, k_sub, SubVolDim);

      DomainCenter[0] = ((i_sub + 0.5) * BoxSize) / SubVolDim;
      DomainCenter[1] = ((j_sub + 0.5) * BoxSize) / SubVolDim;
      DomainCenter[2] = ((k_sub + 0.5) * BoxSize) / SubVolDim;

      DomainDepth[0]  = BoxSize / SubVolDim * 0.5;
      DomainDepth[1]  = BoxSize / SubVolDim * 0.5;
      DomainDepth[2]  = BoxSize / SubVolDim * 0.5;

      printf(" - DomainCenter: (%g %g %g) | DomainDepth: (%g %g %g) \n", DomainCenter[0], 
                                                                         DomainCenter[1], 
                                                                         DomainCenter[2],

                                                                         DomainDepth[0], 
                                                                         DomainDepth[1], 
                                                                         DomainDepth[2]);
    }

  MPI_Bcast(DomainCenter, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(DomainDepth , 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}


void mark_Iter_nums(void)
{
  int i_sub, j_sub, k_sub, n, CataVolTot, inr, i;
  double x, y, z;
  float iterhalfsize;

  CataVolDim   = (((peanokey) 1) << (CATAVOL_BITS_PER_DIMENSION));
  CataVolTot   = CataVolDim * CataVolDim * CataVolDim;
  iterhalfsize = 0.0; //BoxSize / CataVolDim ;//* 0.5;
  NumIterstoLoad = 0;

  CataFileMap = mymalloc_movable(&CataFileMap, "CataFileMap", CataVolTot * sizeof(int));
  memset(CataFileMap, -1, CataVolTot * sizeof(int));

  for(inr = 0; inr < CataVolTot; inr++)
    {
      i_sub =   inr / (CataVolDim * CataVolDim);
      j_sub = ((inr % (CataVolDim * CataVolDim)) / CataVolDim);
      k_sub = ((inr % (CataVolDim * CataVolDim)) % CataVolDim);

      x = ((i_sub + 0.5) * BoxSize) / CataVolDim;
      y = ((j_sub + 0.5) * BoxSize) / CataVolDim;
      z = ((k_sub + 0.5) * BoxSize) / CataVolDim;

      if(ShouldWeLoadIter(x, y, z, iterhalfsize) == 1)
        {
          CataFileMap[inr] = inr;
          NumIterstoLoad ++;
        }
    }

  MPI_Bcast(&NumIterstoLoad, 1, MPI_INT, 0, MPI_COMM_WORLD);
  IterstoLoad = mymalloc_movable(&IterstoLoad, "IterstoLoad", sizeof(int) * NumIterstoLoad);

  if(ThisTask == 0)
    {
      printf("  -- %d cat-files should be loaded ...\n\n", NumIterstoLoad);
      for(i = 0, n=0; i < CataVolTot; i++)
        {
          if(CataFileMap[i] >= 0)
            {
              IterstoLoad[n] = CataFileMap[i];
              n++;
            }
        }

      for(i = 0; i < n; i++)
        printf("FileNum Loaded : %d %d \n", i, IterstoLoad[i]);

      assert(n == NumIterstoLoad);
    }
  MPI_Bcast(IterstoLoad, NumIterstoLoad, MPI_INT, 0, MPI_COMM_WORLD);
}


void load_file_tab(void)
{
  char buf[512];
  FILE *fd;
  float dummy;
  int inum, ifile;
  int LastSnapShotNr;

  LastSnapShotNr = 100;

  FileNumTab = mymalloc_movable(&FileNumTab, "FileNumTab", (LastSnapShotNr + 1) * sizeof(int));

  if(ThisTask == 0)
    {
      sprintf(buf, "./output_schedue.txt");
      if(!(fd = fopen(buf, "r")))
      {
        printf("Can not open snap_hash file \n");
        exit(122);
      }
      printf("%s \n", buf);

      while(!feof(fd))
        {
        fscanf(fd, "%d", &inum);
        fscanf(fd, "%d", &ifile);

        fscanf(fd, "%g", &dummy);
        fscanf(fd, "%g", &dummy);
        fscanf(fd, "%g", &dummy);

        FileNumTab[inum] = ifile;
        }
      if(fd != NULL) fclose(fd);
    }
  MPI_Bcast(FileNumTab, (LastSnapShotNr + 1), MPI_INT, 0, MPI_COMM_WORLD);
}


void load_halo_catalogue_single_hdf5(int nr)
{
  int j=0, i, k, rep;
  unsigned int thisN;
  unsigned long long ngroup, totngroup;
  char buf[512], fname[512];
  FILE *fd;

  float mass, *pos, *vel, *m_crit200;

  hid_t    hdf5_file, attr, hdf5_dataset, hdf5_dataspace_in_file, hdf5_dataspace_in_memory, hdf5_datatype;
  herr_t   ret;
  hsize_t  dimsm[1], dimsm2[2];

  gsl_rng *random_generator;

  gsl_rng_env_setup();
  random_generator = gsl_rng_alloc(gsl_rng_default);

  sprintf(buf, "%s/%s_%03d.%d.hdf5", SimulationDir, SnapshotFileBase, SnapNum, nr);
  if (!(fd = fopen(buf, "r")))
    {
      sprintf(buf, "%s/groups_%03d/%s_%03d.%d.hdf5", SimulationDir, SnapNum, SnapshotFileBase, SnapNum, nr);
      fd = my_fopen(buf,"r");
    }
  fclose(fd);

  hdf5_file = H5Fopen(buf, H5F_ACC_RDONLY, H5P_DEFAULT);
  hid_t ghead = H5Gopen(hdf5_file, "/Header", H5P_DEFAULT);

  attr = H5Aopen(ghead, "Ngroups_ThisFile", H5P_DEFAULT);
  ret  = H5Aread(attr, H5T_STD_U64LE, &ngroup);
  ret =  H5Aclose(attr);

  attr = H5Aopen(ghead, "Ngroups_Total", H5P_DEFAULT);
  ret  = H5Aread(attr, H5T_STD_U64LE, &totngroup);
  ret =  H5Aclose(attr);

  H5Gclose(ghead);

  if(PP == NULL)
    PP = mymalloc_movable(&PP, "PP", (LocNgroups + ngroup) * sizeof(Particle));
  else
    PP = myrealloc_movable(PP, (LocNgroups + ngroup) * sizeof(Particle));


  /* ------------------------------ */ /* ------------------------------ */
  hid_t get_fof_grp = H5Gopen(hdf5_file, "/Group", H5P_DEFAULT);
  dimsm[0] = ngroup;

  /* ------------------------------ */
  /* M_Crit200 */
  hdf5_dataset = H5Dopen(get_fof_grp, "Group_M_Crit200", H5P_DEFAULT);
  m_crit200 = mymalloc("m_crit200", ngroup * sizeof(float));

  hdf5_dataspace_in_memory = H5Screate_simple(1,dimsm,NULL);
  hdf5_dataspace_in_file   = H5Screate_simple(1,dimsm,NULL);

  H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, hdf5_dataspace_in_memory, hdf5_dataspace_in_file, H5P_DEFAULT, m_crit200);

  H5Sclose(hdf5_dataspace_in_file);
  H5Sclose(hdf5_dataspace_in_memory);
  H5Dclose(hdf5_dataset);

  /* ------------------------------ */
  dimsm2[0] = dimsm[0];
  dimsm2[1] = 3;

  hdf5_dataset = H5Dopen(get_fof_grp, "GroupPos", H5P_DEFAULT);
  /* Pos of substructure  */
  pos = mymalloc("pos", 3 * ngroup * sizeof(float));

  hdf5_dataspace_in_memory = H5Screate_simple(2,dimsm2,NULL);
  hdf5_dataspace_in_file   = H5Screate_simple(2,dimsm2,NULL);

  H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, hdf5_dataspace_in_memory, hdf5_dataspace_in_file, H5P_DEFAULT, pos);

  H5Sclose(hdf5_dataspace_in_file);
  H5Sclose(hdf5_dataspace_in_memory);
  H5Dclose(hdf5_dataset);

  /* ------------------------------ */
  hdf5_dataset = H5Dopen(get_fof_grp, "GroupVel", H5P_DEFAULT);
  /* Vel of substructure  */
  vel = mymalloc("vel", 3 * ngroup * sizeof(float));

  hdf5_dataspace_in_memory = H5Screate_simple(2,dimsm2,NULL);
  hdf5_dataspace_in_file   = H5Screate_simple(2,dimsm2,NULL);

  H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, hdf5_dataspace_in_memory, hdf5_dataspace_in_file, H5P_DEFAULT, vel);

  H5Sclose(hdf5_dataspace_in_file);
  H5Sclose(hdf5_dataspace_in_memory);
  H5Dclose(hdf5_dataset);

  H5Gclose(get_fof_grp);

  for(i = 0; i < ngroup; i++)
     {
       for(k = 0; k < 3; k++)
       {
         PP[LocNgroups].pos[k] = pos[3 * i + k]; 
#ifdef Z_SPACE
         if (k == 2)
            PP[LocNgroups].pos[k] +=  vel[3 * i + k] / Hubble / ExpFactor; 
#endif
#ifdef DIVERGENCE
         PP[LocNgroups].vel[k] =  vel[3 * i +k] / Hubble / ExpFactor; 
#endif

         if (PP[LocNgroups].pos[k] >= BoxSize) PP[LocNgroups].pos[k] -= BoxSize;
         if (PP[LocNgroups].pos[k] < 0)  PP[LocNgroups].pos[k] += BoxSize;       
       }

       PP[LocNgroups].mass = m_crit200[i];
       LocNgroups++;
     }

  myfree(vel);
  myfree(pos);
  myfree(m_crit200);

  gsl_rng_free(random_generator);
}


void load_halo_catalogue_single(int nr)
{
  int j=0, i, k, rep, thisN, num_files, ngroup, totngroup, ntrees, *galsPerTree;;
  char buf[512], fname[512];
  float mass, *pos, *vel, *m_crit200;
  FILE *fd;
  struct gal_data_struct *gal_read;
  gsl_rng *random_generator;

  gsl_rng_env_setup();
  random_generator = gsl_rng_alloc(gsl_rng_default);

  sprintf(fname, "%s/groups_%03d/%s_%03d.%d", SimulationDir, SnapNum, "subhalo_tab", SnapNum, nr);
  fd = my_fopen(fname, "r");

  my_fread(&ngroup, sizeof(int), 1, fd);

  if(PP == NULL)
    PP = mymalloc_movable(&PP, "PP", (LocNgroups + ngroup) * sizeof(Particle));
  else
    PP = myrealloc_movable(PP, (LocNgroups + ngroup) * sizeof(Particle));


  pos = mymalloc("pos", 3 * ngroup * sizeof(float));
  vel = mymalloc("vel", 3 * ngroup * sizeof(float));
  m_crit200 = mymalloc("m_crit200", ngroup * sizeof(float));

  fseek(fd, 9 *  sizeof(int), SEEK_CUR);
  fseek(fd, ngroup * 4 * sizeof(int), SEEK_CUR);
  my_fread(pos, ngroup, 3 * sizeof(float), fd);
  my_fread(vel, ngroup, 3 * sizeof(float), fd);
  fseek(fd, ngroup * 4 * sizeof(int), SEEK_CUR);
  my_fread(m_crit200, ngroup, sizeof(float), fd);

  for(i = 0; i < ngroup; i++)
     {
       for(k = 0; k < 3; k++)
       {
         PP[LocNgroups].pos[k] = pos[3 * i + k]; 
#ifdef Z_SPACE
         if (k == 2)
            PP[LocNgroups].pos[k] +=  vel[3 * i + k] / Hubble / ExpFactor; 
#endif
#ifdef DIVERGENCE
         PP[LocNgroups].vel[k] =  vel[3 * i +k] / Hubble / ExpFactor; 
#endif

         if (PP[LocNgroups].pos[k] >= BoxSize) PP[LocNgroups].pos[k] -= BoxSize;
         if (PP[LocNgroups].pos[k] < 0)  PP[LocNgroups].pos[k] += BoxSize;       
       }

       PP[LocNgroups].mass = m_crit200[i];
       //  if ( m_crit200[i] < 1)
       //    continue;
       LocNgroups++;
       }

  myfree(m_crit200);
  myfree(vel);
  myfree(pos);
  fclose(fd);

  gsl_rng_free(random_generator);

}


void load_simple_catalogue_single( int nr)
{
  int j=0, i, k, rep, thisN, num_files, ngroup, totngroup, ntrees, *galsPerTree;;
  char buf[512];
  float mass, box;
  FILE *fd;
  struct simple_data_struct *gal_read;
  gsl_rng *random_generator;

  gsl_rng_env_setup();
  random_generator = gsl_rng_alloc(gsl_rng_default);

  sprintf(buf,"%s/%s.dat",SimulationDir,SnapshotFileBase);
  fd = my_fopen(buf,"r");
  my_fread(&box, sizeof(float),1,fd);
  my_fread(&ngroup, sizeof(int),1,fd);
  printf("%d\n",ngroup);
  BoxSize = box;

  if(ngroup < 1) 
    exit(0);

  if(PP == NULL)
    PP = mymalloc_movable(&PP, "PP", (LocNgroups + ngroup) * sizeof(Particle));
  else
    PP = myrealloc_movable(PP, (LocNgroups + ngroup) * sizeof(Particle));

  float var = 0.003*(1./ExpFactor) * 3e5/ Hubble;

  for (rep=0; rep < CYCLES; rep++)
  {
    if(rep == (CYCLES - 1))
      thisN = ngroup - (CYCLES - 1) * (ngroup / CYCLES);
    else
      thisN = (ngroup / CYCLES);

    gal_read = mymalloc_movable(&gal_read,"gal_read", thisN* sizeof(struct simple_data_struct));
    my_fread(gal_read, thisN*sizeof(struct simple_data_struct),1,fd);

    for(i = 0; i < thisN; i++)
    {
      for(k = 0; k < 3; k++)
      {
        PP[LocNgroups].pos[k] = gal_read[i].Pos[k];
#ifdef Z_SPACE
        if(k == 2)
          PP[LocNgroups].pos[k] +=  gal_read[i].Vel[k] / Hubble / ExpFactor; 
#endif
#ifdef PHOTO_SPACE
        if(k == 2)
          PP[LocNgroups].pos[k] +=  gsl_ran_gaussian(random_generator, var);
#endif


#ifdef DIVERGENCE
        PP[LocNgroups].vel[k] =  gal_read[i].Vel[k] / Hubble / ExpFactor; 
#endif

        if (PP[LocNgroups].pos[k] >= BoxSize) PP[LocNgroups].pos[k] -= BoxSize;
        if (PP[LocNgroups].pos[k] < 0)  PP[LocNgroups].pos[k] += BoxSize;       
        if (PP[LocNgroups].pos[k] < 0 || PP[LocNgroups].pos[k] >= BoxSize)
          //printf("(%f|%f|%f) Box=%f, vmax=%f\n", PP[LocNgroups].pos[0],PP[LocNgroups].pos[1],PP[LocNgroups].pos[2], BoxSize, gal_read[i].Vmax);
          printf("(%f|%f|%f) Box=%f, Mvir=%f\n", PP[LocNgroups].pos[0],PP[LocNgroups].pos[1],PP[LocNgroups].pos[2], BoxSize, gal_read[i].Mvir);
      }

      PP[LocNgroups].mass = gal_read[i].Mvir; 
      //PP[LocNgroups].mass = gal_read[i].Vmax; 
      LocNgroups++;
    }
    myfree_movable(gal_read);
  }
  fclose(fd);

  gsl_rng_free(random_generator);
}


void load_simple_catalogue_single2( int nr)
{
  int j=0, i, k, rep, thisN, num_files, ngroup, totngroup, ntrees, *galsPerTree;;
  char buf[512];
  FILE *fd;
  struct simple_data_struct2 *gal_read;

  sprintf(buf,"%s/%s.%d", SimulationDir, SnapshotFileBase, nr);
  fd = my_fopen(buf,"r");
  my_fread(&ntrees, sizeof(int),1,fd);
  galsPerTree = mymalloc("galsPerTree",sizeof(int)*ntrees);
  my_fread(galsPerTree, sizeof(int),ntrees,fd);

  ngroup=0;
  for (j=0; j < ntrees; j++)
    ngroup += galsPerTree[j];

  myfree(galsPerTree);

  if(ngroup < 1) 
    exit(0);

  if(PP == NULL)
    PP = mymalloc_movable(&PP, "PP", (LocNgroups + ngroup) * sizeof(Particle));
  else
    PP = myrealloc_movable(PP, (LocNgroups + ngroup) * sizeof(Particle));

  for (rep=0; rep < CYCLES; rep++)
  {
    if(rep == (CYCLES - 1))
      thisN = ngroup - (CYCLES - 1) * (ngroup / CYCLES);
    else
      thisN = (ngroup / CYCLES);

    gal_read = mymalloc_movable(&gal_read,"gal_read", thisN* sizeof(struct simple_data_struct2));
    my_fread(gal_read, thisN*sizeof(struct simple_data_struct2),1,fd);

    for(i = 0; i < thisN; i++)
    {
      //PP[LocNgroups].pos[0] = gal_read[i].Pos[0] - 375.;
      //PP[LocNgroups].pos[1] = gal_read[i].Pos[1] - 750.;
      //PP[LocNgroups].pos[2] = gal_read[i].Pos[2] - 375.;
     
      PP[LocNgroups].pos[0] = gal_read[i].Pos[0] - 0.0 ;
      PP[LocNgroups].pos[1] = gal_read[i].Pos[1] - 875.;
      PP[LocNgroups].pos[2] = gal_read[i].Pos[2] - 0.0 ;

      if(PP[LocNgroups].pos[0]  > 500.) PP[LocNgroups].pos[0] -= 1000.;
      if(PP[LocNgroups].pos[2]  > 500.) PP[LocNgroups].pos[2] -= 1000.;
      if(PP[LocNgroups].pos[1]  <-500.) PP[LocNgroups].pos[1] += 1000.;

      for(k = 0; k < 3; k++)
      {
        //PP[LocNgroups].pos[k] = gal_read[i].Pos[k];

        if(PP[LocNgroups].pos[k] >= BoxSize) PP[LocNgroups].pos[k] -= BoxSize;
        if(PP[LocNgroups].pos[k] < 0)        PP[LocNgroups].pos[k] += BoxSize;       
        if(PP[LocNgroups].pos[k] < 0 || PP[LocNgroups].pos[k] >= BoxSize)
          printf("(%f|%f|%f) Box=%f, Mvir=%f type=%d \n", PP[LocNgroups].pos[0], 
                                                          PP[LocNgroups].pos[1], 
                                                          PP[LocNgroups].pos[2], 
                                                          BoxSize, gal_read[i].M_Crit200_atInFall, 
                                                          gal_read[i].Type);
      }

      PP[LocNgroups].mass = gal_read[i].M_Crit200_atInFall; 
      LocNgroups++;
    }
    myfree_movable(gal_read);
  }
  fclose(fd);
}


void load_catalogue_single(int nr)
{
  int j=0, i, k, rep, thisN, num_files, ngroup, totngroup, ntrees, *galsPerTree;;
  char buf[512];
  float mass;
  FILE *fd;
  struct gal_data_struct *gal_read;
  gsl_rng *random_generator;

  gsl_rng_env_setup();
  random_generator = gsl_rng_alloc(gsl_rng_default);

  sprintf(buf,"%s/%s_%d", SimulationDir, SnapshotFileBase, nr);

  fd = my_fopen(buf,"r");
  my_fread(&ntrees, sizeof(int),1,fd);
  my_fread(&ngroup, sizeof(int),1,fd);
  galsPerTree = mymalloc("galsPerTree",sizeof(int)*ntrees);
  my_fread(galsPerTree, sizeof(int),ntrees,fd);
  myfree(galsPerTree);
//  printf("%s %d\n",buf, ngroup);

  if(ngroup < 1) 
    exit(0);

  if(PP == NULL)
    PP = mymalloc_movable(&PP, "PP", (LocNgroups + ngroup) * sizeof(Particle));
  else
    PP = myrealloc_movable(PP, (LocNgroups + ngroup) * sizeof(Particle));

  float var = 0.003*(1./ExpFactor) * 3e5/ Hubble;

  for(rep=0; rep < CYCLES; rep++)
  {
    if(rep == (CYCLES - 1))
      thisN = ngroup - (CYCLES - 1) * (ngroup / CYCLES);
    else
      thisN = (ngroup / CYCLES);

    gal_read = mymalloc_movable(&gal_read,"halo_read", thisN* sizeof(struct gal_data_struct));
    my_fread(gal_read, thisN*sizeof(struct gal_data_struct),1,fd);

    for(i = 0; i < thisN; i++)
    {
      for(k = 0; k < 3; k++)
      {
        PP[LocNgroups].pos[k] = gal_read[i].Pos[k];

#ifdef Z_SPACE
        if (k == 2)
          PP[LocNgroups].pos[k] +=  gal_read[i].Vel[k] / Hubble / ExpFactor; 
#endif

#ifdef PHOTO_SPACE
        if (k == 2)
          PP[LocNgroups].pos[k] +=  gsl_ran_gaussian(random_generator, var);
#endif


#ifdef DIVERGENCE
        PP[LocNgroups].vel[k] =  gal_read[i].Vel[k] / Hubble / ExpFactor; 
#endif

        if (PP[LocNgroups].pos[k] >= BoxSize) PP[LocNgroups].pos[k] -= BoxSize;
        if (PP[LocNgroups].pos[k] < 0)  PP[LocNgroups].pos[k] += BoxSize;       
        if (PP[LocNgroups].pos[k] < 0 || PP[LocNgroups].pos[k] >= BoxSize)
          printf("(%f|%f|%f) Box=%f, mass=%f\n",PP[LocNgroups].pos[0],PP[LocNgroups].pos[1],PP[LocNgroups].pos[2], BoxSize, gal_read[i].CentralMvir);
      }

#ifdef GAL_SHUFFLE
      PP[LocNgroups].HaloMass = gal_read[i].CentralMvir;
      PP[LocNgroups].Type = gal_read[i].Type;
#endif

      if(FileFormat == 1)
      {
        PP[LocNgroups].mass = gal_read[i].BulgeMass + gal_read[i].DiskMass; 
        if (gal_read[i].CentralMvir < 15)// && halo_read[i].Type == 0)
          continue;
      }

      if(FileFormat == 2)
      {
        PP[LocNgroups].mass = gal_read[i].SFR;
        if (gal_read[i].CentralMvir < 15)// && halo_read[i].Type != 0)
          continue;
      }

      if(FileFormat == 3)
      {
        PP[LocNgroups].mass = gal_read[i].Mvir; 
        if ( gal_read[i].CentralMvir < 15 || gal_read[i].Type != 0)
          continue;
      }

      LocNgroups++;
    }
    myfree_movable(gal_read);
  }
  fclose(fd);

  gsl_rng_free(random_generator);
}


void dump_particles(int np)
{
  int i, j;
  FILE *fd;
  char buf[500];

  sprintf(buf,"particle_%d", NTasks);
 
  for(i = 0; i < NTasks; i++ )
  {
    if (ThisTask == i)
    {
      fd = fopen(buf, "a");
      for (j=0; j < np; j++) 
        fprintf(fd, "%f %f %d\n", PP[j].pos[0], PP[j].mass, ThisTask);
    }
    MPI_Barrier(MPI_COMM_WORLD);
   }
 
   if(ThisTask == 0)
     fclose(fd);
}


double periodic_wrap(double x)
{
  if(x > 0.5 * BoxSize)
    x -= BoxSize;
  if(x < -0.5 * BoxSize)
    x += BoxSize;

  return x;
}

