#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <signal.h>
#include <assert.h>

#include "allvars.h"
#include "proto.h"


void load_catalogue()
{
  int nr, num_files, i, ig, readgroups;
  long long buffer;

  num_files = NFiles;
  readgroups = floor(NTask/NumFilesWrittenInParallel);

  if (ThisTask == 0)
    printf("loading catalogue. %d files in %d groups \n",num_files,readgroups);

  LocNgroups = 0;

  for(ig = 0; ig < readgroups; ig++)
  {
    if(ThisTask == 0)
      printf("reading group = %d \n",ig);
    for(nr = 0; nr < num_files; nr++)
      if((nr % NTask) == ThisTask && (ThisTask % readgroups == ig))
      {
		load_halo_catalogue_single(nr);
        //printf("."); fflush(stdout);
      }
     MPI_Barrier(MPI_COMM_WORLD);
  }

  buffer = (long long) LocNgroups;
  MPI_Allreduce(&buffer, &TotNumHalo, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

  if (ThisTask == 0)
    printf("Number of objects Total %lld \n", TotNumHalo);
}

void load_snapshot()
{
  int nr, num_files, it, ig, in, fperread;
  long long buffer, LocNcumu, TotNcumu;
  double local_total = 0.0, total = 0.0;

  num_files = NFiles;
  fperread = floor(num_files/NTask/NumFilesWrittenInParallel);

  if (ThisTask == 0)
    printf("loading snapshot. %d files and %d perread \n",num_files,fperread);

  LocNcumu = 0;
  memset(Grid, 0, maxfftsize * sizeof(fftw_real));

  for(in = 0; in < NumFilesWrittenInParallel; in++)
  {
  LocNgroups = 0;
  for(ig = 0; ig < fperread; ig++)
  {
    nr = in * fperread + ig + ThisTask * (num_files/NTask);
    if (ThisTask == 0)
        printf("file loading this time %d on task %d\n", nr, ThisTask);
    load_halo_catalogue_single(nr);
    MPI_Barrier(MPI_COMM_WORLD);
  }

    TotNumHalo = 0;
    LocNcumu += LocNgroups;
    buffer = (long long) LocNgroups;
    MPI_Allreduce(&buffer, &TotNumHalo, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    if (ThisTask == 0)
        printf("Number of particles load this time %lld \n", TotNumHalo);

    TotNcumu = 0;
    buffer = (long long) LocNcumu;
    MPI_Allreduce(&buffer, &TotNcumu, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    if (ThisTask == 0)
        printf("Number of particles loaded until now  %lld \n", TotNcumu);

    report_memory_usage(&HighMark_run, "MAPPing");
    density(PP, Grid, LocNgroups, 1, 0);

    density_statistics(Grid);
    myfree(PP);
    for(it = 0; it < NTask; it++)
      if(it == ThisTask)
        PP = NULL;
  }

  TotNumHalo = 0;
  buffer = (long long) LocNgroups;
  MPI_Allreduce(&buffer, &TotNumHalo, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

  if (ThisTask == 0)
    printf("Number of particles Total %lld \n", TotNumHalo);

  density_statistics(Grid);
  save_grids(0);
}


void load_halo_catalogue_single(int FileNr)
{
  int j=0, i, k, count;
  char buf[512], fname[512];
  float *pos, *vel;
  FILE *fd_halos;
  io_header header;

#ifdef DMP
  sprintf(buf, "%s/snapshots/snapdir_%03d/%s_%03d.%d", SimulationDir, SnapshotNum, "snapshot", SnapshotNum, FileNr);
  if(!(fd_halos = fopen(buf, "r")))
    {
      sprintf(buf, "can't read file `%s`\n", buf);
      terminate(buf);
    }

  int Npart, dummy;
  long long NpartTotal;
  int cum, rep, thisN;

  my_fread(&dummy, sizeof(int),1,fd_halos);
  my_fread(&header, sizeof(io_header),1,fd_halos);
  my_fread(&dummy, sizeof(int),1,fd_halos);

  Npart = header.npart[1];

  NpartTotal  = header.npartTotal[1];
  NpartTotal += (((long long) header.npartTotalHighWord[1]) << 32);

  if(ThisTask == 0 && FileNr == 0)
	printf("Total particles in the whole snapshot %lld | local Npart %d \n", NpartTotal, Npart);

  if(PP == NULL)
    PP = mymalloc_movable(&PP, "PP", (LocNgroups + Npart) * sizeof(Particle));
  else
    PP = myrealloc_movable(PP, (LocNgroups + Npart) * sizeof(Particle));

  my_fread(&dummy, sizeof(int), 1, fd_halos);
  for (cum = 0, rep=0; rep < CYCLES; rep++)
  {
     if(rep == (CYCLES - 1))
       thisN = Npart- (CYCLES - 1) * (Npart / CYCLES);
     else
       thisN = (Npart / CYCLES);

     pos = (float *) mymalloc("pos", 3 * thisN * sizeof(float));
     my_fread(pos, sizeof(float), 3 * thisN, fd_halos);

     for(i = 0; i < thisN; i++)
     {
       for(k = 0; k < 3; k++)
         PP[LocNgroups].pos[k] = pos[3 * i + k];

       LocNgroups++;
       cum ++;
     }
     myfree_movable(pos);
  }

  my_fread(&dummy, sizeof(int), 1, fd_halos);
  assert(cum == Npart);

  LocNgroups -= Npart;
  my_fread(&dummy, sizeof(int), 1, fd_halos);
  for (cum = 0, rep=0; rep < CYCLES; rep++)
  {
     if(rep == (CYCLES - 1))
       thisN = Npart- (CYCLES - 1) * (Npart / CYCLES);
     else
       thisN = (Npart / CYCLES);

     vel = (float *) mymalloc("vel", 3 * thisN * sizeof(float));
     my_fread(vel, sizeof(float), 3 * thisN, fd_halos);

     for(i = 0; i < thisN; i++)
     {
       for(k = 0; k < 3; k++)
         PP[LocNgroups].vel[k] = vel[3 * i + k];

       LocNgroups++;
       cum ++;
     }
     myfree_movable(vel);
  }

  my_fread(&dummy, sizeof(int), 1, fd_halos);
  assert(cum == Npart);
  fclose(fd_halos);

#else
  sprintf(buf, "%s/groups_tab/groups_%03d/%s_%03d.%d", SimulationDir, SnapshotNum, "subhalo_tab", SnapshotNum, FileNr);
  if(!(fd_halos = fopen(buf, "r")))
    {
      sprintf(buf, "can't read file `%s`\n", buf);
      terminate(buf);
    }

  int nfiles, Nsubgroups, Nids, Ngroups;
  long long TotNgroups, TotNids, TotNsubgroups;

  my_fread(&Ngroups, sizeof(int), 1, fd_halos);
  my_fread(&TotNgroups, sizeof(long long), 1, fd_halos);
  my_fread(&Nids, sizeof(int), 1, fd_halos);
  my_fread(&TotNids, sizeof(long long), 1, fd_halos);
  my_fread(&nfiles, sizeof(int), 1, fd_halos);
  my_fread(&Nsubgroups, sizeof(int), 1, fd_halos);
  my_fread(&TotNsubgroups, sizeof(long long), 1, fd_halos);

#ifndef SUBHALO
  if(PP == NULL)
    PP = mymalloc_movable(&PP, "PP", (LocNgroups + Ngroups) * sizeof(Particle));
  else
    PP = myrealloc_movable(PP, (LocNgroups + Ngroups) * sizeof(Particle));

  if(ThisTask == 0 && FileNr == 0)
	printf("Total Halos in the whole snapshot %llu \n", TotNgroups);

  pos = (float *) mymalloc("pos", 3 * Ngroups * sizeof(float));
  vel = (float *) mymalloc("vel", 3 * Ngroups * sizeof(float));

  /* len and offset into id-list */
  fseek(fd_halos, 2*Ngroups * sizeof(int), SEEK_CUR);

  /* group number */
  fseek(fd_halos, Ngroups * sizeof(long long), SEEK_CUR);

  /* CM */
  fseek(fd_halos, 3 * Ngroups * sizeof(float), SEEK_CUR);

  my_fread(vel, Ngroups, 3 * sizeof(float), fd_halos);
  my_fread(pos, Ngroups, 3 * sizeof(float), fd_halos);

  count = Ngroups;
#else
  if(PP == NULL)
    PP = mymalloc_movable(&PP, "PP", (LocNgroups + Nsubgroups) * sizeof(Particle));
  else
    PP = myrealloc_movable(PP, (LocNgroups + Nsubgroups) * sizeof(Particle));

  if(ThisTask == 0 && FileNr == 0)
	printf("Total Subhalos in the whole snapshot %lld \n", TotNsubgroups);

  pos = (float *) mymalloc("pos", 3 * Nsubgroups * sizeof(float));
  vel = (float *) mymalloc("vel", 3 * Nsubgroups * sizeof(float));

  /* len and offset into id-list */
  fseek(fd_halos, 2*Ngroups * sizeof(int), SEEK_CUR);

  /* group number */
  fseek(fd_halos, Ngroups * sizeof(long long), SEEK_CUR);

  /* CM, vel, location, 3 masses, 4 vel disp, nsub, firstsub */
  fseek(fd_halos, 18* Ngroups * sizeof(float), SEEK_CUR);

  /*---------------------------------*/
  /* Len and offset of substructure  */
  fseek(fd_halos, 2 * Nsubgroups * sizeof(int), SEEK_CUR);

  /* GrNr and SubNr of substructure  */
  fseek(fd_halos, 2 * Nsubgroups * sizeof(long long), SEEK_CUR);

  my_fread(pos, Nsubgroups, 3 * sizeof(float), fd_halos);
  my_fread(vel, Nsubgroups, 3 * sizeof(float), fd_halos);

  count = Nsubgroups;
#endif
  fclose(fd_halos);

  for(i = 0; i < count; i++)
  {
       for(k = 0; k < 3; k++)
       {
         PP[LocNgroups].pos[k] = pos[3 * i + k];
         PP[LocNgroups].vel[k] = vel[3 * i + k];
  
         if (PP[LocNgroups].pos[k] >= BoxSize) PP[LocNgroups].pos[k] -= BoxSize;
         if (PP[LocNgroups].pos[k] < 0)  PP[LocNgroups].pos[k] += BoxSize;
       }
  
       LocNgroups++;
  }

  //if(ThisTask == 0)
  //  printf("Halos read in file %d -- %d \n", FileNr, Ngroups);

  myfree(vel);
  myfree(pos);
#endif
}


void load_grid_field(int dim)
{
  int i,j,k, igrid, nprev;
  int ifile, nread, end, start, dummy, Ndim, Ndim_task, Ndim2, NFilesR;
  float *tmp, *tmp2;
  double box, fac;
  char fname[1024];
  FILE *fd;

  printf("Task=%d nslabs=%d start=%d\n",ThisTask, slabs_per_task[ThisTask],first_slab_of_task[ThisTask]);
  //fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);

  for (nprev= 0, ifile=0; ifile<NFiles; ifile++)
  {
       nread = -1;
       start = -1;
       end = -1;

       if (first_slab_file[ifile] >= first_slab_of_task[ThisTask] ||  
           first_slab_file[ifile] < (first_slab_of_task[ThisTask]+slabs_per_task[ThisTask]))
       {
			start = MAX(first_slab_file[ifile],first_slab_of_task[ThisTask]) - first_slab_file[ifile];
			end   = MIN(first_slab_file[ifile]+nslabs_file[ifile], first_slab_of_task[ThisTask] + slabs_per_task[ThisTask]) - first_slab_file[ifile];
			nread = end-start;
       }

       if (nread > 0)
       {
            //printf("Task %d reading %d slabs (start=%d) from file %d\n",ThisTask, nread, start, ifile);
           // printf("."); fflush(stdout);
          
			sprintf(fname, "%s/halo_fields_%03d/%s_density_%03d.%d", OutputDir, SnapshotNum, FileNameBase, SnapshotNum, ifile);
            fd = my_fopen(fname, "r");

            my_fread(&NTask, sizeof(int), 1, fd);
            my_fread(&NFilesR, sizeof(int), 1, fd);
            my_fread(&BoxSize, sizeof(double), 1, fd);
            my_fread(&Ndim, sizeof(int), 1, fd);
            my_fread(&Ndim_task, sizeof(int), 1, fd);
            my_fread(&dummy, sizeof(int), 1, fd);

            tmp = (float *) mymalloc("tmp", Ndim * Ndim * Ndim_task * sizeof(float));
            if (start != 0) 
               fseek(fd, Ndim * Ndim * (start -1)* sizeof(float), SEEK_CUR);
            my_fread(tmp, sizeof(float), Ndim * Ndim * nread, fd);
            fclose(fd);

#if defined (SMOOTH_HALOFIELD) && defined (SMOOTH_HALOVELFIELD)
			sprintf(fname, "%s/halo_fields_%03d/%s_velocity_%d_%03d.%d", OutputDir, SnapshotNum, FileNameBase, dim, SnapshotNum, ifile);
            fd = my_fopen(fname, "r");

            fseek(fd, 5 * sizeof(int) + sizeof(double), SEEK_CUR);

            tmp2 = (float *) mymalloc("tmp2", Ndim * Ndim * Ndim_task * sizeof(float));
            if (start != 0) 
               fseek(fd, Ndim * Ndim * (start -1)* sizeof(float), SEEK_CUR);
            my_fread(tmp2, sizeof(float), Ndim * Ndim * nread, fd);
            fclose(fd);
#endif
            Ndim2 = (2*(Ndim/2 + 1));
            for(i = 0; i < nread; i++)
                for(j = 0; j < Ndim; j++)
                    for(k = 0; k < Ndim; k++)
                    {
                        igrid = i + first_slab_file[ifile] -  first_slab_of_task[ThisTask];
                        igrid = i + nprev;
#if defined (SMOOTH_HALOFIELD) && defined (SMOOTH_HALOVELFIELD)
                        if (tmp[Ndim * (Ndim * i + j) + k] == 0)
                          Grid[((large_array_offset) Ndim2) * (Ndim * igrid + j) + k] = 0.0;
                        else
						  Grid[((large_array_offset) Ndim2) * (Ndim * igrid + j) + k] = tmp2[Ndim * (Ndim * i + j) + k];
#else
                        Grid[((large_array_offset) Ndim2) * (Ndim * igrid + j) + k] = tmp[Ndim * (Ndim * i + j) + k];
#endif
                    }

	    nprev += nread;

#if defined (SMOOTH_HALOFIELD) && defined (SMOOTH_HALOVELFIELD)
		myfree(tmp2);
#endif
		myfree(tmp);
       }
  }


#if defined (SMOOTH_HALOFIELD) || defined (LINEAR_VELFIELD)
  tmp = (float *) mymalloc("tmp", Ndim * Ndim * slabs_per_task[ThisTask] * sizeof(float));

  for(i = 0; i < slabs_per_task[ThisTask]; i++)
   for(j = 0; j < Ndim; j++)
      for(k = 0; k < Ndim; k++)
        tmp[Ndim * (Ndim * i + j) + k] = Grid[((large_array_offset) Ndim2) * (Ndim * i + j)+k];

  TotNumHalo = density_statistics(tmp);

#ifdef LINEAR_VELFIELD
  fac = TotNumHalo / CUB((double)Ndim);

#ifdef LOGDENS
  double log_mean;

  log_mean = log_density(tmp);
#endif

for(i = 0; i < slabs_per_task[ThisTask]; i++)
 for(j = 0; j < Ndim; j++)
    for(k = 0; k < Ndim; k++)
    {
#ifdef LOGDENS
        Grid[((large_array_offset) Ndim2) * (Ndim * i + j)+k] = log(tmp[Ndim * (Ndim * i + j) + k] / fac) - log_mean;
#else
      Grid[((large_array_offset) Ndim2) * (Ndim * i + j)+k] = tmp[Ndim * (Ndim * i + j) + k] / fac - 1.0;
#endif
    }
#endif
  myfree(tmp);
#endif

}

void get_slabs_file()
{
  int ifile, Ndim_task, NFilesR;
  double box;
  char fname[1000], buf[500];
  FILE *fd;
  int *first_slab_read, first_slab_local; 

  first_slab_file = mymalloc("first_slab_file",NFiles*sizeof(int));
  nslabs_file = mymalloc("nslabs_file", NFiles*sizeof(int));

  first_slab_read = mymalloc("first_slab_read",NFiles*sizeof(int));

  if (ThisTask == 0)
  {
      printf("building nslabs_file \n");
      for (ifile=0; ifile<NFiles; ifile++)
      {
          //printf("."); fflush(stdout);

          sprintf(fname, "%s/halo_fields_%03d/%s_density_%03d.%d", OutputDir, SnapshotNum, FileNameBase, SnapshotNum, ifile);
          fd = my_fopen(fname, "r");

          my_fread(&NTask, sizeof(int), 1, fd);
          my_fread(&NFilesR, sizeof(int), 1, fd);
          my_fread(&BoxSize, sizeof(double), 1, fd);
          my_fread(&Ndim, sizeof(int), 1, fd);
          my_fread(&Ndim_task, sizeof(int), 1, fd);
          my_fread(&first_slab_local, sizeof(int), 1, fd);

          nslabs_file[ifile] = Ndim_task;
          first_slab_read[ifile] = first_slab_local;

          fclose(fd);
      }
      printf("done \n");

      first_slab_file[0] = 0;
      for (ifile=1; ifile<NFiles; ifile++)
	  {
		first_slab_file[ifile] = first_slab_file[ifile-1] + nslabs_file[ifile-1];

		if(first_slab_file[ifile] != first_slab_read[ifile])
		{
		  sprintf(buf, "something strange here\n Task=%d, first_slab_file %d and first_slab_read\n", ThisTask, first_slab_file[ifile], first_slab_read[ifile]);
		  terminate(buf);
		}
	  }
  }

  MPI_Bcast(first_slab_file, sizeof(int)*NFiles, MPI_BYTE, 0, MPI_COMM_WORLD);
  MPI_Bcast(nslabs_file, sizeof(int)*NFiles, MPI_BYTE, 0, MPI_COMM_WORLD);

  if (Ndim != (first_slab_file[NFiles-1]+nslabs_file[NFiles-1]) )
  {
    sprintf(buf, "something strange here\n Task=%d, Ndim=%d \n", ThisTask, Ndim);
	terminate(buf);
  }

  myfree(first_slab_read);
}


void save_grids(int dim)
{
  int i, j, k, Ndim2;
  char buf[1000];
  FILE *fd;

  float *tmp = (float *) mymalloc("tmp", Ndim * Ndim * slabs_per_task[ThisTask] * sizeof(float));

  Ndim2 = (2*(Ndim/2 + 1)); 

  if (ThisTask == 0)
  {
    sprintf(buf, "%s/halo_fields_%03d/", OutputDir, SnapshotNum);
    mkdir(buf, 02755);
#if defined (SMOOTH_HALOFIELD) || defined (LINEAR_VELFIELD)
    sprintf(buf, "%s/halo_fields_%03d/%d/", OutputDir, SnapshotNum, rep);
    mkdir(buf, 02755);
#endif
	printf("Saving grids to %s\n", buf);
  }
 
  MPI_Barrier(MPI_COMM_WORLD); 
  for(i = 0; i < slabs_per_task[ThisTask]; i++)
    for(j = 0; j < Ndim; j++)
      for(k = 0; k < Ndim; k++)
#ifdef GEN_HALOFIELD
        tmp[Ndim * (Ndim * i + j) + k] = Grid[Ndim * (Ndim * i + j) + k];
#else
#ifndef LINEAR_VELFIELD
        tmp[Ndim * (Ndim * i + j) + k] = sGrid[((large_array_offset) Ndim2) * (Ndim * i + j) + k];
#else
        tmp[Ndim * (Ndim * i + j) + k] = Vel[((large_array_offset) Ndim2) * (Ndim * i + j) + k];
#endif
#endif

  if (ThisTask == 0)
      printf("after data assign \n");

#ifdef GEN_HALOFIELD
  if(dim == 0)
	sprintf(buf, "%s/halo_fields_%03d/%s_density_%03d.%d", OutputDir, SnapshotNum, FileNameBase, SnapshotNum, ThisTask);
  if(dim > 0)
	sprintf(buf, "%s/halo_fields_%03d/%s_velocity_%d_%03d.%d", OutputDir, SnapshotNum, FileNameBase, dim-1, SnapshotNum, ThisTask);
#endif

#ifdef SMOOTH_HALOFIELD
  sprintf(buf, "%s/halo_fields_%03d/%d/%s_density_smooth_%ld_%03d_%d.%d", OutputDir, SnapshotNum, rep, FileNameBase, Ndim, SnapshotNum, rep, ThisTask);
#endif
#if defined (SMOOTH_HALOFIELD) && defined (SMOOTH_HALOVELFIELD)
  sprintf(buf, "%s/halo_fields_%03d/%d/%s_velocity_smooth_%d_%ld_%03d_%d.%d", OutputDir, SnapshotNum, rep, FileNameBase, dim, Ndim, SnapshotNum, rep, ThisTask);
#endif

#ifdef LINEAR_VELFIELD
#ifdef LOGDENS
  sprintf(buf, "%s/halo_fields_%03d/%d/%s_rec_vel_log_dens_%d_%ld_%03d_%d.%d", OutputDir, SnapshotNum, rep, FileNameBase, dim, Ndim, SnapshotNum, rep, ThisTask);
#else
  sprintf(buf, "%s/halo_fields_%03d/%d/%s_rec_vel_%d_%ld_%03d_%d.%d", OutputDir, SnapshotNum, rep, FileNameBase, dim, Ndim, SnapshotNum, rep, ThisTask);
#endif
#endif

  if(!(fd = fopen(buf, "w")))
    {
      printf("can't open file `%s'\n", buf);
      exit(8);
    }

  fwrite(&NTask, 1, sizeof(int), fd);
  fwrite(&NFiles, 1, sizeof(int), fd);
  fwrite(&BoxSize, 1, sizeof(double), fd);
  fwrite(&Ndim, 1, sizeof(int), fd);
  fwrite(&slabs_per_task[ThisTask], 1, sizeof(int), fd);
  fwrite(&first_slab_of_task[ThisTask], 1, sizeof(int), fd);
  fwrite(tmp, sizeof(float), Ndim * Ndim * slabs_per_task[ThisTask], fd);
  fclose(fd);

  myfree(tmp);

}

