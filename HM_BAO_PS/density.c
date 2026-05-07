#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <assert.h>
#include <mpi.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_histogram2d.h>

#include <srfftw_mpi.h>

#include "allvars.h"
#include "macros.h"
#include "utils.h"
//#include "gadget.h"
#include "mymalloc.h"
#include "density.h"

//extern int max_size, NProcRead;
//extern double DiluteFactor, invCellSize;
//extern double rms, mean, sd, skew, kurtosis;
//extern time_t t0;

//int local_x_start, local_nx, local_y_start, local_ny, total_local_size;
//int *slab_to_task, *slab_to_task_local, *slabs_per_task, *first_slab_of_task;
//double to_slab_fac;

extern time_t t0;
extern double *MyBoxSize;
extern int ThisTask, NTasks, NumFilesWrittenInParallel;

#define dbg1 0
#define dbg2 0

//grid* density(particle *P, long FirstPart, long NumPart, double *grid)
Grid init_density(int Ngrid, int foldnum)
{
  long i, j;
  Grid grid;

  grid.Ngrid = Ngrid;

  grid_init(&grid, foldnum);

  DBG(" local_x_start %d", dbg1, t0, grid.local_x_start);

//  for(i = 0; i < grid.Ngrid; i++)
//    DBG0("task=%d slab=%d", dbg1, t0, ThisTask, grid.slab_to_task[i], i);


  MSG0(" Allocating %.2f Mb to hold density field (%d,%d,%d)", t0, ThisTask, sizeof(MyFloat) *
       SQR(grid.Ngrid) * grid.local_nx / SQR(1024.), grid.Ngrid, grid.Ngrid, grid.local_nx);

  if(grid.local_nx > 0)
  {
    grid.data = mymalloc_movable(&grid.data,"grid.data", SQR(grid.Ngrid) * grid.local_nx * sizeof(MyFloat));
    //grid.data = mymalloc("grid.data", SQR(grid.Ngrid) * grid.local_nx * sizeof(MyFloat));
#ifdef CENTER_OF_MASS
    grid.pos = mymalloc("grid.pos", 3 * SQR(grid.Ngrid) * grid.local_nx * sizeof(MyFloat));
#endif
  }

  for(i = 0; i < grid.Ngrid * grid.Ngrid * grid.local_nx; i++)
  {
    grid.data[i] = (MyFloat) 0.0;

#ifdef CENTER_OF_MASS
    for(j=0; j<3; j++)
      grid.pos[3*i+j] = (MyFloat) 0.0;
#endif
  }
  return (grid);
}

void grid_init(Grid * grid, int foldnum)
{
  int ii;
  int local_nx, local_x_start, local_ny, local_y_start, total_local_size;
  int *slab_to_task_local;
  rfftwnd_mpi_plan fft_forward_plan, fft_inverse_plan;

  DBG0(" Creating plan for %d cells", dbg1, t0, ThisTask, (int) grid->Ngrid);

  grid->to_slab_fac = grid->Ngrid / MyBoxSize[foldnum];

  fft_forward_plan = rfftw3d_mpi_create_plan(MPI_COMM_WORLD, grid->Ngrid, grid->Ngrid, grid->Ngrid,
					     FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE | FFTW_IN_PLACE | FFTW_THREADSAFE);

  fft_inverse_plan = rfftw3d_mpi_create_plan(MPI_COMM_WORLD, grid->Ngrid, grid->Ngrid, grid->Ngrid,
					     FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE | FFTW_IN_PLACE | FFTW_THREADSAFE);

  //rfftwnd_mpi_local_sizes(fft_forward_plan, &nslab_x, &slabstart_x, &nslab_y, &slabstart_y, &fftsize);

  rfftwnd_mpi_local_sizes(fft_forward_plan, &local_nx, &local_x_start, &local_ny, &local_y_start, &total_local_size);

  rfftwnd_mpi_destroy_plan(fft_forward_plan);
  rfftwnd_mpi_destroy_plan(fft_inverse_plan);

  grid->local_x_start = local_x_start;
  grid->local_nx = local_nx;
  grid->total_local_size = grid->Ngrid * grid->Ngrid * grid->local_nx;
  grid->inv_to_slab_fac = MyBoxSize[foldnum] / grid->Ngrid;
  grid->local_ny = local_ny;
  grid->local_y_start = local_y_start;

  grid->slab_to_task = mymalloc("slab_to_task", grid->Ngrid * sizeof(int));
  grid->slabs_per_task = mymalloc("slabs_per_task", NTasks * sizeof(int));
  grid->first_slab_of_task = mymalloc("first_slab_of_task", NTasks * sizeof(int));
  slab_to_task_local = mymalloc("slab_to_task_local", grid->Ngrid * sizeof(int));

  for(ii = 0; ii < grid->Ngrid; ii++)
    slab_to_task_local[ii] = 0;
  for(ii = 0; ii < local_nx; ii++)
    slab_to_task_local[grid->local_x_start + ii] = ThisTask;

  MPI_Allreduce(slab_to_task_local, grid->slab_to_task, grid->Ngrid, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  MPI_Allgather(&local_nx, 1, MPI_INT, grid->slabs_per_task, 1, MPI_INT, MPI_COMM_WORLD);

  MPI_Allgather(&local_x_start, 1, MPI_INT, grid->first_slab_of_task, 1, MPI_INT, MPI_COMM_WORLD);

  //MPI_Allreduce(&total_local_size, &max_size, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  myfree(slab_to_task_local);
}


void density_load(char fname[], Grid *grid)
{
  FILE *fd;
  int n;
  char buf[255];
  int slabs_per_task, first_slab_of_task, ngrid;
  double box;

  sprintf(buf, "%s_density_field_%d.%d", fname, (int) grid->Ngrid, ThisTask);
  MSG0(" Starting to read the density field", t0, ThisTask);
  int ig, readgroups = floor(NTasks/NumFilesWrittenInParallel);
  for(ig = 0; ig < readgroups; ig++)
  {
    if(ThisTask % readgroups == ig)
    { 
      fd = my_fopen(buf, "r");
      my_fread(&ngrid, sizeof(int), 1, fd);
      my_fread(&n, sizeof(int), 1, fd);

      my_fread(&slabs_per_task, sizeof(int), 1, fd);
      my_fread(&first_slab_of_task, sizeof(int), 1, fd);
      my_fread(&box, sizeof(double), 1, fd);
 
      assert(grid->Ngrid == ngrid);
      assert(sizeof(MyFloat) == n);
      assert(grid->first_slab_of_task[ThisTask] == first_slab_of_task);
      assert(grid->slabs_per_task[ThisTask] == slabs_per_task);

      if(grid->local_nx > 0)
        my_fread(grid->data, sizeof(MyFloat), grid->Ngrid * grid->Ngrid * grid->local_nx, fd);

      fclose(fd);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  MSG0(" Finished reading the field", t0, ThisTask);
}


void density_save(char fname[], Grid grid)
{
  FILE *fd;
  int n;
  char buf[255];

  sprintf(buf, "%s_density_field_%d.%d", fname, (int) grid.Ngrid, ThisTask);
  MSG0(" Starting to dump the density field", t0, ThisTask);

  fd = my_fopen(buf, "w");
  my_fwrite(&grid.Ngrid, sizeof(int), 1, fd);

  n = sizeof(MyFloat);
  my_fwrite(&n, sizeof(int), 1, fd);

  my_fwrite(&grid.slabs_per_task[ThisTask], sizeof(int), 1, fd);
  my_fwrite(&grid.first_slab_of_task[ThisTask], sizeof(int), 1, fd);

  my_fwrite(&MyBoxSize[0], sizeof(double), 1, fd);

  if(grid.local_nx > 0)
    my_fwrite(grid.data, sizeof(MyFloat), grid.Ngrid * grid.Ngrid * grid.local_nx, fd);

  fclose(fd);

  MSG0(" Finished writing the field", t0, ThisTask);
}


void density(Particle * P,long long NumPart, Grid * grid, int nfold, int value)
{
  int rep;
  long long nstart, npart, ip;
  //long i, Nblock = NumPart;
  //long FirstPart = 0;

  if ( value == 0 )   MSG0(" Computing Density Field", t0, ThisTask);
  if ( value == 1 )   MSG0(" Computing Vx Field", t0, ThisTask);
  if ( value == 2 )   MSG0(" Computing Vy Field", t0, ThisTask);
  if ( value == 3 )   MSG0(" Computing Vz Field", t0, ThisTask);

  //#define CYCLES 8
  //int BufferSize = 200;

  for(rep = 0; rep < CYCLES; rep++)
    {
      if(rep == (CYCLES - 1))
        npart = NumPart - (CYCLES - 1) * (NumPart / CYCLES);
      else
        npart = (NumPart / CYCLES);

      nstart = rep * (NumPart / CYCLES);

      //if(ThisTask == 0)
      //  printf("CYCLE %d-%d nstart=%ld npart=%ld, %ld\n", rep, CYCLES, nstart, npart, NumPart);

      mapping(P, nstart, npart, grid, nfold, 0, value);
#ifdef CIC
      mapping(P, nstart, npart, grid, nfold, 1, value);
#endif
    }

//#endif
}

double density_statistics(Grid * grid, const char *label, int foldnum)
{
  long i;
  double delta;
  double local_total = 0.0, total = 0.0;
  double local_var = 0.0, var = 0.0;
  double local_m3 = 0.0, m3 = 0.0;
  double local_min = 1e20, min = 1e20;
  double local_max = -1e20, max = -1e20;
  double local_mean = 0.0, mean = 0.0;
  double local_dm3 = 0.0, dm3 = 0.0;
  double local_dm2 = 0.0, dm2 = 0.0;

  for(i = 0; i < SQR(grid->Ngrid) * (grid->local_nx); i++)
    local_mean += (double) grid->data[i];

  MPI_Allreduce(&local_mean, &mean, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  mean = mean / CUB(grid->Ngrid);

  for(i = 0; i < SQR(grid->Ngrid) * (grid->local_nx); i++)
    {
      //if (grid->data[i] > 15*mean) grid->data[i] = 15*mean;
      if(local_min > grid->data[i])
        local_min = grid->data[i];
      if(local_max < grid->data[i])
        local_max = grid->data[i];

      local_total += (double) grid->data[i];
      local_var += SQR(grid->data[i]);
      local_m3 += CUB(grid->data[i]);
    }

  MPI_Allreduce(&local_m3, &m3, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&local_var, &var, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&local_total, &total, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&local_min, &min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&local_max, &max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  mean = total / CUB(grid->Ngrid);
  grid->m1 = mean;
  grid->m2 = var / CUB(grid->Ngrid);
  grid->m3 = m3 / CUB(grid->Ngrid);
  grid->total = total;

  for(i = 0; i < SQR(grid->Ngrid) * (grid->local_nx); i++)
    {
      delta = grid->data[i] / grid->m1 - 1;
      local_dm2 += SQR(delta);
      local_dm3 += CUB(delta);
    }

  MPI_Allreduce(&local_dm2, &dm2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&local_dm3, &dm3, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  grid->dm2 = dm2 / CUB(grid->Ngrid);
  grid->dm3 = dm3 / CUB(grid->Ngrid);

  if(ThisTask == 0 && foldnum == 0)
    {
      printf(" ------------------------------------------------\n");
      printf(" Grid statistics for %s\n", label);
      printf(" ------------------------------------------------\n");
      printf(" Size: %d\n", (int) grid->Ngrid);
      printf(" Local Size: %d starting at %d \n", grid->local_nx, grid->local_x_start);
      printf(" Minimum Value: %f\n", min);
      printf(" Maximum Value: %f\n", max);
      printf(" Mean: %f\n", grid->m1);
      printf(" <S^2>: %f\n", grid->m2);
      printf(" <S^3>: %f\n", grid->m3);
      printf(" <(S/<S>-1)^2>: %f\n", grid->dm2);
      printf(" <(S/<S>-1)^3>: %f\n", grid->dm3);
      printf(" Total: %f\n", total);
      printf(" Variance: %f\n", sqrt(var / grid->Ngrid - mean * mean));
      printf(" ------------------------------------------------\n");

    }
  return (total);
}

void mapping(Particle * P, long FirstPart, long NumPart, Grid * grid, int nfold, int mode, int value)
{
  int i, j, k, itask, send_task, jpos;
  long ii, kk;
  long *counter, *send_size, *rec_size;
  long **send_buffer, **rec_buffer;
  MyFloat **send_mass_buffer, **rec_mass_buffer, **send_pos_buffer, **rec_pos_buffer;
  MyFloat mass, dx,dy,dz;
  int i1,j1,k1;
  char buf[255];
  MPI_Status mpi_status;

  counter = mymalloc("counter", NTasks * sizeof(long));
  send_size = mymalloc("send_size", NTasks * sizeof(long));
  rec_size = mymalloc("rec_size", NTasks * sizeof(long));

  DBG0("Counting how many particles we have to send. to_slab_fac=%g", dbg1, t0, ThisTask, grid->to_slab_fac);

  for(ii = 0; ii < NTasks; ii++)
    {
      rec_size[ii] = 0;
      send_size[ii] = 0;
      counter[ii] = 0;
    }

//      for (ii = 0; ii < grid->Ngrid; ii++)    
  //    printf("%d ii=%d slab_to_task[ii]=%d\n", ThisTask, ii, grid->slab_to_task[ii]);

  for(ii = FirstPart; ii < FirstPart + NumPart; ii++)
    {
      i = (grid->to_slab_fac * P[ii].pos[0]);
      if (mode == 0) 
        i = (i + grid->Ngrid) % grid->Ngrid;
      if (mode == 1) 
        i = (i + 1 + grid->Ngrid) % grid->Ngrid;

     if ( i < 0 || i > grid->Ngrid)
        printf("ii=%lld pos=%f grid=%d\n",ii,P[ii].pos[0], (int) grid->Ngrid);
      assert(i >= 0 && i < grid->Ngrid);
      itask = grid->slab_to_task[i];


      assert(itask >= 0 && itask < NTasks);
      send_size[itask]++;
    }

  for(ii = 0; ii < NTasks; ii++)
    if(send_size[ii] > 0)
      DBG("Task:%d sending:%d particles to task %d", t0, dbg2, ThisTask, send_size[ii], ii);

  DBG(" [Task=%d] We allocate a buffer to send to the processors", dbg2, t0, ThisTask);


  rec_buffer = (long **) mymalloc("rec_buffer", NTasks * sizeof(*rec_buffer));
  send_buffer = (long **) mymalloc("send_buffer", NTasks * sizeof(*send_buffer));

  rec_mass_buffer = (MyFloat **) mymalloc("rec_mass_buffer", NTasks * sizeof(*rec_mass_buffer));
  send_mass_buffer = (MyFloat **) mymalloc("send_mass_buffer", NTasks * sizeof(*send_mass_buffer));

  rec_pos_buffer = (MyFloat **) mymalloc("rec_pos_buffer", NTasks * sizeof(*rec_pos_buffer));
  send_pos_buffer = (MyFloat **) mymalloc("send_pos_buffer", NTasks * sizeof(*send_pos_buffer));

  for(ii = 0; ii < NTasks; ii++)
    send_buffer[ii] = (long *) mymalloc("send_buffer[ii]", MAX(send_size[ii], 1) * sizeof(long));

  for(ii = 0; ii < NTasks; ii++)
    send_mass_buffer[ii] = (MyFloat *) mymalloc("send_mass_buffer[ii]", MAX(send_size[ii], 1) * sizeof(MyFloat));

  for(ii = 0; ii < NTasks; ii++)
    send_pos_buffer[ii] = (MyFloat *) mymalloc("send_pos_buffer[ii]", 3 * MAX(send_size[ii], 1) * sizeof(MyFloat));

  DBG(" [Task=%d] Filling that buffer", dbg2, t0, ThisTask);
  for(ii = FirstPart; ii < FirstPart + NumPart; ii++)
    {
      if (value == 0 ) mass = P[ii].mass;

      if (FileFormat == 0 || FileFormat == 3)
          if (value == 0 ) mass = 1.0;
#ifdef DIVERGENCE
      if (value == 1) mass = P[ii].vel[0];
      if (value == 2) mass = P[ii].vel[1];
      if (value == 3) mass = P[ii].vel[2];
#endif
      i = (grid->to_slab_fac * P[ii].pos[0]);
      j = (grid->to_slab_fac * P[ii].pos[1]);
      k = (grid->to_slab_fac * P[ii].pos[2]);

      i = (i + grid->Ngrid) % grid->Ngrid;
      j = (j + grid->Ngrid) % grid->Ngrid;
      k = (k + grid->Ngrid) % grid->Ngrid;

      if (mode == 0) 
        i = (i + grid->Ngrid) % grid->Ngrid;
      else
        i = (i + 1 + grid->Ngrid) % grid->Ngrid;

      itask = grid->slab_to_task[i];

      assert(itask >= 0 && itask < NTasks);
      i -= grid->first_slab_of_task[itask];
      //assert(i >= 0 && i < grid->slabs_per_task[itask]);

      send_buffer[itask][counter[itask]] = k + grid->Ngrid * (j + grid->Ngrid * i);
      send_mass_buffer[itask][counter[itask]] = mass;
      for (jpos=0; jpos<3; jpos++)
        send_pos_buffer[itask][3*counter[itask] + jpos] = P[ii].pos[jpos];
      counter[itask]++;
    }

  for(ii = 0; ii < NTasks; ii++)
    assert(counter[ii] == send_size[ii]);

  // Sending and rceiving buffers
  DBG(" [Task=%d] Now we can send the buffers", dbg2, t0, ThisTask);
  for(ii = 0; ii < NTasks; ii++)
    {
      send_task = ii;
      if(send_task == ThisTask)
	{
	  for(kk = 0; kk < NTasks; kk++)
	    {
	      if(kk != ThisTask)
		{
		  MPI_Send(&send_size[kk], 1, MPI_LONG, kk, TAG_N, MPI_COMM_WORLD);
		  if(send_size[kk] > 0)
		    {
		      MPI_Send(send_buffer[kk], send_size[kk], MPI_LONG, kk, TAG_PDATA, MPI_COMM_WORLD);
		      MPI_Send(send_mass_buffer[kk],  send_size[kk], SEND_TYPE, kk, TAG_MDATA, MPI_COMM_WORLD);
		      MPI_Send(send_pos_buffer[kk], 3 * send_size[kk], SEND_TYPE, kk, TAG_POSDATA, MPI_COMM_WORLD);
		    }
		}
	    }
	  // to be consistent
	  //rec_buffer[ThisTask] = (long*) mymalloc( "rec_buffer[ii]",1 * sizeof(long) );
	}
      else
	{
	  MPI_Recv(&rec_size[send_task], 1, MPI_LONG, send_task, TAG_N, MPI_COMM_WORLD, &mpi_status);
	  rec_buffer[ii] = (long *) mymalloc("rec_buffer[ii]", MAX(rec_size[send_task], 1) * sizeof(long));
	  rec_mass_buffer[ii] = (MyFloat *) mymalloc("rec_mass_buffer[ii]", MAX(rec_size[send_task], 1) * sizeof(MyFloat));
	  rec_pos_buffer[ii] = (MyFloat *) mymalloc("rec_pos_buffer[ii]", 3 * MAX(rec_size[send_task], 1) * sizeof(MyFloat));
	  if(rec_size[ii] > 0)
	    {
	      DBG("Task:%d received %d particles from:%d", dbg2, t0, ThisTask, rec_size[ii], ii);
	      MPI_Recv(rec_buffer[ii], rec_size[ii], MPI_LONG, send_task, TAG_PDATA, MPI_COMM_WORLD, &mpi_status);
	      MPI_Recv(rec_mass_buffer[ii], rec_size[ii], SEND_TYPE, send_task, TAG_MDATA, MPI_COMM_WORLD, &mpi_status);
	      MPI_Recv(rec_pos_buffer[ii], 3 * rec_size[ii], SEND_TYPE, send_task, TAG_POSDATA, MPI_COMM_WORLD, &mpi_status);
	      for(kk = 0; kk < rec_size[ii]; kk++)
              {
#ifdef CIC
                 i = (int) (rec_pos_buffer[ii][3*kk+0] * grid->to_slab_fac);
                 j = (int) (rec_pos_buffer[ii][3*kk+1] * grid->to_slab_fac);
                 k = (int) (rec_pos_buffer[ii][3*kk+2] * grid->to_slab_fac);

                 dx = rec_pos_buffer[ii][3*kk+0]  * grid->to_slab_fac - i;
                 dy = rec_pos_buffer[ii][3*kk+1]  * grid->to_slab_fac - j;
                 dz = rec_pos_buffer[ii][3*kk+2]  * grid->to_slab_fac - k;


                 i = (i + grid->Ngrid) % grid->Ngrid;
                 j = (j + grid->Ngrid) % grid->Ngrid;
                 k = (k + grid->Ngrid) % grid->Ngrid;

                 i -= grid->first_slab_of_task[ThisTask];
                 i1 = i + 1;
                 j1 = j + 1;
                 k1 = k + 1;
                 if(i1 >= grid->Ngrid) i1 -= grid->Ngrid;
                 if(j1 >= grid->Ngrid) j1 -= grid->Ngrid;
                 if(k1 >= grid->Ngrid) k1 -= grid->Ngrid;

                mass = rec_mass_buffer[ii][kk];
                if ( mode == 0)
                {
                  grid->data[k + grid->Ngrid * (j +  grid->Ngrid * i)]    += mass*(1.0 - dx) * (1.0 - dy) * (1.0 - dz);
                  grid->data[k +  grid->Ngrid * (j1 +  grid->Ngrid * i)]  += mass*(1.0 - dx) * dy * (1.0 - dz);
                  grid->data[k1 +  grid->Ngrid * (j +  grid->Ngrid * i)]  += mass*(1.0 - dx) * (1.0 - dy) * dz;
                  grid->data[k1 +  grid->Ngrid * (j1 +  grid->Ngrid * i)] += mass*(1.0 - dx) * dy * dz;
                }

                if ( mode == 1)
                {
                  grid->data[k +  grid->Ngrid * (j +  grid->Ngrid * i1)]  += mass * dx * (1.0 - dy) * (1.0 - dz);
                  grid->data[k +  grid->Ngrid * (j1 +  grid->Ngrid * i1)] += mass * dx * dy * (1.0 - dz);
                  grid->data[k1 + grid->Ngrid * (j + grid->Ngrid * i1)]   += mass * dx * (1.0 - dy) * dz;
                  grid->data[k1 + grid->Ngrid * (j1 + grid->Ngrid * i1)]  += mass * dx * dy * dz;
                }
#else
		grid->data[rec_buffer[ii][kk]] += (MyFloat) rec_mass_buffer[ii][kk];
#endif
#ifdef CENTER_OF_MASS
                for (jpos=0; jpos<3; jpos++)
	         grid->pos[3*rec_buffer[ii][kk]+jpos] += ((MyFloat) rec_mass_buffer[ii][kk])*rec_pos_buffer[ii][3*kk+jpos];
#endif
              }
	    }
	  myfree(rec_pos_buffer[ii]);
	  myfree(rec_mass_buffer[ii]);
	  myfree(rec_buffer[ii]);
	}
    }

  for(kk = 0; kk < send_size[ThisTask]; kk++)
  {
      ii = ThisTask;
#ifdef CIC
      i = (int) (send_pos_buffer[ii][3*kk+0] * grid->to_slab_fac);
      j = (int) (send_pos_buffer[ii][3*kk+1] * grid->to_slab_fac);
      k = (int) (send_pos_buffer[ii][3*kk+2] * grid->to_slab_fac);

      dx = send_pos_buffer[ii][3*kk+0]  * grid->to_slab_fac - i;
      dy = send_pos_buffer[ii][3*kk+1]  * grid->to_slab_fac - j;
      dz = send_pos_buffer[ii][3*kk+2]  * grid->to_slab_fac - k;

      i = (i + grid->Ngrid) % grid->Ngrid;
      j = (j + grid->Ngrid) % grid->Ngrid;
      k = (k + grid->Ngrid) % grid->Ngrid;

      i -= grid->first_slab_of_task[ThisTask];
      i1 = i + 1;
      if(i1 >= grid->Ngrid) i1 -= grid->Ngrid;
      j1 = j + 1;
      if(j1 >= grid->Ngrid) j1 -= grid->Ngrid;
      k1 = k + 1;
      if(k1 >= grid->Ngrid) k1 -= grid->Ngrid;
                 
      mass = send_mass_buffer[ii][kk];
      if ( mode == 0)
      {
        grid->data[k + grid->Ngrid * (j +  grid->Ngrid * i)]    += mass*(1.0 - dx) * (1.0 - dy) * (1.0 - dz);
        grid->data[k +  grid->Ngrid * (j1 +  grid->Ngrid * i)]  += mass*(1.0 - dx) * dy * (1.0 - dz);
        grid->data[k1 +  grid->Ngrid * (j +  grid->Ngrid * i)]  += mass*(1.0 - dx) * (1.0 - dy) * dz;
        grid->data[k1 +  grid->Ngrid * (j1 +  grid->Ngrid * i)] += mass*(1.0 - dx) * dy * dz;
      }
      if ( mode == 1)
      {
         grid->data[k +  grid->Ngrid * (j +  grid->Ngrid * i1)]  += mass * dx * (1.0 - dy) * (1.0 - dz);
         grid->data[k +  grid->Ngrid * (j1 +  grid->Ngrid * i1)] += mass * dx * dy * (1.0 - dz);
         grid->data[k1 + grid->Ngrid * (j + grid->Ngrid * i1)]   += mass * dx * (1.0 - dy) * dz;
         grid->data[k1 + grid->Ngrid * (j1 + grid->Ngrid * i1)]  += mass * dx * dy * dz;
      }
#else
    grid->data[send_buffer[ii][kk]] += (MyFloat) send_mass_buffer[ii][kk];
#endif
#ifdef CENTER_OF_MASS
    for (jpos=0; jpos<3; jpos++)
        grid->pos[3*send_buffer[ii][kk]+jpos] += ((MyFloat) send_mass_buffer[ii][kk]) * send_pos_buffer[ii][3*kk+jpos];
#endif
  }

  //report_memory_usage(&HighMark, "DENSITY");

  DBG0(" Freeing  buffers", dbg1, t0, ThisTask);

  for(ii = NTasks - 1; ii >= 0; ii--)
    myfree(send_pos_buffer[ii]);

  for(ii = NTasks - 1; ii >= 0; ii--)
    myfree(send_mass_buffer[ii]);

  for(ii = NTasks - 1; ii >= 0; ii--)
    myfree(send_buffer[ii]);
  
  myfree(send_pos_buffer);
  myfree(rec_pos_buffer);

  myfree(send_mass_buffer);
  myfree(rec_mass_buffer);

  myfree(send_buffer);
  myfree(rec_buffer);

  myfree(rec_size);
  myfree(send_size);
  myfree(counter);
}

void grid_clean(Grid grid)
{
  DBG0(" Cleaning plan for cells", dbg1, t0, ThisTask);
  if(grid.local_nx > 0)
  {
#ifdef CENTER_OF_MASS
    myfree(grid.pos);
#endif
    myfree(grid.data);
  }
  myfree(grid.first_slab_of_task);
  myfree(grid.slabs_per_task);
  myfree(grid.slab_to_task);

}

/*
Particle *create_density(Grid N_grid, float *bias, float N2, Grid DM_grid, long *N_NPart)
{
  int i, j, k, iseed;
  long n, Total = 0;
  int xx, yy, zz;
  double delta_h, delta_dm;
  double sigma, S;
  gsl_rng *random_generator;
  Particle *NH;

  iseed = 1223443;
  random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);
  gsl_rng_set(random_generator, iseed);

  MSG0(" Creating Density Field", t0, ThisTask);
  //printf("%d %d %ld\n",N_grid.Ngrid,N_grid.local_nx,SQR(N_grid.Ngrid)*N_grid.local_nx);
  sigma = DM_grid.m2 / SQR(DM_grid.m1) - 1;
  S = DM_grid.m3 / CUB(DM_grid.m1) - 3 * DM_grid.m2 / SQR(DM_grid.m1) + 2;
  sigma = DM_grid.dm2;
  S = DM_grid.dm3;

  for(i = 0; i < N_grid.Ngrid * N_grid.Ngrid * N_grid.local_nx; i++)
    {
      delta_dm = DM_grid.data[i] / DM_grid.m1 - 1;
      //if (delta_dm > 15) delta_dm = 15;
      //delta_h = bias[0] + bias[1]*delta_dm + bias[2]/2.*SQR(delta_dm) + bias[3]/6.*CUB(delta_dm);
      delta_h = bias[1] * delta_dm + bias[2] / 2. * (SQR(delta_dm) - sigma) + bias[3] / 6. * (CUB(delta_dm) - S);
      //delta_h = bias[1]*delta_dm+ bias[2]/2.*(SQR(delta_dm)-sigma);// + bias[3]/6.*(CUB(delta_dm)-S);
      N_grid.data[i] = (delta_h + 1.0) * N2;
      //N_grid.data[i] = delta_h;
      //N_grid.data[i] = ( delta_dm + 1.0 )*N2;

      if(N_grid.data[i] < 0.0)
	N_grid.data[i] = 0.0;
      else
	N_grid.data[i] = gsl_ran_poisson(random_generator, N_grid.data[i]);

    }

  for(i = 0; i < N_grid.Ngrid * N_grid.Ngrid * N_grid.local_nx; i++)
    Total += N_grid.data[i];

  *N_NPart = Total;
  NH = (Particle *) my_malloc(Total * sizeof(Particle));

  long ii = 0;
  for(i = 0; i < N_grid.local_nx; i++)
    for(j = 0; j < N_grid.Ngrid; j++)
      for(k = 0; k < N_grid.Ngrid; k++)
	{
	  for(n = 0; n < (int) N_grid.data[k + N_grid.Ngrid * (j + i * N_grid.Ngrid)]; n++)
	    {
	      NH[ii].pos[0] = (i + N_grid.local_x_start + gsl_rng_uniform(random_generator)) * N_grid.inv_to_slab_fac;
	      NH[ii].pos[1] = (j + gsl_rng_uniform(random_generator)) * N_grid.inv_to_slab_fac;
	      NH[ii].pos[2] = (k + gsl_rng_uniform(random_generator)) * N_grid.inv_to_slab_fac;

	      xx = (N_grid.to_slab_fac * NH[ii].pos[0]);
	      yy = (N_grid.to_slab_fac * NH[ii].pos[1]);
	      zz = (N_grid.to_slab_fac * NH[ii].pos[2]);

	      xx = (xx + N_grid.Ngrid) % N_grid.Ngrid;	// to guarantee the result is positive
	      yy = (yy + N_grid.Ngrid) % N_grid.Ngrid;
	      zz = (zz + N_grid.Ngrid) % N_grid.Ngrid;
	      xx -= N_grid.local_x_start;
	      //NH[counter].mass = exp(log(m1) + gsl_rng_uniform(random_generator)*(log(m2)-log(m1)));
	      if((xx != i || yy != j || zz != k) && dbg1 == 1)
		{
		  printf("%d %d - %d %f\n", xx, i, N_grid.local_x_start, NH[ii].pos[0]);
		  printf("%d %d\n", yy, j);
		  printf("%d %d\n", zz, k);
		}
	      ii++;
	    }
	}


  DBG("totales %ld %ld", dbg1, t0, ii, Total);

  gsl_rng_free(random_generator);

  return (NH);
}
*/
