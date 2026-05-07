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

#ifdef NOTYPEPREFIX_FFTW
#include      <rfftw_mpi.h>
#else
#ifdef DOUBLEPRECISION_FFTW
#include     <drfftw_mpi.h> /* double precision FFTW */
#else
#include     <srfftw_mpi.h>
#endif
#endif

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
  MPI_Comm_size(MPI_COMM_WORLD, &NTask);
  int i;
  for(PTask = 0; NTask > (1 << PTask); PTask++);

  if(argc != 3)
  {
      if (ThisTask == 0)
      {
          printf("\n  usage: L-VelocityField-Halo-MXXL  <parameterfile>  <outputnum>\n");
          printf("  <parameterfile>    see readparmeterfile.c\n");
          printf("  <outputnum>        snapshot number\n\n");
      }
      exit(1);
    }

  read_parameter_file(argv[1]);
  SnapshotNum = atoi(argv[2]);

  HighMark_run = 0;

  mymalloc_init();

  init_grid();
 
#ifdef GEN_HALOFIELD
#ifndef DMP
  /* --- load in cata data --- */
  TotNumHalo = 0;
  load_catalogue();

  for (Dim = 0; Dim < 4; Dim ++)
  {
	memset(Grid, 0, maxfftsize * sizeof(fftw_real));
	density(PP, &Grid, LocNgroups, 1, Dim);
	save_grids(Dim);
  }

  myfree(PP);
#else
  load_snapshot();
#endif
#else
  /* --- load in mesh data --- */
  get_slabs_file();
#if defined (SMOOTH_HALOFIELD) && defined (SMOOTH_HALOVELFIELD)
  for (Dim = 0; Dim < 3; Dim ++)
  {
    if (ThisTask == 0)
    {
      if(Dim == 0) printf("Smoothing Vx Field \n");
      if(Dim == 1) printf("Smoothing Vy Field \n");
      if(Dim == 2) printf("Smoothing Vz Field \n");
    }
#endif
  memset(Grid, 0, maxfftsize * sizeof(fftw_real));
  load_grid_field(Dim);

  /* --- fft dens data --- */
  fft_Grid = (fftw_complex *) &Grid[0];
  if(ThisTask == 0)
      printf("fft forward \n");
  rfftwnd_mpi(fft_forward_plan, 1, Grid, workspace, FFTW_TRANSPOSED_ORDER);

  for (rep=0; rep<=0; rep++)
  {  
      if ( rep == 0 ){ 
          Smoothing = 0.0000001;
      } else {
          Smoothing = 1.25* pow(2, rep-1);
      }

      if (ThisTask == 0)
          printf("Smoothing :%g\n",Smoothing);
#ifdef SMOOTH_HALOFIELD
	  smooth();
	  save_grids(Dim);
	  myfree(sGrid);
#endif
#ifdef LINEAR_VELFIELD
	  for (Dim=0; Dim<3; Dim++)
	  {
        velocity_field(Dim);
		save_grids(Dim);
		myfree(Vel);
	  }
#endif

#if defined (SMOOTH_HALOFIELD) && defined (SMOOTH_HALOVELFIELD)
  }
#endif
  }
#endif

  clean_grid();

 
  MPI_Barrier(MPI_COMM_WORLD);
  if(ThisTask == 0)
      printf("Done\n"); 

  MPI_Finalize();		/* clean up & finalize MPI */

  return 0;
}


double density_statistics(fftw_real *grid)
{
    long long i;
    double delta;
    double local_total = 0.0, total = 0.0;
    double local_var = 0.0, var = 0.0;
    double local_m3 = 0.0, m3 = 0.0, m1=0.0, m2=0.0;
    double local_min = 1e20, min = 1e20;
    double local_max = -1e20, max = -1e20;
    double local_mean = 0.0, mean = 0.0;
    double local_dm3 = 0.0, dm3 = 0.0;
    double local_dm2 = 0.0, dm2 = 0.0;

    for (i=0; i < SQR(Ndim) * (nslab_x); i++)
        local_mean += (double) grid[i];
    MPI_Allreduce (&local_mean, &mean, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    mean = mean/CUB(Ndim);

    for (i=0; i < SQR(Ndim) * (nslab_x); i++)
    {
        //if (grid->data[i] > 15*mean) grid->data[i] = 15*mean;
        if (local_min > grid[i])
            local_min = grid[i];
        if (local_max < grid[i])
            local_max = grid[i];

        local_total += (double) grid[i];
        local_var += SQR(grid[i]);
        local_m3 += CUB(grid[i]);
    }

    MPI_Allreduce (&local_m3, &m3, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce (&local_var, &var, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce (&local_total, &total, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce (&local_min, &min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce (&local_max, &max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    mean = total/CUB(Ndim);
    m1 = mean;
    m2 = var/CUB(Ndim);
    m3 = m3/CUB(Ndim);
    total = total;

    for (i=0; i < SQR(Ndim) * (nslab_x); i++)
    {
//        delta = grid[i]/m1 - 1;
        delta = grid[i];
        local_dm2 += SQR(delta);
        local_dm3 += CUB(delta);
    }

    MPI_Allreduce (&local_dm2, &dm2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce (&local_dm3, &dm3, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    dm2 = dm2/CUB(Ndim);
    dm3 = dm3/CUB(Ndim);

    if (ThisTask == 0)
    {
        printf(" ------------------------------------------------\n");
        printf(" Grid statistics\n");
        printf(" Size: %ld\n", Ndim);
        printf(" Local Size: %d starting at %d \n", nslab_x, slabstart_x);
        printf(" Minimun Value: %f\n", min);
        printf(" Maximun Value: %f\n", max);
        printf(" Mean: %f\n", m1);
        printf(" <S^2>: %f\n", m2);
        printf(" <S^3>: %f\n", m3);
        printf(" <(S/<S>-1)^2>: %f\n", dm2);
        printf(" <(S/<S>-1)^3>: %f\n", dm3);
        printf(" Total: %f\n", total);
        printf(" Variance: %f\n", sqrt(var/Ndim - mean*mean));
        printf(" ------------------------------------------------\n");

    }

  return total;

}

double log_density(float *grid)
{
    long long i;
    double local_total = 0.0, total = 0.0;
    double local_log_total = 0.0, log_total = 0.0;
    double log_mean, mean, true_mean;

    true_mean = (double) TotNumHalo / CUB(Ndim);

    for (i=0; i < SQR(Ndim) * (nslab_x); i++)
    {
        local_log_total += (double) log(grid[i]/true_mean);
        local_total     += (double) (grid[i]);
    }

    MPI_Allreduce (&local_total, &total, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce (&local_log_total, &log_total, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    mean = total/CUB(Ndim);
    log_mean = log_total/CUB(Ndim);

    if (ThisTask == 0)
    {
        printf(" ------------------------------------------------\n");
        printf(" Grid statistics\n");
        printf(" Size: %ld\n", Ndim);
        printf(" Local Size: %d starting at %d \n", nslab_x, slabstart_x);
        printf(" Mean: %f\n", mean);
        printf(" Log-Mean: %f\n", log_mean);
        printf(" True Mean: %f\n", true_mean);
        printf(" ------------------------------------------------\n");
    }

  return log_mean;

}


