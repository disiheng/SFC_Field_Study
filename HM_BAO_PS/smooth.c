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

#include "../utils/macros.h"
#include "../utils/utils.h"
#include "../utils/gadget.h"
#include "../utils/mymalloc.h"
#include "density.h"

#define TAG_N          10
#define TAG_PDATA      11

//extern int max_size, NProcRead;
//extern double DiluteFactor, invCellSize;
//extern double rms, mean, sd, skew, kurtosis;
//extern time_t t0;

//int local_x_start, local_nx, local_y_start, local_ny, total_local_size;
//int *slab_to_task, *slab_to_task_local, *slabs_per_task, *first_slab_of_task;
typedef long long large_array_offset;

extern time_t t0;
extern double BoxSize;
extern int ThisTask, NTasks;

#define dbg1 0
#define dbg2 0

static rfftwnd_mpi_plan fft_forward_plan, fft_inverse_plan;

void smooth(Grid * orig, Grid * smooth, float smoothing)
{
  int ndim;
  int *slab_to_task_local, *slab_to_task, *slabs_per_task, *first_slab_of_task;
  fftw_complex *fft_tmpgrid;
  fftw_real *tmpgrid;
  fftw_real *workspace;

  double k2, kx, ky, kz, smth, kvec, hubble;
  double dx, dy, dz;
  double fx, fy, fz, ff;
  double asmth2, fac, fac2;
  int i, j, k, ndim2, sendTask, recvTask;
  int x, y, z, ip, dim;
  int task0, task1, nimport, nexport;
  int nslab_x, nslab_y, slabstart_x, slabstart_y;
  int slab_x, slab_y, slab_z;
  int slab_xx, slab_yy, slab_zz;
  int smallest_slab;
  int fftsize, maxfftsize;
  int ngrp;
  int *send_count_bak;
  char buf[500];
  float to_slab_fac;
  double Norm_local = 0.0, Normalization = 0.0, win;
  double total;


  /* Workspace out the ranges on each processor. */
  ndim = orig->Ngrid;
  fft_forward_plan = rfftw3d_mpi_create_plan(MPI_COMM_WORLD, ndim, ndim, ndim, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE | FFTW_IN_PLACE | FFTW_THREADSAFE);

  fft_inverse_plan = rfftw3d_mpi_create_plan(MPI_COMM_WORLD, ndim, ndim, ndim, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE | FFTW_IN_PLACE | FFTW_THREADSAFE);

  rfftwnd_mpi_local_sizes(fft_forward_plan, &nslab_x, &slabstart_x, &nslab_y, &slabstart_y, &fftsize);


  slab_to_task = (int *) mymalloc("slab_to_task", ndim * sizeof(int));
  slab_to_task_local = (int *) mymalloc("slab_to_task_local", ndim * sizeof(int));
  for(i = 0; i < ndim; i++)
    slab_to_task_local[i] = 0;

  for(i = 0; i < nslab_x; i++)
    slab_to_task_local[slabstart_x + i] = ThisTask;

  MPI_Allreduce(slab_to_task_local, slab_to_task, ndim, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  myfree(slab_to_task_local);

  MPI_Allreduce(&nslab_x, &smallest_slab, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

  slabs_per_task = (int *) mymalloc("slabs_per_task", NTasks * sizeof(int));
  MPI_Allgather(&nslab_x, 1, MPI_INT, slabs_per_task, 1, MPI_INT, MPI_COMM_WORLD);

  first_slab_of_task = (int *) mymalloc("first_slab_of_task", NTasks * sizeof(int));
  MPI_Allgather(&slabstart_x, 1, MPI_INT, first_slab_of_task, 1, MPI_INT, MPI_COMM_WORLD);

  to_slab_fac = (float) (ndim / BoxSize);
  MPI_Allreduce(&fftsize, &maxfftsize, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  fac = density_statistics(orig, "Smooth grid");
  //fac = CUB((double) ndim);


  //just to be sure
  assert(nslab_x == orig->local_nx);
  assert(slabstart_x == orig->local_x_start);
  assert(nslab_y == orig->local_ny);
  assert(slabstart_y == orig->local_y_start);
  assert(slabs_per_task[ThisTask] == orig->slabs_per_task[ThisTask]);

  if(ThisTask == 0)
    {
      printf("Starting Smoothing calculation.  (presently allocated=%g MB)\n", AllocatedBytes / (1024.0 * 1024.0));
      printf("Normalization: %g\n", fac);
      fflush(stdout);
    }

  if(ThisTask == 0)
    printf("fft forward (%d) \n", maxfftsize);

  workspace = (fftw_real *) mymalloc("workspace", maxfftsize * sizeof(fftw_real));

  tmpgrid = (fftw_real *) mymalloc("tmpgrid", maxfftsize * sizeof(fftw_real));
  memset(tmpgrid, 0, maxfftsize * sizeof(fftw_real));

  report_memory_usage(&HighMark, "SMOOTH");


  ndim2 = (2 * (ndim / 2 + 1));

  for(i = 0; i < slabs_per_task[ThisTask]; i++)
    for(j = 0; j < ndim; j++)
      for(k = 0; k < ndim; k++)
	tmpgrid[((large_array_offset) ndim2) * (ndim * i + j) + k] = orig->data[ndim * (ndim * i + j) + k];

  asmth2 = (2 * M_PI) * smoothing / BoxSize;
  asmth2 *= asmth2;

  if(ThisTask == 0)
    printf("FFT backwards\n");

  fft_tmpgrid = (fftw_complex *) & tmpgrid[0];
  rfftwnd_mpi(fft_forward_plan, 1, tmpgrid, workspace, FFTW_TRANSPOSED_ORDER);


  for(y = slabstart_y; y < slabstart_y + nslab_y; y++)
    for(x = 0; x < ndim; x++)
      for(z = 0; z < ndim / 2 + 1; z++)
	{
	  if(x > ndim / 2)
	    kx = x - ndim;
	  else
	    kx = x;
	  if(y > ndim / 2)
	    ky = y - ndim;
	  else
	    ky = y;
	  if(z > ndim / 2)
	    kz = z - ndim;
	  else
	    kz = z;

	  k2 = kx * kx + ky * ky + kz * kz;

	  if(k2 > 0)
	    {

	      ip = ndim * (ndim / 2 + 1) * (y - slabstart_y) + (ndim / 2 + 1) * x + z;

	      win = exp(-k2 * asmth2 / 2);
	      Norm_local += win;

	      fft_tmpgrid[ip].re *= (fftw_real) win;
	      fft_tmpgrid[ip].im *= (fftw_real) win;

	    }

	}

  //Norm_local = 1./NTask;


  MPI_Allreduce(&Norm_local, &Normalization, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  if(ThisTask == 0)
    printf("Normalization = %g\n", Normalization);

  for(y = slabstart_y; y < slabstart_y + nslab_y; y++)
    for(x = 0; x < ndim; x++)
      for(z = 0; z < ndim / 2 + 1; z++)
	{

	  ip = ndim * (ndim / 2 + 1) * (y - slabstart_y) + (ndim / 2 + 1) * x + z;

	  //fft_tmpgrid[ip].re *= 1./fac;///TotNumParticles;//1/fac;
	  //fft_tmpgrid[ip].im *= 1./fac;///TotNumParticles;//1/fac;

	}

  if(ThisTask == 0)
    printf("FFT backwards\n");
  rfftwnd_mpi(fft_inverse_plan, 1, tmpgrid, workspace, FFTW_TRANSPOSED_ORDER);

  for(i = 0; i < slabs_per_task[ThisTask]; i++)
    for(j = 0; j < ndim; j++)
      for(k = 0; k < ndim; k++)
	smooth->data[ndim * (ndim * i + j) + k] = tmpgrid[((large_array_offset) ndim2) * (ndim * i + j) + k] / fac;

  rfftwnd_mpi_destroy_plan(fft_forward_plan);
  rfftwnd_mpi_destroy_plan(fft_inverse_plan);

  myfree(tmpgrid);
  myfree(workspace);
  myfree(first_slab_of_task);
  myfree(slabs_per_task);
  myfree(slab_to_task);

}
