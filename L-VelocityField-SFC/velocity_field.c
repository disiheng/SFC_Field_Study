#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <signal.h>

#include "allvars.h"
#include "proto.h"

#define WORKSIZE 100000

#ifdef NOTYPEPREFIX_FFTW
#include      <rfftw_mpi.h>
#else
#ifdef DOUBLEPRECISION_FFTW
#include     <drfftw_mpi.h> /* double precision FFTW */
#else
#include     <srfftw_mpi.h>
#endif
#endif

//#define  Ndim2 (2*(Ndim/2 + 1))


void init_grid(void)
{
  int i;
  int *slab_to_task_local;

 
  /* Workspace out the ranges on each processor. */
  fft_forward_plan = rfftw3d_mpi_create_plan(MPI_COMM_WORLD, Ndim, Ndim, Ndim,
                         FFTW_REAL_TO_COMPLEX,
                         FFTW_ESTIMATE | FFTW_IN_PLACE | FFTW_THREADSAFE);

  fft_inverse_plan = rfftw3d_mpi_create_plan(MPI_COMM_WORLD, Ndim, Ndim, Ndim,
                         FFTW_COMPLEX_TO_REAL,
                         FFTW_ESTIMATE | FFTW_IN_PLACE | FFTW_THREADSAFE);

  rfftwnd_mpi_local_sizes(fft_forward_plan, &nslab_x, &slabstart_x, &nslab_y, &slabstart_y, &fftsize);

  if(ThisTask == 0)
    printf("fftsize=%d nslab_y=%d slabstart_y=%d, nslab_y=%d slabstart_y=%d \n", fftsize, nslab_x,
       slabstart_x, nslab_y, slabstart_y);

  slab_to_task = (int *) mymalloc("slab_to_task", Ndim * sizeof(int));
  slab_to_task_local = (int *) mymalloc("slab_to_task_local", Ndim * sizeof(int));
  for(i = 0; i < Ndim; i++)
    slab_to_task_local[i] = 0;

  for(i = 0; i < nslab_x; i++)
    slab_to_task_local[slabstart_x + i] = ThisTask;

  MPI_Allreduce(slab_to_task_local, slab_to_task, Ndim, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  myfree(slab_to_task_local);


  MPI_Allreduce(&nslab_x, &smallest_slab, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

  slabs_per_task = (int *) mymalloc("slabs_per_task", NTask * sizeof(int));
  MPI_Allgather(&nslab_x, 1, MPI_INT, slabs_per_task, 1, MPI_INT, MPI_COMM_WORLD);

  first_slab_of_task = (int *) mymalloc("first_slab_of_task", NTask * sizeof(int));
  MPI_Allgather(&slabstart_x, 1, MPI_INT, first_slab_of_task, 1, MPI_INT, MPI_COMM_WORLD);

  to_slab_fac = (float) (Ndim / BoxSize);
  MPI_Allreduce(&fftsize, &maxfftsize, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

#ifdef GEN_HALOFIELD
  Grid = (fftw_real *) mymalloc("Grid", Ndim * Ndim * nslab_x * sizeof(fftw_real));
#else
  workspace = (fftw_real *) mymalloc("workspace", maxfftsize * sizeof(fftw_real));
  Grid = (fftw_real *) mymalloc("Grid", maxfftsize * sizeof(fftw_real));
#endif
}

void clean_grid()
{
  rfftwnd_mpi_destroy_plan(fft_forward_plan);
  rfftwnd_mpi_destroy_plan(fft_inverse_plan);

#ifdef GEN_HALOFIELD
  myfree(Grid);
#else
  myfree(nslabs_file);
  myfree(first_slab_file);
  myfree(Grid);
  myfree(workspace);
#endif
  myfree(first_slab_of_task);
  myfree(slabs_per_task);
  myfree(slab_to_task);
}


void velocity_field(int dim)
{
  double k2, kx, ky, kz, smth, kvec, hubble;
  double fx, fy, fz, ff;
  double asmth2, fac, fac2, delk;
  int i, j, k, Ndim2, sendTask, recvTask;
  int x, y, z, xx, yy, zz, ip;
  int task0, task1, nimport, nexport;
  int slab_x, slab_y, slab_z;
  int slab_xx, slab_yy, slab_zz;
  int ngrp;
  int *send_count_bak;
  char buf[500];
  double Norm_local=0.0,Normalization=0.0, win;
  double beta, Hubble;
  fftw_complex div;

  fac = 1./CUB((double)Ndim);
  delk = 2 * M_PI / BoxSize;

  Hubble = hubble_function(Time);
  beta = F_Omega(Time);
  fac2 = Hubble * beta;


  if(ThisTask == 0)
    {
      printf("Starting Reconstructing Velocity Field.  (presently allocated=%g MB)\n", AllocatedBytes / (1024.0 * 1024.0));

      if(dim == 0) printf("computing Vx Field \n");
      if(dim == 1) printf("computing Vy Field \n");
      if(dim == 2) printf("computing Vz Field \n");

      printf("Normalization: %g\n",fac);
    }

  Vel = (fftw_real *) mymalloc("Vel", maxfftsize * sizeof(fftw_real));

  report_memory_usage(&HighMark_run, "VEL_RECON");

  asmth2 = Smoothing;
  asmth2 *= asmth2;

  fft_Vel = (fftw_complex *) &Vel[0];

  for(y = slabstart_y; y < slabstart_y + nslab_y; y++)
    for(x = 0; x < Ndim; x++)
      for(z = 0; z < Ndim / 2 + 1; z++)
        {
          if(x > Ndim / 2)
             xx = x - Ndim;
          else
             xx = x;
          if(y > Ndim / 2)
             yy = y - Ndim;
          else
             yy = y;
          if(z > Ndim / 2)
             zz = z - Ndim;
          else
             zz = z;

		  kx = (float) xx * delk;
		  ky = (float) yy * delk;
		  kz = (float) zz * delk;

          k2 = kx * kx + ky * ky + kz * kz;

		  if(k2 == 0)
			fft_Vel[0].re = fft_Vel[0].im = 0.0;
          if(k2 > 0)
          {
              //do deconvolution
              fx = fy = fz = 1;
              if(kx != 0)
              {
                  fx = (M_PI * xx) / (float)Ndim;
                  fx = sin(fx) / fx;
              }
              if(ky != 0)
              {
                  fy = (M_PI * yy) / (float)Ndim;
                  fy = sin(fy) / fy;
              }
              if(kz != 0)
              {
                  fz = (M_PI * zz) / (float)Ndim;
                  fz = sin(fz) / fz;
              }
              ff = 1. / (fx * fy * fz);

			  ip = Ndim * (Ndim / 2 + 1) * (y - slabstart_y) + (Ndim / 2 + 1) * x + z;

			  win = exp(-k2 * asmth2 / 2) * ff * ff * ff * ff; // for CIC Method
			  Norm_local += win; 

              div = fft_Grid[ip];

			  // v_x(k) = -i (H*f) * k_x / |k|^2 * delta(k)
			  if (dim == 0){
				fft_Vel[ip].re = - div.im * kx / k2 * fac2 * fac * win;
				fft_Vel[ip].im =   div.re * kx / k2 * fac2 * fac * win;
			  }
			  if (dim == 1){
				fft_Vel[ip].re = - div.im * ky / k2 * fac2 * fac * win;
				fft_Vel[ip].im =   div.re * ky / k2 * fac2 * fac * win;
			  }
			  if (dim == 2){
				fft_Vel[ip].re = - div.im * kz / k2 * fac2 * fac * win;
				fft_Vel[ip].im =   div.re * kz / k2 * fac2 * fac * win;
			  }
		  }
   }
 
  MPI_Allreduce(&Norm_local, &Normalization, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  if(ThisTask == 0)
	printf("Normalization = %g\n",Normalization*fac);

  if(ThisTask == 0)
      printf("FFT backwards\n");
  rfftwnd_mpi(fft_inverse_plan, 1, Vel, workspace, FFTW_TRANSPOSED_ORDER);

}


double D1(double astart, double aend)
{
  double Dplus;
  Dplus = growth(astart) / growth(aend);
  return Dplus;
}

double D2(double astart, double aend)
{
  double omega_a, Dplus, a;

  a = astart;
  omega_a = Omega0 / (Omega0 + a * a * a * (1 - Omega0));
  Dplus = growth(astart) / growth(aend);
  return -3./7. * Dplus * Dplus * pow(omega_a, -1./143.);
}


double growth(double a)
{
  double hubble_a;

  hubble_a = sqrt(Omega0 / (a * a * a) + (1 - Omega0));

  double result, abserr;
  gsl_integration_workspace *workspace;
  gsl_function F;

  workspace = gsl_integration_workspace_alloc(WORKSIZE);

  F.function = &growth_int;

  gsl_integration_qag(&F, 0, a, 0, 1.0e-8, WORKSIZE, GSL_INTEG_GAUSS41, workspace, &result, &abserr);

  gsl_integration_workspace_free(workspace);

  return hubble_a * result;
}


double growth_int(double a, void *param)
{
  return pow(a / (Omega0 + (1 - Omega0) * a * a * a), 1.5);
}


double F_Omega(double a)
{
  double omega_a;

  omega_a = Omega0 / (Omega0 + a * a * a * (1 - Omega0));

  return pow(omega_a, 0.55);
}


double F2_Omega(double a)
{
  double omega_a;

  omega_a = Omega0 / (Omega0 + a * a * a * (1 - Omega0));

  return 2 * pow(omega_a, 4./7.);
}


double hubble_function(double a)
{
  double hubble_a;

  hubble_a = Omega0 / (a * a * a) + (1 - Omega0);
  hubble_a = 100. * sqrt(hubble_a);

  return hubble_a;
}

