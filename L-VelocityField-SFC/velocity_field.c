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

/* FFTW3 MPI: fftw3-mpi.h is included via allvars.h */


void init_grid(void)
{
  int i;
  int *slab_to_task_local;
  ptrdiff_t alloc_local;

  /* FFTW3 MPI: create plans for 3D real-to-complex and complex-to-real.
   * Use transposed output for efficient MPI data distribution.
   * In-place: pass same pointer for in and out (FFTW3 handles it). */
  fft_forward_plan = fftw_mpi_plan_dft_r2c_3d(
      Ndim, Ndim, Ndim,
      (double *)Grid, (fftw_complex *)Grid,
      MPI_COMM_WORLD, FFTW_ESTIMATE);

  fft_inverse_plan = fftw_mpi_plan_dft_c2r_3d(
      Ndim, Ndim, Ndim,
      (fftw_complex *)Grid, (double *)Grid,
      MPI_COMM_WORLD, FFTW_ESTIMATE);

  /* Get local data distribution */
  alloc_local = fftw_mpi_local_size_3d_transposed(
      Ndim, Ndim, Ndim, MPI_COMM_WORLD,
      &local_n0, &local_0_start, &local_n1, &local_1_start);

  nslab_x = (int)local_n0;
  slabstart_x = (int)local_0_start;
  nslab_y = (int)local_n1;
  slabstart_y = (int)local_1_start;
  fftsize = alloc_local;

  if(ThisTask == 0)
    printf("fftsize=%td nslab_x=%d slabstart_x=%d nslab_y=%d slabstart_y=%d\n",
           fftsize, nslab_x, slabstart_x, nslab_y, slabstart_y);

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
  {
    ptrdiff_t local_max = fftsize;
    MPI_Allreduce(&local_max, &maxfftsize, 1, MPI_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
  }

  /* Allocate Grid buffer: must hold both real input (local_n0*Ndim*Ndim doubles)
   * and complex output (fftsize * 2 doubles, since fftsz=alloc_local in fftw_complex units).
   * In FFTW3 MPI, no external workspace is needed (handled internally). */
  {
    size_t real_bytes  = (size_t)local_n0 * (size_t)Ndim * (size_t)Ndim * sizeof(double);
    size_t cplx_bytes  = (size_t)alloc_local * 2 * sizeof(double);
    size_t grid_bytes  = (real_bytes > cplx_bytes) ? real_bytes : cplx_bytes;

#ifdef GEN_HALOFIELD
    Grid = (fftw_real *) mymalloc("Grid", grid_bytes);
#else
    Grid = (fftw_real *) mymalloc("Grid", grid_bytes);
    workspace = (fftw_real *) mymalloc("workspace", grid_bytes);
#endif
  }
}

void clean_grid()
{
  fftw_destroy_plan(fft_forward_plan);
  fftw_destroy_plan(fft_inverse_plan);

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

  Vel = (fftw_real *) mymalloc("Vel", (size_t)fftsize * sizeof(fftw_complex));

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
			fft_Vel[0][0] = fft_Vel[0][1] = 0.0;
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
				fft_Vel[ip][0] = - div[1] * kx / k2 * fac2 * fac * win;
				fft_Vel[ip][1] =   div[0] * kx / k2 * fac2 * fac * win;
			  }
			  if (dim == 1){
				fft_Vel[ip][0] = - div[1] * ky / k2 * fac2 * fac * win;
				fft_Vel[ip][1] =   div[0] * ky / k2 * fac2 * fac * win;
			  }
			  if (dim == 2){
				fft_Vel[ip][0] = - div[1] * kz / k2 * fac2 * fac * win;
				fft_Vel[ip][1] =   div[0] * kz / k2 * fac2 * fac * win;
			  }
		  }
   }
 
  MPI_Allreduce(&Norm_local, &Normalization, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  if(ThisTask == 0)
	printf("Normalization = %g\n",Normalization*fac);

  if(ThisTask == 0)
      printf("FFT backwards\n");
  fftw_execute(fft_inverse_plan);

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

