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

void smooth(void)
{
  double k2, kx, ky, kz, smth, delk;
  double fx, fy, fz, ff;
  double asmth2, fac, fac2;
  int i, j, k, Ndim2, sendTask, recvTask;
  int x, y, z, xx, yy, zz, ip;
  int task0, task1, nimport, nexport;
  int slab_x, slab_y, slab_z;
  int slab_xx, slab_yy, slab_zz;
  int ngrp;
  int *send_count_bak;
  char buf[500];
  double Norm_local=0.0,Normalization=0.0, win;

  delk = 2 * M_PI / BoxSize;
  fac2 = 1./CUB((double)Ndim);
  fac  = TotNumHalo / CUB((double)Ndim);

  if(ThisTask == 0)
    {
      printf("Starting Smoothing calculation. (presently allocated=%g MB)\n", AllocatedBytes / (1024.0 * 1024.0));
      printf("Normalization: %g\n",fac);
      fflush(stdout);
    }

  sGrid = (fftw_real *) mymalloc("sGrid", maxfftsize * sizeof(fftw_real));

  report_memory_usage(&HighMark_run, "SMOOTH");

  asmth2 = Smoothing;
  asmth2 *= asmth2;

  fft_sGrid = (fftw_complex *) &sGrid[0];

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

		  ip = Ndim * (Ndim / 2 + 1) * (y - slabstart_y) + (Ndim / 2 + 1) * x + z;

		  win = exp(-k2 * asmth2 / 2);
		  Norm_local += win; 

		  fft_sGrid[ip].re = fft_Grid[ip].re * win * fac2 / fac;
		  fft_sGrid[ip].im = fft_Grid[ip].im * win * fac2 / fac;
   }
 
  MPI_Allreduce(&Norm_local, &Normalization, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  if(ThisTask == 0)
	printf("Normalization = %g\n",Normalization*fac2);

  if(ThisTask == 0)
      printf("FFT backwards\n");
  rfftwnd_mpi(fft_inverse_plan, 1, sGrid, workspace, FFTW_TRANSPOSED_ORDER);

}


