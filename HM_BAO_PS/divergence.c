
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <mpi.h>

#include <srfftw_mpi.h>

#include "allvars.h"
#include "mymalloc.h"
//#include "gadget.h"
#include "utils.h"
#include "power.h"

#define ASMTH 1.25

//#define CUB(x) (x*x*x)
#define dbg1 0
#define dbg2 0

//#ifndef PARTICLE_TYPE
//#define PARTICLE_TYPE 1+2+4+8+16+32 
//#endif
//extern int Ngrid, Ngrid2;
extern time_t t0;
extern int ThisTask, NTasks;
extern double BoxSize;
//extern int *slab_to_task, *slab_to_task_local, *slabs_per_task, *first_slab_of_task;
double c_pecvel(double a, double Omega_M, double Omega_L, double H0);

double c_pecvel(double a, double Omega_M, double Omega_L, double H0)
{

  double Omega_c=1.-Omega_M-Omega_L;

  double H=H0*sqrt(Omega_M/a/a/a+Omega_L+Omega_c/a/a); // Beware only LCDM

  // Om=Om_M/(H/H0)^2/a^3
  double Omega=Omega_M/(((H/H0)*(H/H0))*(a*a*a));

  // f=dlnD_dlna
  double f=pow(Omega,5./9.);

  double out=f*H;

  return(out);

}

void divergence(Grid grid, Grid grid_vx, Grid grid_vy, Grid grid_vz, Grid * theta_grid, float boxsize)
{
  long long i, j, k, ii, jj, kk, x, y, z, zz;
  int bin, bin_per, bin_par, dim, ijk_fft;
  double fx, fy, fz, ff, po, kvec;
  long long ip;
  long l, m, n, nbin;
  long long ijk, le, me, ne;
  unsigned long *nx;
  long long id = 0, id2;
  double norm;
  double pk, k2, kx, ky, kz, deltak;
  double r2, rx, ry, rz, rper, rpar;
  double radius, radius_par, radius_per;
  char buf[255];
  double ktot, Sk, ak, bk, k_nyquist;
  double kpar, kper;
  double mode, mode_par, mode_per, fac, smth, Asmth;
  double K0, K1, binfac;
  double R0, R1, rbinfac;
  int rep;
  int BINS_PS;

  static fftw_complex *fft_of_rhogrid;
  static fftw_real *rhogrid, *forcegrid, *workspace, *workspace2;

  static fftw_complex *fft_of_div;
  static fftw_real *div;

  static rfftwnd_mpi_plan fft_forward_plan, fft_inverse_plan;
  static int fftsize, maxfftsize, dimprodmax;

  FILE *fd_fft;
  char buf_fft[255];
  float *data_fft;
  FILE *fd;


  int local_x_start, local_nx, local_y_start, local_ny, total_local_size;

  MSG0(" Creating plan for %d cells", t0, ThisTask, grid.Ngrid);

  fft_forward_plan = rfftw3d_mpi_create_plan(MPI_COMM_WORLD, grid.Ngrid, grid.Ngrid, grid.Ngrid, 
                                             FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE | FFTW_IN_PLACE | FFTW_THREADSAFE);

  fft_inverse_plan = rfftw3d_mpi_create_plan(MPI_COMM_WORLD, grid.Ngrid, grid.Ngrid, grid.Ngrid, 
                                             FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE | FFTW_IN_PLACE | FFTW_THREADSAFE);

  rfftwnd_mpi_local_sizes(fft_forward_plan, &local_nx, &local_x_start, &local_ny, &local_y_start, &fftsize);


  MPI_Allreduce(&fftsize, &maxfftsize, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  rhogrid = (fftw_real *) mymalloc_movable(&rhogrid,"rhogrid",fftsize * sizeof(fftw_real));
  workspace = (fftw_real *) mymalloc_movable(&workspace,"workspace",maxfftsize * sizeof(fftw_real));

  div = (fftw_real *) mymalloc_movable(&div,"div",fftsize * sizeof(fftw_real));
  workspace2 = (fftw_real *) mymalloc_movable(&workspace2,"workspace",maxfftsize * sizeof(fftw_real));
  fft_of_div = (fftw_complex *) &div[0];

  for(i = 0; i < maxfftsize; i++)	/* clear local density field */
    workspace2[i] = 0;
  for(i = 0; i < fftsize; i++)	/* clear local density field */
  {
    rhogrid[i] = 0;
    div[i] = 0.;
    fft_of_div[i].re = 0.;
    fft_of_div[i].im = 0.;
  }

  double NumPart = 0, totNumPart = 0;
  ExpFactor = 1.;
  Omega_m = 0.25; 
  // Om=Om_M/(H/H0)^2/a^3
  double omegaz = Omega_m * 1/CUB(ExpFactor)/ (Omega_m/CUB(ExpFactor) + (1 - Omega_m));
  // f=dlnD_dlna
  double f=pow(omegaz,5./9.);

  fac = 1/f;
 
  for(i = 0; i < grid.local_nx; i++)
      for(j = 0; j < grid.Ngrid; j++)
	  for(k = 0; k < grid.Ngrid; k++)
	    {
	      ijk = k + grid.Ngrid * (j + grid.Ngrid * i);
              if (grid.data[ijk] > 0)
              {
                grid_vx.data[ijk] = fac * grid_vx.data[ijk]/grid.data[ijk];  
                grid_vy.data[ijk] = fac * grid_vy.data[ijk]/grid.data[ijk];  
                grid_vz.data[ijk] = fac * grid_vz.data[ijk]/grid.data[ijk];  
              }
	    }

  for (dim = 0; dim < 3; dim++)
  {
    MSG0(" Computing", t0, ThisTask);
    if (dim == 0)  MSG0(" Kx", t0, ThisTask);
    if (dim == 1)  MSG0(" Ky", t0, ThisTask);
    if (dim == 2)  MSG0(" Kz", t0, ThisTask);
 
    for(i = 0; i < grid.local_nx; i++)
      for(j = 0; j < grid.Ngrid; j++)
	  for(k = 0; k < grid.Ngrid; k++)
	    {
	      ijk = k + grid.Ngrid * (j + grid.Ngrid * i);
              ijk_fft = (i * grid.Ngrid + j) * (2 * (grid.Ngrid / 2 + 1)) + k;
	      if (dim == 0 ) rhogrid[ijk_fft] = (fftw_real) grid_vx.data[ijk];
	      if (dim == 1 ) rhogrid[ijk_fft] = (fftw_real) grid_vy.data[ijk];
	      if (dim == 2 ) rhogrid[ijk_fft] = (fftw_real) grid_vz.data[ijk];
            }            

    fft_of_rhogrid = (fftw_complex *) & rhogrid[0];
    for(i = 0; i < fftsize; i++)	/* clear local density field */
      workspace[i] = 0;

    MSG0(" FFT", t0, ThisTask);
    rfftwnd_mpi(fft_forward_plan, 1, rhogrid, workspace, FFTW_TRANSPOSED_ORDER);
    MSG0(" FFT done", t0, ThisTask);

    for(y = grid.local_y_start; y < grid.local_y_start + grid.local_ny; y++)
     for(x = 0; x < grid.Ngrid; x++)
      for(z = 0; z < grid.Ngrid/2 + 1; z++)
      {
	  if(x > grid.Ngrid / 2)
	    kx = x - grid.Ngrid;
	  else
	    kx = x;
	  if(y > grid.Ngrid / 2)
	    ky = y - grid.Ngrid;
	  else
	    ky = y;
	  if(z > grid.Ngrid / 2)
	    kz = z - grid.Ngrid;
	  else
	    kz = z;

	  k2 = kx * kx + ky * ky + kz * kz;

	  if(k2 > 0)
	    {

	      ip = grid.Ngrid * (grid.Ngrid / 2 + 1) * (y - grid.local_y_start) + (grid.Ngrid / 2 + 1) * x + z;

	     // if(k2 < (grid.Ngrid / 2.0))
		{
                  switch (dim)
                    {
                    case 0:
                      kvec = kx;
                      break;
                    case 1:
                      kvec = ky;
                      break;
                    case 2:
                      kvec = kz;
                      break;
                    }

		  fft_of_div[ip].re += (fftw_real) -kvec * fft_of_rhogrid[ip].im * ( 2 * M_PI / boxsize);
		  fft_of_div[ip].im += (fftw_real)  kvec * fft_of_rhogrid[ip].re * ( 2 * M_PI / boxsize);

		}
	    }
	}

  }

  MSG0(" FFT Back", t0, ThisTask);
  rfftwnd_mpi(fft_inverse_plan, 1, div, workspace2, FFTW_TRANSPOSED_ORDER);
  MSG0(" FFT Back done", t0, ThisTask);

  for(i = 0; i < grid.local_nx; i++)
    for(j = 0; j < grid.Ngrid; j++)
      for(k = 0; k < grid.Ngrid; k++)
      {
         ijk = k + grid.Ngrid * (j + grid.Ngrid * i);
         ijk_fft = (i * grid.Ngrid + j) * (2 * (grid.Ngrid / 2 + 1)) + k;
	 theta_grid->data[ijk] = (float) div[ijk_fft]/CUB(grid.Ngrid);
      }            

  myfree_movable(workspace2);
  myfree_movable(workspace);
  myfree_movable(rhogrid);
  myfree_movable(div);

  rfftwnd_mpi_destroy_plan(fft_forward_plan);
  rfftwnd_mpi_destroy_plan(fft_inverse_plan);

}


