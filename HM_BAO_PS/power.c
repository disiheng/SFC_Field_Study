
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <mpi.h>
#include <string.h>

#ifdef DOUBLEPRECISION_FFTW
#include     <drfftw_mpi.h>
#else
#include     <srfftw_mpi.h>
#endif

#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_sf_legendre.h>

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

/*
double bispectrum (double *qq, size_t dim, void *params)
{

  //Randomnly chooses the orientation of q1 on a sphere, adn the orientation of 
  //q2 on a circle centered at q1, with tetha fixed (3 angles). 
  
   gsl_ran_dir_3d (const gsl_rng * r, &x, &y, &z);
   gsl_ran_dir_2d (const gsl_rng * r, &x, &y);
}
*/

void overdensity(Grid * grid)
{
  int i;
  for(i = 0; i < grid->Ngrid * grid->Ngrid * grid->local_nx; i++)
     grid->data[i] = grid->data[i]/grid->m1 - 1;

  density_statistics(grid, "Overdensity", 0);
}


void power(char fname[], Grid grid, Grid grid_x, int nfold, float boxsize)
{
  long long i, j, k, ii, jj, kk, x, y, z, zz;
  int binr, bin, bin_per, bin_par, LOGR_bin, bin_angle, sign;
  double fx, fy, fz, ff, kpo,kpoi, kpo_x, kpo_r;
  long long ip;
  long l, m, n, nbin;
  long long ijk, le, me, ne;
  unsigned long *nx;
  long long id = 0, id2;
  double norm;
  double pk, k2, kx, ky, kz, deltak;
  double r2, rx, ry, rz, rper, rpar;
  double radius, radius_par, radius_per, angle;
  char buf[255];
  double ktot, Sk, ak, bk, k_nyquist;
  double kpar, kper;
  double mode, mode_par, mode_per, fac, fac2, smth, Asmth;
  double data_corr=0, data_corr_x=0, data_r = 0;
  double MU0, MU1, mufac;
  double K0, K1, binfac, bbinfac;
  double R0, R1, rbinfac, LOGR_R0, LOGR_R1, LOGR_rbinfac;
  double *ftmp;
  int rep;
  int BINS_PS, BINS_BS, BINS_XI, LOGR_BINS_XI, BINS_MU;

  static fftw_complex *fft_of_rhogrid;
  static fftw_real *rhogrid, *workspace;

#ifdef CROSS
  static fftw_complex *fft_of_rhogrid_x;
  static fftw_real *rhogrid_x, *workspace_x;
#ifdef RATIO
  static fftw_complex *fft_of_rhogrid_r;
  static fftw_real *rhogrid_r, *workspace_r;
#endif
#endif

  static rfftwnd_mpi_plan fft_forward_plan, fft_inverse_plan;
  static int fftsize, maxfftsize, dimprodmax;

  FILE *fd_fft;
  char buf_fft[255];
  FILE *fd;
  int local_x_start, local_nx, local_y_start, local_ny, total_local_size;

  result *Correl, *Power, *LogCorrel;
  result *Correl2D, *Power2D, *CorrelAngle2D; 
#ifdef CROSS
  result *CrossCorrel, *CrossPower;
  result *CrossCorrelAngle2D; 
#ifdef RATIO
  result *Ratio;
#endif
#endif

  result *Bispec, *P1, *P2, *P3, *Bispec2D;
  result *buffer, *buffer2d;

  K0 = 2 * M_PI / boxsize;	/* minimum k */
  BINS_PS = grid.Ngrid/16;
  K1 = BINS_PS * K0;		/* maximum k */
  //improve parallelisatio. For large grids, one needs lot of memory
  //BINS_PS = 100;
  binfac = BINS_PS / (K1 - K0);
  //binfac = BINS_PS / (K1);

  BINS_BS = 20;
  bbinfac = (BINS_BS-1)/ M_PI; 

  BINS_XI = 100;
  R0 = 0.;
  R1 = 200.;
  rbinfac = BINS_XI / (R1 - R0);

  BINS_MU = 20;
  MU0 = 0;
  MU1 = 1.;
  mufac = BINS_MU / (MU1 - MU0);

  LOGR_BINS_XI = 100;
  LOGR_R0 = log10(0.005);
  LOGR_R1 = log10(200.);
  LOGR_rbinfac = LOGR_BINS_XI / (LOGR_R1 - LOGR_R0);

  Power = mymalloc_movable(&Power,"Power", BINS_PS * sizeof(result));
#ifdef CROSS
  CrossPower = mymalloc_movable(&CrossPower,"CrossPower", BINS_PS * sizeof(result));
#endif
  Power2D = mymalloc_movable(&Power2D,"Power2D", SQR(BINS_PS) * sizeof(result));
  Bispec = mymalloc_movable(&Bispec,"Bispec", BINS_BS * sizeof(result));
  P1 = mymalloc_movable(&P1,"P1", 1 * sizeof(result));
  P2 = mymalloc_movable(&P2,"P2", 1 * sizeof(result));
  P3 = mymalloc_movable(&P3,"P3", BINS_BS * sizeof(result));
  Bispec2D = mymalloc_movable(&Bispec2D,"Bispec2D", SQR(BINS_BS) * sizeof(result));

  Correl = mymalloc_movable(&Correl,"Correl", BINS_XI * sizeof(result));
  LogCorrel = mymalloc_movable(&LogCorrel,"LogCorrel", LOGR_BINS_XI * sizeof(result));
#ifdef CROSS
  CrossCorrel = mymalloc_movable(&CrossCorrel,"CrossCorrel", BINS_XI * sizeof(result));
#endif

  Correl2D = mymalloc_movable(&Correl2D,"Correl2D", SQR(BINS_XI)*sizeof(result));
  CorrelAngle2D = mymalloc_movable(&CorrelAngle2D,"CorrelAngle2D", SQR(BINS_XI) * sizeof(result));

#ifdef CROSS
  CrossCorrelAngle2D = mymalloc_movable(&CrossCorrelAngle2D,"CrossCorrelAngle2D", SQR(BINS_XI) * sizeof(result));
#endif

#if defined(CROSS) && defined(RATIO)
  Ratio = mymalloc_movable(&Ratio,"Ratio", BINS_XI * sizeof(result));
#endif

  memset(Power,        0, BINS_PS * sizeof(result));
#ifdef CROSS
  memset(CrossPower,   0, BINS_PS * sizeof(result));
#endif
  memset(Power2D,      0, SQR(BINS_PS) * sizeof(result));
  memset(Bispec,       0, BINS_BS * sizeof(result));
  memset(P1,           0, 1 * sizeof(result));
  memset(P2,           0, 1 * sizeof(result));
  memset(P3,           0, BINS_BS * sizeof(result));
  memset(Bispec2D,     0, SQR(BINS_BS) * sizeof(result));
  memset(Correl,       0, BINS_XI * sizeof(result));
  memset(LogCorrel,    0, LOGR_BINS_XI * sizeof(result));
#ifdef CROSS
  memset(CrossCorrel,  0, BINS_XI * sizeof(result));
#endif
  memset(Correl2D,     0, SQR(BINS_XI) * sizeof(result));
  memset(CorrelAngle2D,0, BINS_MU*BINS_XI* sizeof(result));
#ifdef CROSS
  memset(CrossCorrelAngle2D,0, BINS_MU*BINS_XI* sizeof(result));
#ifdef RATIO
  memset(Ratio,        0, BINS_XI* sizeof(result));
#endif
#endif


  MSG0(" Creating plan for %d cells", t0, ThisTask, grid.Ngrid);

  fft_forward_plan = rfftw3d_mpi_create_plan(MPI_COMM_WORLD, grid.Ngrid, grid.Ngrid, grid.Ngrid, 
                                             FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE | FFTW_IN_PLACE | FFTW_THREADSAFE);

  fft_inverse_plan = rfftw3d_mpi_create_plan(MPI_COMM_WORLD, grid.Ngrid, grid.Ngrid, grid.Ngrid, 
                                             FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE | FFTW_IN_PLACE | FFTW_THREADSAFE);

  rfftwnd_mpi_local_sizes(fft_forward_plan, &local_nx, &local_x_start, &local_ny, &local_y_start, &fftsize);

  MPI_Allreduce(&fftsize, &maxfftsize, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  rhogrid = (fftw_real *) mymalloc_movable(&rhogrid,"rhogrid",fftsize * sizeof(fftw_real));
  workspace = (fftw_real *) mymalloc_movable(&workspace,"workspace",maxfftsize * sizeof(fftw_real));
  memset(rhogrid, 0, fftsize * sizeof(fftw_real));
  memset(workspace, 0 , maxfftsize * sizeof(fftw_real)); 
  fft_of_rhogrid = (fftw_complex *) &rhogrid[0];

#ifdef CROSS
  rhogrid_x = (fftw_real *) mymalloc_movable(&rhogrid_x, "rhogrid_x", fftsize * sizeof(fftw_real));
  workspace_x = (fftw_real *) mymalloc_movable(&workspace_x, "workspace_x", maxfftsize * sizeof(fftw_real));

  memset(rhogrid_x, 0, fftsize * sizeof(fftw_real));
  memset(workspace_x, 0 , maxfftsize * sizeof(fftw_real)); 
  fft_of_rhogrid_x = (fftw_complex *) &rhogrid_x[0];

#ifdef RATIO
    rhogrid_r = (fftw_real *) mymalloc_movable(&rhogrid_r, "rhogrid_r", fftsize * sizeof(fftw_real));
    workspace_r = (fftw_real *) mymalloc_movable(&workspace_r, "workspace_r", maxfftsize * sizeof(fftw_real));

    memset(rhogrid_r, 0, fftsize * sizeof(fftw_real));
    memset(workspace_r, 0 , maxfftsize * sizeof(fftw_real)); 
    fft_of_rhogrid_r = (fftw_complex *) &rhogrid_r[0];
#endif
#endif

  for( i = 0; i < grid.local_nx; i++)
    for(j = 0; j < grid.Ngrid; j++)
	  for(k = 0; k < grid.Ngrid; k++)
	    {
	      ijk = k + grid.Ngrid * (j + grid.Ngrid * i);
	      rhogrid[(i * grid.Ngrid + j) * (2 * (grid.Ngrid / 2 + 1)) + k] = (fftw_real) grid.data[ijk];
#ifdef CROSS
	      rhogrid_x[(i * grid.Ngrid + j) * (2 * (grid.Ngrid / 2 + 1)) + k] = (fftw_real) grid_x.data[ijk];
#endif
	    }


  MSG0(" FFT", t0, ThisTask);
  rfftwnd_mpi(fft_forward_plan, 1, rhogrid, workspace, FFTW_TRANSPOSED_ORDER);
#ifdef CROSS
  rfftwnd_mpi(fft_forward_plan, 1, rhogrid_x, workspace_x, FFTW_TRANSPOSED_ORDER);
#endif
  MSG0(" FFT done", t0, ThisTask);

#ifdef BISPECTRUM
  MSG0(" Bispectrum Mode Counting", t0, ThisTask);
  float kk0, kk1, kk2, kk3, alpha, delta_k, delta_k2, theta;
  int ip1,ip2,ip3, zz1, zz2, zz3;
  float a1,a2,a3,a4,b1,b2,b3,b4, mode1, mode2, mode3;
  float k1x,k2x,k3x, k1y,k2y,k3y, k1z,k2z,k3z, dx, dy, dz;
  int x1,x2,x3, y1,y2,y3, z1,z2,z3;
  float kp1, kp2, kp3;
  double arg;

  // choose a |k1|,|k2|,|k3|, and |delta_k| or equivalently, |k1|,|k2|, and theta
  
  delta_k = 0.0125/10;
  delta_k = M_PI/boxsize;
  delta_k2 = SQR(delta_k);
  kk1 = 0.05;
  alpha = 2; 
  kk2 = alpha*kk1;

  //for (i=0; i<8; i++) 
  //  printf("%d %f %d\n",i, modulus(i, 8), inv_modulus(modulus(i, 8), 8)); 
/*
  for(y1 = grid.local_y_start; y1 < grid.local_y_start + grid.local_ny; y1++)
    for(x1 = 0; x1 < grid.Ngrid; x1++)
     for(z1 = 0; z1 < grid.Ngrid; z1++)
    {

          k1x = modulus(x1, grid.Ngrid);
          k1y = modulus(y1, grid.Ngrid);
          k1z = modulus(z1, grid.Ngrid);

          ip1 = fftw_index(k1x,k1y,k1z, grid.Ngrid, grid.local_y_start, &sign);
          a3 = fft_of_rhogrid[ip1].re;
          b3 = sign * fft_of_rhogrid[ip1].im;

          ip2 = fftw_index(-k1x,-k1y,-k1z, grid.Ngrid, grid.local_y_start, &sign);
          a4 = fft_of_rhogrid[ip2].re;
          b4 = sign * fft_of_rhogrid[ip2].im;

//          if (a3 != a4 || b3 != -b4)
 //           printf("%d %d (%g,%g) (%g,%g)\n", ip1,ip2, a3,a4,b3,b4);
          //assert(a3 == a4);
          //assert(b3 == -b4);
     } 
*/
  MSG0(" Bispectrum Communication", t0, ThisTask);
  fftw_real *global_rhogrid = mymalloc("globlal_rhogrid", SQR(grid.Ngrid) * 2 * (grid.Ngrid/2 + 1) * sizeof(fftw_real)); 
  MPI_Allgather(rhogrid, 2 * (grid.Ngrid/2 + 1) * grid.Ngrid * grid.local_ny, MPI_FLOAT, 
          global_rhogrid, 2 * (grid.Ngrid/2+1) * grid.Ngrid * grid.local_ny, MPI_FLOAT, MPI_COMM_WORLD);
  fftw_complex *global_fft_of_rhogrid = (fftw_complex *) &global_rhogrid[0];

  MSG0(" Bispectrum Mode Counting", t0, ThisTask);

  for(y1 = 0; y1 < grid.Ngrid; y1++)
    for(x1 = 0; x1 < grid.Ngrid; x1++)
      for(z1 = 0; z1 < grid.Ngrid/2+1; z1++)
      {
          k1x = modulus(x1, grid.Ngrid);
          k1y = modulus(y1, grid.Ngrid);
          k1z = modulus(z1, grid.Ngrid);
          ip = fftw_index(k1x,k1y,k1z, grid.Ngrid, 0, &sign);
          
          global_fft_of_rhogrid[ip].re *= deconvolution_factor(k1x,k1y,k1z, grid.Ngrid);
          global_fft_of_rhogrid[ip].im *= deconvolution_factor(k1x,k1y,k1z, grid.Ngrid);
      }

  fac = 2 * M_PI / boxsize;
  fac2 = SQR(fac); 

  for(y1 = grid.local_y_start; y1 < grid.local_y_start + grid.local_ny; y1++)
    for(x1 = 0; x1 < grid.Ngrid; x1++)
      for(z1 = 0; z1 < grid.Ngrid; z1++)
      {
          k1x = modulus(x1, grid.Ngrid);
          k1y = modulus(y1, grid.Ngrid);
          k1z = modulus(z1, grid.Ngrid);
          mode1 = sqrt( SQR(k1x) + SQR(k1y) + SQR(k1z)) * fac;

          if (mode1 > (kk1 - delta_k) && mode1 < (kk1 + delta_k) )
          {
            ip = fftw_index(k1x,k1y,k1z, grid.Ngrid, 0, &sign);
            a1 = global_fft_of_rhogrid[ip].re;
            b1 = sign * global_fft_of_rhogrid[ip].im;
            kp1 = SQR(a1) + SQR(b1); 

            //for(y2 = grid.local_y_start; y2 < grid.local_y_start + grid.local_ny; y2++)
            for(y2 = 0; y2 < grid.Ngrid; y2++)
              for(x2 = 0; x2 < grid.Ngrid; x2++)
                for(z2 = 0; z2 < grid.Ngrid; z2++)
                {
                  k2x = modulus(x2, grid.Ngrid);
                  k2y = modulus(y2, grid.Ngrid);
                  k2z = modulus(z2, grid.Ngrid);
                  mode2 = sqrt( SQR(k2x) + SQR(k2y) + SQR(k2z)) * fac;

                  //if (abs(mode2/mode1 - alpha) < 0.05)
                  if (mode2 > (kk2 - delta_k) && mode2 < (kk2 + delta_k))
                  {
                    ip = fftw_index(k2x,k2y,k2z, grid.Ngrid, 0, &sign);
                    a2 = global_fft_of_rhogrid[ip].re;
                    b2 = sign * global_fft_of_rhogrid[ip].im;
                    kp2 = SQR(a2) + SQR(b2); 

                    arg = (k1x*k2x + k1y*k2y + k1z*k2z)/(mode1 * mode2) * fac2;
                    if (arg > 1.f) arg = 1.f;
                    if (arg < -1.f) arg = -1.f;
                    mode = acos(arg);
                    bin = bbinfac * mode; 
                 
 
                    k3x =  -(k2x + k1x);   
                    k3y =  -(k2y + k1y);   
                    k3z =  -(k2z + k1z);   
                    mode3 = sqrt( SQR(k3x) + SQR(k3y) + SQR(k3z)) * fac; 
             //     bin = mode3/0.2*BINS_BS;

                    ip = fftw_index(k3x,k3y,k3z, grid.Ngrid, 0, &sign);
                    a3 = global_fft_of_rhogrid[ip].re;
                    b3 = sign * global_fft_of_rhogrid[ip].im;
                    kp3 = SQR(a3) + SQR(b3); 

                    kpo  = a1*a2*a3 - a1*b2*b3 - b1*a2*b3 - b1*b2*a3; 
                    kpoi = a1*a2*b3 + a1*b2*a3 + b1*a2*a3 - b1*b2*b3;
                    if (z1 == 0 || z2 == 0 || k3z == 0)
                       continue;

                    //ip = fftw_index(-k3x,-k3y,-k3z, grid.Ngrid, grid.local_y_start, &sign);
                    //a4 = fft_of_rhogrid[ip].re;
                    //b4 = sign * fft_of_rhogrid[ip].im;

                    if ( bin >= 0 && bin < BINS_BS)
                    {
                      Bispec[bin].r += mode;
                      //if (bin == 2)     
                      //   printf("%f\n",kpo);
                      Bispec[bin].data += kpo;
                      Bispec[bin].data2 += kpoi;
                      Bispec[bin].n += 1;

                      P1[0].data += kp1; 
                      P2[0].data += kp2; 
                      P3[bin].data += kp3; 
                      P3[bin].r += mode3; 
                      P1[0].n += 1; 
                      P2[0].n += 1; 
                      P3[bin].n += 1; 
                    } 
                    else
                      printf("%f %d\n",arg, bin);
                  }
               }
           }
	}
  myfree(global_rhogrid);
#endif

  MSG0(" Mode Counting", t0, ThisTask);
  for(y = grid.local_y_start; y < grid.local_y_start + grid.local_ny; y++)
    for(x = 0; x < grid.Ngrid; x++)
      for(z = 0; z < grid.Ngrid / 2 + 1; z++)
      {
          zz = z;
          if(z >= grid.Ngrid / 2 + 1)
            zz = grid.Ngrid - z;

          kx = modulus(x, grid.Ngrid);
          ky = modulus(y, grid.Ngrid);
          kz = modulus(z, grid.Ngrid);

          k2 = kx * kx + ky * ky + kz * kz;
          kper =  kx * kx + ky * ky;
          kpar =  kz * kz;

          if(k2 > 0)
            {

              ip = fftw_index(kx,ky,kz, grid.Ngrid, grid.local_y_start, &sign);
#ifdef CROSS
#ifdef RATIO
              kpo_r = (fft_of_rhogrid_x[ip].re * fft_of_rhogrid_x[ip].re + fft_of_rhogrid_x[ip].im * fft_of_rhogrid_x[ip].im);
              fft_of_rhogrid_r[ip].re = (fftw_real) kpo_r;
              fft_of_rhogrid_r[ip].im = (fftw_real) 0.0;
#endif
              kpo_x = (fft_of_rhogrid[ip].re * fft_of_rhogrid_x[ip].re + fft_of_rhogrid[ip].im * fft_of_rhogrid_x[ip].im);
              fft_of_rhogrid_x[ip].re = (fftw_real) kpo_x;
              fft_of_rhogrid_x[ip].im = (fftw_real) 0.0;
#endif
              kpo = (fft_of_rhogrid[ip].re * fft_of_rhogrid[ip].re + fft_of_rhogrid[ip].im * fft_of_rhogrid[ip].im);
              fft_of_rhogrid[ip].re = (fftw_real) kpo;
              fft_of_rhogrid[ip].im = (fftw_real) 0.0;

              if(k2 < (grid.Ngrid)*(grid.Ngrid))
                {
                  // do deconvolution 
                  //smth = deconvolution_factor(kx,ky,kz, grid.Ngrid);
                  //kpo *=  smth;

                  mode = sqrt(k2) * 2 * M_PI / boxsize;
                  mode_per = sqrt(kper) * 2 * M_PI / boxsize;
                  mode_par = sqrt(kpar) * 2 * M_PI / boxsize;

                  bin     = floor((mode - K0) * binfac);
                  bin_per = floor((mode_per - K0) * binfac);
                  bin_par = floor((mode_par - K0) * binfac);
                  //bin = log(mode / K0) * binfac;

                  if(bin >= 0 && bin < BINS_PS)
                    {
                      Power[bin].r += mode;
                      Power[bin].data += kpo;
                      Power[bin].n += 1;
                    }

#ifdef CROSS
                  if(bin >= 0 && bin < BINS_PS)
                    {
                      CrossPower[bin].r += mode;
                      CrossPower[bin].data += kpo_x;
                      CrossPower[bin].n += 1;
                    }
#endif

                  if(bin_par >= 0 && bin_par < BINS_PS && bin_per >= 0 && bin_per < BINS_PS)
                    {
                      bin = bin_par + BINS_PS*bin_per;
                      Power2D[bin].x += mode_par;
                      Power2D[bin].y += mode_per;
                      Power2D[bin].r += mode;
                      Power2D[bin].data += kpo;
                      Power2D[bin].n += 1;
                    }
                }
            }
	}

  if(grid.local_y_start == 0)
    fft_of_rhogrid[0].re = fft_of_rhogrid[0].im = 0.0;

  MSG0(" FFT Back", t0, ThisTask);
  rfftwnd_mpi(fft_inverse_plan, 1, rhogrid, workspace, FFTW_TRANSPOSED_ORDER);
  MSG0(" FFT Back done", t0, ThisTask);

#ifdef CROSS
  if(grid_x.local_y_start == 0)
    fft_of_rhogrid_x[0].re = fft_of_rhogrid_x[0].im = 0.0;
  MSG0(" FFT Back", t0, ThisTask);
  rfftwnd_mpi(fft_inverse_plan, 1, rhogrid_x, workspace_x, FFTW_TRANSPOSED_ORDER);
  MSG0(" FFT Back done", t0, ThisTask);
#ifdef RATIO
  if(grid_x.local_y_start == 0)
    fft_of_rhogrid_r[0].re = fft_of_rhogrid_r[0].im = 0.0;
  MSG0(" FFT Back", t0, ThisTask);
  rfftwnd_mpi(fft_inverse_plan, 1, rhogrid_r, workspace_r, FFTW_TRANSPOSED_ORDER);
  MSG0(" FFT Back done", t0, ThisTask);
#endif
#endif

  for(x = local_x_start; x < (local_x_start+grid.local_nx); x++)
    for(y = 0; y < grid.Ngrid; y++)
	  for(z = 0; z < grid.Ngrid; z++)
	  {

	    if(x > grid.Ngrid / 2)
	      rx = x - grid.Ngrid;
        else
          rx = x;

        if(y > grid.Ngrid / 2)
          ry = y - grid.Ngrid;
        else
          ry = y;

        if(z > grid.Ngrid / 2)
          rz = z - grid.Ngrid;
        else
          rz = z;

/*
          rx = x;
          ry = y;
          rz = z;
*/

        r2   = SQR(rx) + SQR(ry) + SQR(rz);
        rper = SQR(rx) + SQR(ry);
        rpar = SQR(rz);

        ijk = z + grid.Ngrid * (y + grid.Ngrid * (x-local_x_start));
        data_corr = rhogrid[((x-local_x_start) * grid.Ngrid + y) * (2 * (grid.Ngrid / 2 + 1)) + z];
#ifdef CROSS
        data_corr_x = rhogrid_x[((x-local_x_start) * grid.Ngrid + y) * (2 * (grid.Ngrid / 2 + 1)) + z];
#ifdef RATIO
        data_r = rhogrid_r[((x-local_x_start) * grid.Ngrid + y) * (2 * (grid.Ngrid / 2 + 1)) + z];
        data_r = data_corr_x/data_r;
#endif
#endif

        if(r2 > 0)
          {
            radius = sqrt(r2) *  boxsize/grid.Ngrid;
            radius_per = sqrt(rper) * boxsize/grid.Ngrid;
            radius_par = sqrt(rpar) * boxsize/grid.Ngrid;
            angle = radius_par/radius;

            binr = (radius - R0) * rbinfac;
            LOGR_bin = (log10(radius) - LOGR_R0) * LOGR_rbinfac;
            bin_per = (radius_per - R0) * rbinfac;
            bin_par = (radius_par - R0) * rbinfac;
            bin_angle = (angle - MU0) * mufac;

            if(binr >= 0 && binr < BINS_XI)
              {
                bin = binr;
                Correl[bin].r += radius;
                Correl[bin].data += data_corr;
                Correl[bin].n += 1;
              }

#ifdef CROSS
#ifdef RATIO
            if(binr >= 0 && binr < BINS_XI)
              {
                bin = binr;
                Ratio[bin].r += radius;
                Ratio[bin].data += data_r;
                Ratio[bin].n += 1;
              }
#endif
            if(binr >= 0 && binr < BINS_XI)
              {
                bin = binr;
                CrossCorrel[bin].r += radius;
                CrossCorrel[bin].data += data_corr_x;
                CrossCorrel[bin].n += 1;
              }

            if(bin_angle >= 0 && bin_angle < BINS_MU && binr >= 0 && binr < BINS_XI)
              {
                bin = bin_angle + BINS_MU * binr;
                CrossCorrelAngle2D[bin].x += radius;
                CrossCorrelAngle2D[bin].y += angle;
                CrossCorrelAngle2D[bin].r += radius;
                CrossCorrelAngle2D[bin].data += data_corr_x;
                CrossCorrelAngle2D[bin].n += 1;
              }
#endif
            if(LOGR_bin >= 0 && LOGR_bin < LOGR_BINS_XI)
              {
                bin = LOGR_bin;
                LogCorrel[bin].r += radius;
                LogCorrel[bin].data += data_corr;
                LogCorrel[bin].n += 1;
              }

            if(bin_par >= 0 && bin_par < BINS_XI && bin_per >= 0 && bin_per < BINS_XI)
 		      {
                bin = bin_par + BINS_XI * bin_per;
		        Correl2D[bin].x += radius_per;
		        Correl2D[bin].y += radius_par;
		        Correl2D[bin].r += radius;
		        Correl2D[bin].data += data_corr;
		        Correl2D[bin].n += 1;
		      }

            if(bin_angle >= 0 && bin_angle < BINS_MU && binr >= 0 && binr < BINS_XI)
 		      {
                bin = bin_angle + BINS_MU * binr;
		        CorrelAngle2D[bin].x += radius;
		        CorrelAngle2D[bin].y += angle;
		        CorrelAngle2D[bin].r += radius;
		        CorrelAngle2D[bin].data += data_corr;
		        CorrelAngle2D[bin].n += 1;
		      }
          }
      }

#ifdef RATIO
  myfree_movable(workspace_r);
  myfree_movable(rhogrid_r);
#endif

#ifdef CROSS
  myfree_movable(workspace_x);
  myfree_movable(rhogrid_x);
#endif

  myfree_movable(workspace);
  myfree_movable(rhogrid);

  if(ThisTask == 0)
    printf(" binning results\n");

  //norm = 4 * M_PI * CUB(Power[i].r) * CUB(boxsize / 2 / M_PI) / CUB(SQR((double) grid.Ngrid)) ;
  norm = CUB(boxsize / 2 / M_PI) / CUB(SQR((double) grid.Ngrid)) ;
  gather_results(Power, BINS_PS, norm);
#ifdef CROSS
  gather_results(CrossPower, BINS_PS, norm);
#endif
  gather_results(Power2D, SQR(BINS_PS), norm);

#ifdef BISPECTRUM
  norm = CUB(boxsize / 2 / M_PI) / SQR(CUB((double) grid.Ngrid)) ;
  gather_results(P1, 1, norm);
  gather_results(P2, 1, norm);
  gather_results(P3, BINS_BS, norm);
   //check normalization
  norm = SQR(CUB(boxsize / 2 / M_PI)) / CUB(CUB((double) grid.Ngrid));
  gather_results(Bispec, BINS_BS, norm);
  gather_results(Bispec2D, SQR(BINS_BS), norm);
#endif

  norm = CUB(BoxSize/boxsize) / CUB(SQR((double) grid.Ngrid)); 
  gather_results(Correl, BINS_XI, norm);
#ifdef CROSS
  gather_results(CrossCorrel, BINS_XI, norm);
#endif
  gather_results(LogCorrel, LOGR_BINS_XI, norm);
  gather_results(Correl2D, SQR(BINS_XI), norm);
  gather_results(CorrelAngle2D, BINS_XI*BINS_MU, norm);
#ifdef CROSS
  gather_results(CrossCorrelAngle2D, BINS_XI*BINS_MU, norm);
#ifdef RATIO
  gather_results(Ratio, BINS_XI, 1.0);
#endif
#endif

//  for (i=0;i<BINS_XI; i++)
//    Ratio[i].data = CrossCorrel[i].data/Ratio[i].data;

  sprintf(buf, "%s_pk.%d.%d.fold%d", fname, (int)grid.Ngrid, flag_redshift, nfold);
  write_results(Power,      BINS_PS, 1,       0, buf);

#ifdef BISPECTRUM
  for (i = 0; i < BINS_BS; i++)
     P3[i].data = P1[0].data*P2[0].data + P2[0].data*P3[i].data + P1[0].data*P3[i].data;

  printf("%f %f\n",P1[0].data,P2[0].data);
  sprintf(buf, "%s_redbispec.%d.%d.fold%d", fname, (int)grid.Ngrid, flag_redshift, nfold);
  write_results(P3,      BINS_BS, 1,       0, buf);
  sprintf(buf, "%s_bispec.%d.%d.fold%d", fname, (int)grid.Ngrid, flag_redshift, nfold);
  write_results(Bispec,      BINS_BS, 1,       0, buf);
  sprintf(buf, "%s_bispec2d.%d.%d.fold%d", fname, (int)grid.Ngrid, flag_redshift, nfold);
  write_results(Bispec2D,      BINS_BS, BINS_BS,       1, buf);
#endif

#ifdef CROSS
  sprintf(buf, "%s_crosspk.%d.%d.fold%d", fname, (int)grid.Ngrid, flag_redshift, nfold);
  write_results(CrossPower, BINS_PS, 1,       0, buf);
#ifdef RATIO
  sprintf(buf, "%s_ratio.%d.%d.fold%d", fname, (int)grid.Ngrid, flag_redshift, nfold);
  write_results(Ratio, BINS_XI, 1,       0, buf);
#endif
#endif

  sprintf(buf, "%s_linear_pk2d.%d.%d.fold%d", fname, (int)grid.Ngrid, flag_redshift, nfold);
  write_results(Power2D,    BINS_PS, BINS_PS, 1, buf);

  sprintf(buf, "%s_linear_xi.%d.%d.fold%d", fname, (int)grid.Ngrid, flag_redshift, nfold);
  write_results(Correl,       BINS_XI, 1, 0, buf);

#ifdef CROSS
  sprintf(buf, "%s_linear_crossxi.%d.%d.fold%d", fname, (int)grid.Ngrid, flag_redshift, nfold);
  write_results(CrossCorrel,  BINS_XI, 1, 0, buf);
#endif

  sprintf(buf, "%s_log_xi.%d.%d.fold%d", fname, (int)grid.Ngrid, flag_redshift, nfold);
  write_results(LogCorrel,    LOGR_BINS_XI, 1, 0, buf);
  sprintf(buf, "%s_linear_xi2d.%d.%d.fold%d", fname, (int)grid.Ngrid, flag_redshift, nfold);
  write_results(Correl2D,     BINS_XI, BINS_XI, 1, buf);
  sprintf(buf, "%s_linear_xiangle2d.%d.%d.fold%d", fname, (int)grid.Ngrid, flag_redshift, nfold);
  write_results(CorrelAngle2D, BINS_XI, BINS_MU, 1, buf);

  if(ThisTask == 0)
    {
      sprintf(buf, "%s_linear_leg.%d.%d.fold%d", fname, (int)grid.Ngrid, flag_redshift, nfold);
      double *Pleg; 
      Pleg = mymalloc("Pleg", 5 * sizeof(double) * BINS_XI);
      memset(Pleg,0, 5 * BINS_XI * sizeof(double));

      for (l = 0; l < 5 ; l++) 
        for (i = 0; i < BINS_XI; i++)
          for (j = 0; j < BINS_MU; j++)
          {
            bin = j + BINS_MU * i;
            Pleg[i + BINS_XI * l ] += (2 * l + 1) / mufac * CorrelAngle2D[bin].data*gsl_sf_legendre_Pl(l, CorrelAngle2D[bin].y) ;
          }

      printf("%s\n", buf);

      fd = my_fopen(buf, "w");
      for(i = 0; i < BINS_XI; i++)
        fprintf(fd, "%e %e %e %e\n", Correl[i].r, Pleg[i + BINS_XI * 0],  Pleg[i + BINS_XI * 2],  Pleg[i + BINS_XI * 4]);
      fclose(fd);

      myfree(Pleg);

#ifdef CROSS
      sprintf(buf, "%s_linear_crossleg.%d.%d.fold%d", fname, (int)grid.Ngrid, flag_redshift, nfold);
      Pleg = mymalloc("Pleg", 5 * sizeof(double) * BINS_XI);
      memset(Pleg,0, 5 * BINS_XI * sizeof(double));
 
      for (l = 0; l < 5 ; l++) 
        for (i = 0; i < BINS_XI; i++)
          for (j = 0; j < BINS_MU; j++)
          {
            bin = j + BINS_MU * i;
            Pleg[i + BINS_XI * l ] += (2 * l + 1) / mufac * CrossCorrelAngle2D[bin].data * gsl_sf_legendre_Pl(l, CrossCorrelAngle2D[bin].y) ;
          }

      printf("%s\n", buf);

      fd = my_fopen(buf, "w");
      for(i = 0; i < BINS_XI; i++)
	    fprintf(fd, "%e %e %e %e\n", Correl[i].r, Pleg[i + BINS_XI * 0],  Pleg[i + BINS_XI * 2],  Pleg[i + BINS_XI * 4]);
      fclose(fd);

      myfree(Pleg);
#endif
    }

#if defined(CROSS) && defined(RATIO)
  myfree(Ratio);
#endif
#ifdef CROSS
  myfree(CrossCorrelAngle2D);
#endif
  myfree(CorrelAngle2D);
  myfree(Correl2D);
#ifdef CROSS
  myfree(CrossCorrel);
#endif
  myfree(LogCorrel);
  myfree(Correl);
  myfree(Bispec2D);
  myfree(P3);
  myfree(P2);
  myfree(P1);
  myfree(Bispec);
  myfree(Power2D);
#ifdef CROSS
  myfree(CrossPower);
#endif
  myfree(Power);

  rfftwnd_mpi_destroy_plan(fft_forward_plan);
  rfftwnd_mpi_destroy_plan(fft_inverse_plan);
}


void gather_results(result *in, int size, double norm)
{
  int i, n;
  long long *nbuf, *nin;
  double *xbuf, *ybuf, *xin, *yin, *rin, *din, *dbuf, *d2in, *d2buf, *rbuf;

  nbuf = mymalloc("nbuf", NTasks * size * sizeof(long long));
  xbuf = mymalloc("xbuf", NTasks * size * sizeof(double));
  ybuf = mymalloc("ybuf", NTasks * size * sizeof(double));
  rbuf = mymalloc("rbuf", NTasks * size * sizeof(double));
  dbuf = mymalloc("dbuf", NTasks * size * sizeof(double));
  d2buf = mymalloc("d2buf", NTasks * size * sizeof(double));

  nin = mymalloc("nin", NTasks * size * sizeof(long long));
  xin = mymalloc("xin", NTasks * size * sizeof(double));
  yin = mymalloc("yin", NTasks * size * sizeof(double));
  rin = mymalloc("rin", NTasks * size * sizeof(double));
  din = mymalloc("din", NTasks * size * sizeof(double));
  d2in = mymalloc("d2in", NTasks * size * sizeof(double));

  for(i = 0; i < size; i++)
  {
      nin[i] = in[i].n;
      xin[i] = in[i].x;
      yin[i] = in[i].y;
      rin[i] = in[i].r;
      din[i] = in[i].data;
      d2in[i] = in[i].data2;
  }

  MPI_Allgather(xin, size, MPI_DOUBLE,    xbuf, size, MPI_DOUBLE,    MPI_COMM_WORLD);
  MPI_Allgather(yin, size, MPI_DOUBLE,    ybuf, size, MPI_DOUBLE,    MPI_COMM_WORLD);
  MPI_Allgather(rin, size, MPI_DOUBLE,    rbuf, size, MPI_DOUBLE,    MPI_COMM_WORLD);
  MPI_Allgather(din, size, MPI_DOUBLE,    dbuf, size, MPI_DOUBLE,    MPI_COMM_WORLD);
  MPI_Allgather(d2in,size, MPI_DOUBLE,   d2buf, size, MPI_DOUBLE,    MPI_COMM_WORLD);
  MPI_Allgather(nin, size, MPI_LONG_LONG, nbuf, size, MPI_LONG_LONG, MPI_COMM_WORLD);

  for(i = 0; i < size; i++)
    {
      in[i].n = in[i].x = in[i].y = in[i].r = in[i].data = 0;
      for(n = 0; n < NTasks; n++)
      {
        in[i].x  += xbuf[n * size + i];
        in[i].y  += ybuf[n * size + i];
        in[i].r  += rbuf[n * size + i];
        in[i].data  += dbuf[n * size + i];
        in[i].data2  += d2buf[n * size + i];
        in[i].n  += nbuf[n * size + i];
      }
    }

  for(i = 0; i < size; i++)
    {
      if(in[i].n > 0)
      {
        in[i].x  = in[i].x / in[i].n;
        in[i].y  = in[i].y / in[i].n;
        in[i].r  = in[i].r / in[i].n;
        in[i].data  = norm * in[i].data / in[i].n;
        in[i].data2 = norm * in[i].data2 / in[i].n;
      } 
    }

  myfree(d2in);
  myfree(din);
  myfree(rin);
  myfree(yin);
  myfree(xin);
  myfree(nin);

  myfree(d2buf);
  myfree(dbuf);
  myfree(rbuf);
  myfree(ybuf);
  myfree(xbuf);
  myfree(nbuf);
}  


void write_results(result *in, int sizex, int sizey, int type, char buf[])
{
  int i, size, nbins_x, nbins_y;
  double *ftmp;
  FILE *fd;
  
  size = sizex*sizey; 
  if(ThisTask == 0)
    {
      printf("%s\n", buf);
      fd = my_fopen(buf, "w");

      if (type == 0)
      {
        for(i = 0; i < size; i++)
	  fprintf(fd, "%g %g %g %d\n", in[i].r, in[i].data, in[i].data2, in[i].n);
      }

      if (type == 1)
      {
        ftmp = mymalloc("ftmp",size * sizeof(double));
        fwrite(&sizex, sizeof(int),1,fd);
        fwrite(&sizey, sizeof(int),1,fd);
        for (i = 0; i < size; i++) ftmp[i] = in[i].x;
        fwrite(ftmp, size * sizeof(double),1,fd);
        for (i = 0; i < size; i++) ftmp[i] = in[i].y;
        fwrite(ftmp, size * sizeof(double),1,fd);
        for (i = 0; i < size; i++) ftmp[i] = in[i].r;
        fwrite(ftmp, size * sizeof(double),1,fd);
        for (i = 0; i < size; i++) ftmp[i] = in[i].data;
        fwrite(ftmp, size * sizeof(double),1,fd);
        for (i = 0; i < size; i++) ftmp[i] = in[i].n;
        fwrite(ftmp, size * sizeof(double),1,fd);
        myfree(ftmp);
      }
      fclose(fd);
   }
}

int fftw_index(double kx, double ky, double kz, int ngrid, int local_y_start, int *sign)
{
  int x, y, z;
  int ip;

  if (kz < 0) 
  {
     x = inv_modulus(-kx,ngrid);
     y = inv_modulus(-ky,ngrid);
     z = inv_modulus(-kz,ngrid);
     *sign = -1;
  }
  else
  {
     x = inv_modulus(kx,ngrid);
     y = inv_modulus(ky,ngrid);
     z = inv_modulus(kz,ngrid);
     *sign = 1;
   }

  if(z >= ngrid / 2 + 1)
    z = ngrid - z;
  ip = (ngrid / 2 + 1) * ( ngrid * (y - local_y_start) + x) + z;
  return ip;
}


int inv_modulus(double kx, int ngrid)
{
  int x;
  if (kx < 0)
    x = ngrid + kx;
  else
    x = kx;
  return x;
}


double modulus(int x, int ngrid)
{
  double kx;
  if(x > ngrid / 2)
     kx = x - ngrid;
   else
    kx = x;
  return kx;
}

double deconvolution_factor(double kx, double ky, double kz, int ngrid)
{
    double smth, ff, fx, fy, fz; 
    fx = fy = fz = 1;
    if(kx != 0)
    {
      fx = (M_PI * kx) / ngrid;
      fx = sin(fx) / fx;
    }
    if(ky != 0)
    {
      fy = (M_PI * ky) / ngrid;
      fy = sin(fy) / fy;
    }
    if(kz != 0)
    {
      fz = (M_PI * kz) / ngrid;
      fz = sin(fz) / fz;
    }
    ff = 1 / (fx * fy * fz);
#ifdef CIC                 
    ff *= ff;
#endif 
		  // end deconvolution 
   return ff;
}


