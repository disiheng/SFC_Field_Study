
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include <srfftw.h>

#include "nrutil.h"
#include "allvars.h"
#include "proto.h"
#include "utils.h"

#define ASMTH 1.25

extern int Ngrid, Ngrid2;
extern time_t t0;

void power(fftw_real * rho, fftw_complex * fft_of_rho, float *kmode, float *power_spectrum, long *nmode, int rec)
{
  int i, j, k, ii, jj, kk;
  int l, m, n, nbin;
  int ijk, le, me, ne;
  unsigned long *nx;
  long id = 0, id2;
  float *data, norm;
  double pk, k2, kx, ky, kz, deltak;
  char buf[255];
  float ktot, Sk, ak, bk, k_nyquist;
  float kpar, kper;
  double fac, asmth2, smth, Asmth;

  Asmth = ASMTH * BoxSize / Ngrid;
  asmth2 = (2 * M_PI) * Asmth / BoxSize;
  asmth2 *= asmth2;

  rfftwnd_plan fft_forward_plan, fft_reverse_plan;
  fftw_real scale = 1.0 / (Ngrid * Ngrid * Ngrid);

  float **pow2d, **k2d, **nk2d;
  float **k2dper, **k2dpar;
  int npow, nbinx, nbiny;
  npow = (int) NCell / 2;
  k2d = matrix(1, npow, 1, npow);
  pow2d = matrix(1, npow, 1, npow);
  nk2d = matrix(1, npow, 1, npow);
  k2dpar = matrix(1, npow, 1, npow);
  k2dper = matrix(1, npow, 1, npow);


  for(i = 0; i < npow; i++)
    {
      kmode[i] = 0.0;
      power_spectrum[i] = 0.0;
      nmode[i] = 0;
    }
  for(i = 1; i <= npow; i++)
    for(j = 1; j <= npow; j++)
      {
	k2d[i][j] = 0.;
	pow2d[i][j] = 0.;
	nk2d[i][j] = 0.;
	k2dpar[i][j] = 0.;
	k2dper[i][j] = 0.;
      }

  msg("- Creating FFT plans", t0);

  fft_forward_plan = rfftw3d_create_plan(Ngrid, Ngrid, Ngrid, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE | FFTW_IN_PLACE);

  fft_reverse_plan = rfftw3d_create_plan(Ngrid, Ngrid, Ngrid, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE | FFTW_IN_PLACE);

  // The output is an Ngrid x Ngrid x (Ngrid/2+1) array of fftw_complex
  msg("- FFT ", t0);
  rfftwnd_one_real_to_complex(fft_forward_plan, rho, fft_of_rho);

  msg("- Computing the power", t0);

  CellSize = BoxSize / NCell;
  k_nyquist = 2. * Pi / (2. * CellSize);
  deltak = 2. * Pi / BoxSize;	//k_nyquist/(NCell/2.); //

  //power_spectrum[0] = out[0]*out[0];  /* DC component */
  for(i = 0; i < Ngrid; i++)
    for(j = 0; j < Ngrid; j++)
      for(k = 0; k < Ngrid / 2 + 1; k++)	/* (k < Ngrid/2 rounded up) */
	{
	  //for( k=0; k<Ngrid; k++)  /* (k < Ngrid/2 rounded up) */

	  ii = i;
	  if(i > Ngrid / 2)
	    ii -= Ngrid;
	  jj = j;
	  if(j > Ngrid / 2)
	    jj -= Ngrid;
	  kk = k;
	  if(k > Ngrid / 2)
	    kk -= Ngrid;

	  kx = (double) ii *deltak;
	  ky = (double) jj *deltak;
	  kz = (double) kk *deltak;

	  ktot = sqrt(SQR(kx) + SQR(ky) + SQR(kz));
	  kpar = sqrt(SQR(kx) + SQR(ky));
	  kper = sqrt(SQR(kz));

	  nbin = (int) (ktot / deltak);
	  nbinx = (int) (kpar / deltak / 2.);
	  nbiny = (int) (kper / deltak / 2.);

	  //if (ktot < k_nyquist && ktot > deltak && kz < 2.*Hubble/(photoz*3.e5))
	  if(ktot < k_nyquist && ktot >= deltak)
	    {
	      //Sk   = wk_assgn(kx) * wk_assgn(ky) * wk_assgn(kz);
	      Sk = wk_assgn(ii) * wk_assgn(jj) * wk_assgn(kk);

	      kk = k;
	      if(k >= Ngrid / 2 + 1)
		kk = Ngrid - k;

	      ijk = kk + (Ngrid / 2 + 1) * (j + Ngrid * i);
	      pk = (SQR(fft_of_rho[ijk].re) + SQR(fft_of_rho[ijk].im)) / SQR(Sk);
	      power_spectrum[nbin] += pk;
	      kmode[nbin] += ktot;
	      nmode[nbin] += 1;

	      if(nbinx > npow || nbiny > npow)
		{
		  printf("nbx:%d  nby: %d npow:%d\n", nbinx, nbiny, npow);
		  exit(1);
		}
	      if(nbiny > 0 && nbinx > 0)
		{
		  k2d[nbinx][nbiny] += ktot;
		  k2dper[nbinx][nbiny] += kper;
		  k2dpar[nbinx][nbiny] += kpar;
		  pow2d[nbinx][nbiny] += pk;
		  nk2d[nbinx][nbiny] += 1.;
		}

	      // multiply with Green's function for the potential
	      // do deconvolution
	    }
	}

  msg(" - finished computing power", t0);


  rfftwnd_destroy_plan(fft_forward_plan);
  rfftwnd_destroy_plan(fft_reverse_plan);

  msg(" writing 2d Pk", t0);
  FILE *fd;
  sprintf(buf, "%s.%d.2d", OutFile, (int) NCell);
  fd = fopen(buf, "w");
  for(i = 1; i <= npow / 2 - 1; i++)
    for(j = 1; j <= npow / 2 - 1; j++)
      {
	if(nk2d[i][j] == 0)
	  nk2d[i][j] = 1;
	fprintf(fd, "%f %f %f %f\n", k2dper[i][j] / nk2d[i][j], k2dpar[i][j] / nk2d[i][j],
		(pow2d[i][j] / nk2d[i][j]) * CUB(BoxSize / SQR(NCell)), nk2d[i][j]);
      }
  fclose(fd);

  //if (N % 2 == 0) /* N is even */
  //     power_spectrum[N/2] = out[N/2]*out[N/2];  /* Nyquist freq. */
  //...

}


float wk_assgn(int id)
{
  float d, kd, wk;

  //kd = 0.5*k*CellSize;
  kd = Pi * id / NCell;
  wk = 1.;
  if(kd != 0.0)
    {
      switch (AssgnScheme)
	{
	case 1:
	  wk = (sin(kd) / (kd));
	  break;
	case 2:
	  wk = SQR(sin(kd) / (kd));
	  break;
	}
    }

  return wk;
}
