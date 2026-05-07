
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include <srfftw_mpi.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_histogram.h>

#include "allvars.h"
#include "utils.h"
#include "correl.h"
#include "histograms.h"


#define ASMTH 1.25

extern int Ngrid, Ngrid2;
extern time_t t0;
extern double BoxSize;

extern int local_x_start, local_nx, local_y_start, local_ny, total_local_size;
extern int *slab_to_task, *slab_to_task_local, *slabs_per_task, *first_slab_of_task;
extern double to_slab_fac;

static rfftwnd_mpi_plan fft_forward_plan, fft_inverse_plan;
static int fftsize, maxfftsize, dimprodmax;

static fftw_complex *fft_of_rhogrid;
static fftw_real *rhogrid, *forcegrid, *workspace;

#define dbg1 0
#define dbg2 0

#define ABSN(a) ((a > 0) ? a : -a )

int wrapn(int ii, int Ngrid)
{
  if(ii >= Ngrid)
    ii -= Ngrid;
  if(ii < 0)
    ii += Ngrid;
  return ii;
}

void xcorrel(char fname[], Grid grid, Grid grid2)
{
  long i, j, k;
  unsigned long ibin;
  long i_ng, j_ng, k_ng, nn;
  int nbins;
  long ijk, id_ng, id;
  int ii, jj, kk, ai, aj, ak;
  float CellSize, MaxDistance;
  char buf[255];
  MPI_Status status;
  float *grid_nb;

  int *xi_bin;
  double *distance;
  FILE *fd;


  CellSize = BoxSize / (float) grid.Ngrid;
  MaxDistance = 40;
  nn = MaxDistance / CellSize;

  nbins = 50;
  gsl_histogram *xi = gsl_histogram_calloc(nbins);
  gsl_histogram *npairs = gsl_histogram_calloc(nbins);
  gsl_histogram_set_ranges_uniform(xi, 0, MaxDistance * sqrt(3));
  gsl_histogram_set_ranges_uniform(npairs, 0, MaxDistance * sqrt(3));

  distance = my_malloc(CUB(nn + 1) * sizeof(double));
  xi_bin = my_malloc(CUB(nn + 1) * sizeof(int));

  MSG0(" creating distance mask (%d)", t0, ThisTask, nn);

  for(i_ng = -nn; i_ng <= nn; i_ng++)
    for(j_ng = -nn; j_ng <= nn; j_ng++)
      for(k_ng = -nn; k_ng <= nn; k_ng++)
	{
	  id = ABSN(k_ng) + (nn + 1) * (ABSN(j_ng) + (nn + 1) * ABSN(i_ng));
	  distance[id] = CellSize * sqrt(SQR(i_ng) + SQR(j_ng) + SQR(k_ng));
	  gsl_histogram_find(xi, distance[id], &ibin);
	  xi_bin[id] = ibin;
	}

  MSG0(" Allocating extra %fMb to hold extra grids", t0, ThisTask, SQR(grid.Ngrid) * nn * sizeof(float) / 1024. / 1024.);
  float *grid_extra = my_malloc(SQR(grid.Ngrid) * nn * sizeof(float));
  float *buffer = my_malloc(SQR(grid.Ngrid) * nn * sizeof(float));

  for(i = 0; i < nn; i++)
    for(j = 0; j < grid.Ngrid; j++)
      for(k = 0; k < grid.Ngrid; k++)
	{
	  ijk = k + grid.Ngrid * (j + grid.Ngrid * i);
	  buffer[ijk] = grid.data[ijk];
	}

  if(NTasks == 1)
    {
      for(i = 0; i < nn * SQR(grid.Ngrid); i++)
	grid_extra[i] = buffer[i];
    }
  else
    {
      for(i = 0; i < NTasks; i++)
	{
	  if(ThisTask == i)
	    {
	      DBG("Sending extra grid", dbg1, t0, ThisTask);
	      MPI_Send(buffer, sizeof(float) * SQR(grid.Ngrid) * nn, MPI_BYTE, wrapn(i - 1, NTasks), TAG_GRID, MPI_COMM_WORLD);
	    }
	  if(ThisTask == wrap(i - 1, NTasks))
	    {
	      DBG("Receiving extra grid", dbg1, t0, ThisTask);
	      MPI_Recv(grid_extra, sizeof(float) * SQR(grid.Ngrid) * nn, MPI_BYTE, i, TAG_GRID, MPI_COMM_WORLD, &status);
	    }
	}
    }

  free(buffer);


  MSG0(" computing pairs", t0, ThisTask);
  for(i = 0; i < local_nx; i++)
    {
      MSG0(" cell %d", t0, ThisTask, i);
      for(j = 0; j < grid.Ngrid; j++)
	for(k = 0; k < grid.Ngrid; k++)
	  {
	    ijk = k + grid.Ngrid * (j + grid.Ngrid * i);

	    for(i_ng = 0; i_ng <= nn; i_ng++)
	      {
		//ii = wrapn(i_ng+i, Ngrid);
		ii = i_ng + i;
		grid_nb = grid.data;
		if(ii >= local_nx)
		  {
		    ii -= local_nx;
		    grid_nb = &grid_extra[0];
		  }

		ai = ABSN(i_ng);
		for(j_ng = -nn; j_ng <= nn; j_ng++)
		  {
		    jj = wrapn(j_ng + j, grid.Ngrid);
		    aj = ABSN(j_ng);
		    for(k_ng = -nn; k_ng <= nn; k_ng++)
		      {
			id_ng = wrapn(k_ng + k, grid.Ngrid) + grid.Ngrid * (jj + grid.Ngrid * ii);
			id = ABSN(k_ng) + (nn + 1) * (aj + (nn + 1) * ai);

			xi->bin[xi_bin[id]] += (grid2.data[ijk] - grid2.m1) * (grid_nb[id_ng] - grid.m1);
			npairs->bin[xi_bin[id]] += 1.0;
		      }
		  }
	      }
	  }
    }

  double *pairsbuf, *xibuf;
  double *pairsTot, *xiTot;

  pairsbuf = my_malloc(xi->n * sizeof(double));
  xibuf = my_malloc(xi->n * sizeof(double));

  pairsTot = my_malloc(xi->n * sizeof(double));
  xiTot = my_malloc(xi->n * sizeof(double));

  for(i = 0; i < xi->n; i++)
    {
      xibuf[i] = xi->bin[i];
      pairsbuf[i] = npairs->bin[i];
    }

  MPI_Allreduce(xibuf, xiTot, xi->n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(pairsbuf, pairsTot, npairs->n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  for(i = 0; i < xi->n; i++)
    if(pairsTot[i] > 0)
      xi->bin[i] = 1000. * xiTot[i] / pairsTot[i] / grid.m1 / grid2.m1;

  MSG0(" writing results", t0, ThisTask);

  if(ThisTask == 0)
    {
      sprintf(buf, "%s_xcorrel.%d", fname, grid.Ngrid);
      fd = my_fopen(buf, "w");
      gsl_histogram_fprintf(fd, xi, "%f", "%f");
      fclose(fd);
    }

  free(pairsbuf);
  free(xibuf);
  free(xiTot);
  free(pairsTot);

  free(distance);
  free(xi_bin);
  gsl_histogram_free(xi);
  gsl_histogram_free(npairs);
  free(grid_extra);


  rfftwnd_mpi_destroy_plan(fft_forward_plan);
}
