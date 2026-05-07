#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <string.h>

#include <srfftw_mpi.h>
#include <gsl/gsl_histogram.h>

#include "allvars.h"
#include "utils.h"
#include "correl.h"
#include "mymalloc.h"

#define ASMTH 1.25


#define dbg1 0
#define dbg2 0

#define ABSN(a) ((a > 0) ? a : -a )

//typedef struct result_data
//{
//  double x; 
//  double y;
//  double r;
//  double data;
//  long long n; 
//} result;

int wrapn(int ii, int ngrid)
{
  if(ii >= ngrid)
    ii -= ngrid;
  if(ii < 0)
    ii += ngrid;
  return ii;
}				//wrapn

void correl(char fname[], Grid grid, int nfold, float boxsize)
{
  int i, j, k, itask, jtask;
  unsigned long ibin;
  int i_ng, j_ng, k_ng, nn, jpos;
  int upperklim;
  int nbins;
  int ijk, id_tmp1, id_tmp2, id_ng, id;
  int ii, jj, kk;
  int ai, aj;
  float CellSize;
  MPI_Status status;
  GridFloat *grid_nb, *grid_nb_pos;

  int *xi_bin;
  double *distance;
  int *dist_ok, *j_given_i, *k_given_ij;
  FILE *fd;
  char buf[1000];

  double R0, R1, rbinfac, LOGR_R0, LOGR_R1, LOGR_rbinfac;
  int BINS_XI, LOGR_BINS_XI;

  result *Correl, *LogCorrel;

  CellSize = boxsize / (float) grid.Ngrid;

  BINS_XI = 70;
  R0 = 0.;
  R1 = 200;
  rbinfac = BINS_XI / (R1 - R0);

  LOGR_BINS_XI = 50;
  LOGR_R0 = log10(0.001);
  LOGR_R1 = log10(200.);
  LOGR_rbinfac = LOGR_BINS_XI / (LOGR_R1 - LOGR_R0);

  Correl = mymalloc_movable(&Correl,"Correl", BINS_XI * sizeof(result));
  //LogCorrel = mymalloc_movable(&LogCorrel,"LogCorrel", LOGR_BINS_XI * sizeof(result));

  memset(Correl, 0, BINS_XI * sizeof(result));
  //memset(LogCorrel, 0, BINS_XI * sizeof(result));

  nn = MIN(R1 / CellSize, 0.5 * (grid.Ngrid + 1));	//maximum number of lattice points between i and j in any direction

  xi_bin     = mymalloc("xi_bin",     CUB(nn + 1) * sizeof(int));
  distance   = mymalloc("distance",   CUB(nn + 1) * sizeof(double));
  dist_ok    = mymalloc("dist_ok",    CUB(nn + 1) * sizeof(int));
  j_given_i  = mymalloc("j_given_i",  (nn + 1) * sizeof(int));
  k_given_ij = mymalloc("k_given_ij", SQR(nn + 1) * sizeof(int));
  memset(dist_ok, 0, CUB(nn + 1) * sizeof(int));

  MSG0(" creating distance mask (%d)", t0, ThisTask, nn);

  for(i_ng = -nn; i_ng <= nn; i_ng++)
    {
      id_tmp1 = ABSN(i_ng);
      j_given_i[id_tmp1] = 0;
      for(j_ng = -nn; j_ng <= nn; j_ng++)
        {
          id_tmp2 = ABSN(j_ng) + (nn + 1) * id_tmp1;
          k_given_ij[id_tmp2] = 0;
          for(k_ng = -nn; k_ng <= nn; k_ng++)
            {
              id = ABSN(k_ng) + (nn + 1) * id_tmp2;
              distance[id] = CellSize * sqrt(SQR(i_ng) + SQR(j_ng) + SQR(k_ng));
              if(distance[id] < R1 && distance[id] > R0)
                  {
                      dist_ok[id] = 1;
                      ibin = (distance[id]-R0)*rbinfac;
                      xi_bin[id] = ibin;	//assign physical distance bin from the centre to all grid points in a block of size 2nn x 2nn x 2nn around
                      j_given_i[id_tmp1] = j_ng;
                      k_given_ij[id_tmp2] = k_ng;
                  }
            }	
        }
    }	

  int slabsNeed = 0, TaskNeed = 0, TasksNeedingMe = 0;
  int count_slabs = 0, count_slabs_send;
  while(count_slabs < nn)
    {
      TaskNeed++;
      count_slabs += grid.slabs_per_task[wrapn(ThisTask + TaskNeed, NTasks)];
    }

  for(itask = 0; itask < NTasks; itask++)
    {
      count_slabs_send = 0;
      for(jtask = 1; count_slabs_send < nn; jtask++)
	{
	  if(wrapn(itask + jtask, NTasks) == ThisTask)
	    {
	      TasksNeedingMe++;
	      break;
	    }			//if
	  count_slabs_send += grid.slabs_per_task[wrapn(itask + jtask, NTasks)];
	}			//for
    }				//for

  MSG0(" Allocating extra %fMb to hold extra grids (%dx%d)", t0, ThisTask,
       SQR(grid.Ngrid) * nn * sizeof(GridFloat) / 1024. / 1024., grid.Ngrid, grid.local_nx);

  GridFloat *grid_extra = mymalloc("grid_extra", SQR(grid.Ngrid) * count_slabs * sizeof(GridFloat));
  GridFloat *buffer = mymalloc("buffer", SQR(grid.Ngrid) * grid.local_nx * sizeof(GridFloat));

  double gc_x, gc_y, gc_z;
  double dx = 0.0, dy = 0.0, dz = 0.0;
  double maxdd = (boxsize / grid.Ngrid) * 0.5;
  double dim = grid.local_nx * SQR(grid.Ngrid);

  for(i = 0; i < grid.local_nx; i++)
      for(j = 0; j < grid.Ngrid; j++)
	  for(k = 0; k < grid.Ngrid; k++)
	    {
	      ijk = k + grid.Ngrid * (j + grid.Ngrid * i);
	      buffer[ijk] = grid.data[ijk];
	    }			//for

  MPI_Barrier(MPI_COMM_WORLD);
  if(NTasks == 1)
      for(i = 0; i < grid.local_nx * SQR(grid.Ngrid); i++)
	  grid_extra[i] = buffer[i];
  else
    {
      int start = 0;
      for(itask = 0; itask < NTasks; itask++)
	{
	  for(jtask = 1; jtask <= MAX(TaskNeed, TasksNeedingMe); jtask++)
	    {
	      if(ThisTask == itask && jtask <= TasksNeedingMe)
		{
		  DBG("Sending extra grid", dbg1, t0, ThisTask);
		    MPI_Send(buffer, sizeof(GridFloat) * SQR(grid.Ngrid) * grid.local_nx, MPI_BYTE,
			   wrapn(itask - jtask, NTasks), TAG_GRID, MPI_COMM_WORLD);
		}		//if
	      if(ThisTask == wrap(itask - jtask, NTasks) && jtask <= TaskNeed)
		{
		  DBG("Receiving extra grid", dbg1, t0, ThisTask);
		  for(start = 0, i = 1; i < jtask; i++)
		    start += grid.slabs_per_task[wrap(itask - jtask + i, NTasks)] * SQR(grid.Ngrid);
		  MPI_Recv(&grid_extra[start],
			   sizeof(GridFloat) * SQR(grid.Ngrid) * grid.slabs_per_task[itask], MPI_BYTE, itask,
			   TAG_GRID, MPI_COMM_WORLD, &status);
		}		//if
	    }			//for
	}			//for
    }				//else

  double dist_rr, xx;
  myfree(buffer);

  MSG0(" computing pairs", t0, ThisTask);
  for(i = 0; i < grid.local_nx; i++)
    {				//each processor has its own slab of data
      MSG0(" cell %d out of %d", t0, ThisTask, i, grid.local_nx - 1);
      for(j = 0; j < grid.Ngrid; j++)
        {
          for(k = 0; k < grid.Ngrid; k++)
            {
              ijk = k + grid.Ngrid * (j + grid.Ngrid * i);	//determine cell id, then loop over the nn x 2nn x 2nn block above it

              for(i_ng = 0; i_ng <= nn; i_ng++)
                {
                  //ii = wrapn(i_ng+i, Ngrid);
                  ii = i_ng + i;
                  grid_nb = grid.data;
                  if(ii >= grid.local_nx)
                    {
                      ii -= grid.local_nx;
                      grid_nb = &grid_extra[0];
                    }//if

                  ai = ABSN(i_ng);
                  for(j_ng = -j_given_i[ai]; j_ng <= j_given_i[ai]; j_ng++)
                    {
                      jj = wrapn(j_ng + j, grid.Ngrid);
                      if(jj >= 0 && jj < grid.Ngrid)
                        {
                          aj = ABSN(j_ng);
                          id_tmp1 = aj + (nn + 1) * ai;
                          upperklim = k_given_ij[id_tmp1];

                          if(i_ng == 0)
                            {
                              upperklim = 0;
                              if(j_ng >= 0)
                                upperklim = -1;
                            }//if
                         for(k_ng = -k_given_ij[id_tmp1]; k_ng <= upperklim; k_ng++)
                           {
                             id = ABSN(k_ng) + (nn + 1) * id_tmp1;
                             if(dist_ok[id])
                               {
                                 id_ng = wrapn(k_ng + k, grid.Ngrid) + grid.Ngrid * (jj + grid.Ngrid * ii);
                                             
                                 Correl[xi_bin[id]].data += grid.data[ijk] * grid_nb[id_ng];
                                 Correl[xi_bin[id]].n += 1.0;
                                 Correl[xi_bin[id]].r += distance[id];
                               }	
                           }	
                         }
                    }	
                }
            }	
        }
    }	

  gather_results(Correl, BINS_XI);

  if(ThisTask == 0)
    {
#ifdef PARTICLE_TYPE
      sprintf(buf, "%s_linear_xigrid.pt_%d.%d.%d.fold%d", fname, PARTICLE_TYPE, (int)grid.Ngrid, flag_redshift, nfold);
#else
      sprintf(buf, "%s_linear_xigrid.%d.%d.fold%d", fname, (int)grid.Ngrid, flag_redshift, nfold);
#endif
      printf("%s\n", buf);
      
      fd = my_fopen(buf, "w");
      //gsl_histogram_fprintf(fd, xi, "%e", "%e");
      for(i = 0; i < BINS_XI; i++)
	  fprintf(fd, "%e %e \n", Correl[i].r, Correl[i].data);
      fclose(fd);
      printf(" writing results to %s\n", buf);
    }//if

  myfree(grid_extra);
  myfree(k_given_ij);
  myfree(j_given_i);
  myfree(dist_ok);
  myfree(distance);
  myfree(xi_bin);

  //myfree(LogCorrel);
  myfree(Correl);
}


//void gather_results(result *in, int size, double norm)
// {        
//   int i, n;
//   long long *nbuf, *nin;
//   double *xbuf, *ybuf, *xin, *yin, *rin, *din, *dbuf, *d2in, *d2buf, *rbuf;
//          
//   nbuf = mymalloc("nbuf", NTasks * size * sizeof(long long));
//   xbuf = mymalloc("xbuf", NTasks * size * sizeof(double));
//   ybuf = mymalloc("ybuf", NTasks * size * sizeof(double));
//   rbuf = mymalloc("rbuf", NTasks * size * sizeof(double));
//   dbuf = mymalloc("dbuf", NTasks * size * sizeof(double));
//   d2buf = mymalloc("d2buf", NTasks * size * sizeof(double));
//          
//   nin = mymalloc("nin", NTasks * size * sizeof(long long));
//   xin = mymalloc("xin", NTasks * size * sizeof(double));
//   yin = mymalloc("yin", NTasks * size * sizeof(double));                                                                                                             
//   rin = mymalloc("rin", NTasks * size * sizeof(double));
//   din = mymalloc("din", NTasks * size * sizeof(double));
//   d2in = mymalloc("d2in", NTasks * size * sizeof(double));
//          
//   for(i = 0; i < size; i++)
//     {      
//         nin[i] = in[i].n;
//         xin[i] = in[i].x;
//         yin[i] = in[i].y;
//         rin[i] = in[i].r;
//         din[i] = in[i].data;
//         d2in[i] = in[i].data2;
//     }      
//          
//   MPI_Allgather(xin, size, MPI_DOUBLE,    xbuf, size, MPI_DOUBLE,    MPI_COMM_WORLD);
//   MPI_Allgather(yin, size, MPI_DOUBLE,    ybuf, size, MPI_DOUBLE,    MPI_COMM_WORLD);
//   MPI_Allgather(rin, size, MPI_DOUBLE,    rbuf, size, MPI_DOUBLE,    MPI_COMM_WORLD);
//   MPI_Allgather(din, size, MPI_DOUBLE,    dbuf, size, MPI_DOUBLE,    MPI_COMM_WORLD);
//   MPI_Allgather(d2in,size, MPI_DOUBLE,   d2buf, size, MPI_DOUBLE,    MPI_COMM_WORLD);
//   MPI_Allgather(nin, size, MPI_LONG_LONG, nbuf, size, MPI_LONG_LONG, MPI_COMM_WORLD);
//
//   for(i = 0; i < size; i++)
//     {  
//       in[i].n = in[i].x = in[i].y = in[i].r = in[i].data = 0;
//       for(n = 0; n < NTasks; n++)
//         {
//           in[i].x  += xbuf[n * size + i];
//           in[i].y  += ybuf[n * size + i];
//           in[i].r  += rbuf[n * size + i];
//           in[i].data  += dbuf[n * size + i];
//           in[i].data2  += d2buf[n * size + i];
//           in[i].n  += nbuf[n * size + i];
//         }
//     }
//      
//   for(i = 0; i < size; i++)
//     {  
//       if(in[i].n > 0)
//         {
//           in[i].x  = in[i].x / in[i].n;
//           in[i].y  = in[i].y / in[i].n;
//           in[i].r  = in[i].r / in[i].n;
//           in[i].data  = norm * in[i].data / in[i].n;
//           in[i].data2 = norm * in[i].data2 / in[i].n;
//        } 
//     }
//      
//   myfree(d2in);
//   myfree(din);
//   myfree(rin);
//   myfree(yin);
//   myfree(xin);
//   myfree(nin);
//      
//   myfree(d2buf);
//   myfree(dbuf);
//   myfree(rbuf);
//   myfree(ybuf);
//   myfree(xbuf);
//   myfree(nbuf);
//}
