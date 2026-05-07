
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include <gsl/gsl_histogram.h>
#include <gsl/gsl_histogram2d.h>

#include "allvars.h"
#include "utils.h"

/*
gsl_histogram **init_histograms(double minx, double maxx)
{
	gsl_histogram **gsl_hist = (gsl_histogram **) my_malloc(sizeof(gsl_histogram*)*NSnaps);

	for ( int i = 0; i < NSnaps; i++ )
	{
		gsl_hist[i] = gsl_histogram_calloc (nbins);
		gsl_histogram_set_ranges_uniform (gsl_hist[i], minx, maxx);
	}
	return(gsl_hist);
}

gsl_histogram2d **init_histograms2d(double minx, double maxx, double miny, double maxy)
{
	gsl_histogram2d **gsl_hist = (gsl_histogram2d **) my_malloc(sizeof(gsl_histogram2d*)*NSnaps);

	for ( int i = 0; i < NSnaps; i++ )
	{
		gsl_hist[i] = gsl_histogram2d_calloc (nbins,nbins);
		gsl_histogram2d_set_ranges_uniform (gsl_hist[i], minx, maxx, miny, maxy );
	}
	return(gsl_hist);
}
*/
/*
gsl_histogram3d **init_histograms3d(double minx, double maxx, double miny, double maxy, double minz, double maxz)
{
	gsl_histogram3d **gsl_hist = (gsl_histogram3d **) my_malloc(sizeof(gsl_histogram3d*)*NSnaps);

	for ( int i = 0; i < NSnaps; i++ )
	{
		gsl_hist[i] = gsl_histogram3d_calloc (nbins,nbins,nbins);
		gsl_histogram3d_set_ranges_uniform (gsl_hist[i], minx, maxx, miny, maxy, minz, maxz );
	}
	return(gsl_hist);
}

gsl_histogram4d **init_histograms4d(double minx, double maxx, double miny, double maxy, double minz, double maxz, double mint, double maxt)
{
	gsl_histogram4d **gsl_hist = (gsl_histogram4d **) my_malloc(sizeof(gsl_histogram4d*)*NSnaps);

	for ( int i = 0; i < NSnaps; i++ )
	{
		gsl_hist[i] = gsl_histogram4d_calloc (nbins,nbins,nbins,nbins);
		gsl_histogram4d_set_ranges_uniform (gsl_hist[i], minx, maxx, miny, maxy, minz, maxz, mint, maxt );
	}
	return(gsl_hist);
}
*/
void save_histogram(gsl_histogram * h, char *filename)
{
  FILE *fd;
  int i;
  float *tmp, *tmpTot;
  char buf[255];

  tmp = my_malloc(sizeof(float) * h->n);
  tmpTot = my_malloc(sizeof(float) * h->n);
  for(i = 0; i < h->n; i++)
    tmp[i] = h->bin[i];

  MPI_Allreduce(tmp, tmpTot, h->n, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);

  for(i = 0; i < h->n; i++)
    h->bin[i] = tmpTot[i];

  if(ThisTask == 0)
    {
      fd = my_fopen(filename, "w");
      gsl_histogram_fprintf(fd, h, "%f", "%f");
      fclose(fd);
    }
  my_free(tmpTot);
  my_free(tmp);
}

void save_histogram2d(gsl_histogram2d * h, char *filename)
{
  FILE *fd;
  int i;
  float *tmp, *tmpTot;
  char buf[255];

  tmp = my_malloc(sizeof(float) * h->nx * h->ny);
  tmpTot = my_malloc(sizeof(float) * h->nx * h->ny);
  for(i = 0; i < h->nx * h->ny; i++)
    tmp[i] = h->bin[i];

  MPI_Allreduce(tmp, tmpTot, h->nx * h->ny, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);

  for(i = 0; i < h->nx * h->ny; i++)
    h->bin[i] = tmpTot[i];

  if(ThisTask == 0)
    {
      fd = my_fopen(filename, "w");
      gsl_histogram2d_fprintf(fd, h, "%f", "%f");
      fclose(fd);
    }
  my_free(tmpTot);
  my_free(tmp);

}

/*
void save_histograms3d(gsl_histogram3d *h, char *filename,int snap)
{
	FILE *fd;
	int i;
	float *tmp, *tmpTot;
	char buf[255];

	tmp = my_malloc(sizeof(float) * h->nx * h->ny * h->nz);
	tmpTot = my_malloc(sizeof(float) * h->nx * h->ny * h->nz);
	for (i = 0; i < h->nx * h->ny * h->nz; i++)
		tmp[i] = h->bin[i];

	MPI_Allreduce( tmp, tmpTot, h->nx * h->ny * h->nz, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);

	for (i = 0; i < h->nx * h->ny * h->nz; i++)
		h->bin[i] = tmpTot[i];

	if (ThisTask == 0 )
	{
		sprintf(buf,"%s/%s%s.%03d",base,filename,sufix,snap);
		fd = my_fopen(buf,"w");
		gsl_histogram3d_fprintf(fd, h, "%f","%f");
		fclose(fd);
	}
	my_free(tmpTot);
	my_free(tmp);

}

void save_histograms4d(gsl_histogram4d *h, char *filename, int snap)
{
	FILE *fd;
	int i;
	float *tmp, *tmpTot;
	char buf[255];

	tmp = my_malloc(sizeof(float) * h->nx * h->ny * h->nz * h->nt);
	tmpTot = my_malloc(sizeof(float) * h->nx * h->ny * h->nz * h->nt);
	for (i = 0; i < h->nx * h->ny * h->nz * h->nt; i++)
		tmp[i] = h->bin[i];

	MPI_Allreduce( tmp, tmpTot, h->nx * h->ny * h->nz * h->nt, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);

	for (i = 0; i < h->nx * h->ny * h->nz * h->nt; i++)
		h->bin[i] = tmpTot[i];

	if (ThisTask == 0 )
	{
		sprintf(buf,"%s/%s%s.%03d",base,filename,sufix,snap);
		fd = my_fopen(buf,"w");
		gsl_histogram4d_fprintf(fd, h, "%f","%f");
		fclose(fd);
	}
	my_free(tmpTot);
	my_free(tmp);

}
*/
