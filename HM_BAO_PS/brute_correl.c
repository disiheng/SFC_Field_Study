#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <srfftw_mpi.h>
#include <gsl/gsl_histogram.h>

#include "allvars.h"
#include "../utils/utils.h"
#include "brute_correl.h"
#include "../utils/mymalloc.h"

#ifdef SUBREGION
#define ABSN(a) ((a > 0) ? a : -a )
#endif

float myperiodic (float x) {
  float xx;
  xx = x;
  if (x > 0.5*BoxSize)
    xx = BoxSize - x; 
  if (x < -0.5*BoxSize)
    xx = x + BoxSize; 
  return xx;
} //myperiodic

#ifdef GRIDSHELLS
float gridconv (float x) {
  return floor(x*InvCellSize+0.5)*CellSize;
} //gridconv
#endif

#ifdef SUBREGION
int inside_sub(float x, float y, float z) {
  if (ABSN(x-SubCenX)<0.5*BoxSize && ABSN(y-SubCenY)<0.5*BoxSize && ABSN(z-SubCenZ)<0.5*BoxSize) return 1;
  else return 0;
} //inside_sub
#endif

void brute_correl(char fname[], long Npart, Particle *P) {
  long i,j,jj,k;
  int ibin, firsti, lasti;
  float slabsize;
  float dist_sq;
  float rlow, rhigh;
  double *DDbuf, *DD;
  double MassTot=0.0;
  double RR;
  FILE *fd;
  char buf[1000];
#ifdef GRIDSHELLS
  int *shellcounts;
  CellSize=BoxSize/NCell;
  InvCellSize=NCell/BoxSize;
#endif
  //float MaxDistance=MIN(150.,0.5*BoxSize*(NCell+1.0)/NCell); //conforms with correl's MaxDistance
  float MaxDistance=150.;
#ifdef SMALLBRUTE
  float ValidDistance=SMALLBRUTE;
#else
  float ValidDistance=MIN(MaxDistance,0.5*BoxSize);
#endif
#ifndef LOGBINS
  float MinDistance=0.0;
  float InvRange=1.0/(MaxDistance-MinDistance);
#else
  float MinDistance=0.01;
  float LogMinDistance=log10(MinDistance);
  float LogMaxDistance=log10(MaxDistance);
  float InvLogRange=1.0/(log10(MaxDistance)-log10(MinDistance));
#endif
  int nbins=60;
  gsl_histogram *corr=gsl_histogram_calloc(nbins);
#ifndef LOGBINS
  gsl_histogram_set_ranges_uniform(corr, MinDistance, MaxDistance);
#else
  gsl_histogram_set_ranges_uniform(corr, LogMinDistance, LogMaxDistance);
#endif
  DDbuf=mymalloc("DDbuf",nbins*sizeof(double));
  DD=mymalloc("DD",nbins*sizeof(double));
  for (i=0; i<nbins; ++i) {
    DDbuf[i]=0.0;
    DD[i]=0.0;
  } //for
  if (ThisTask==0) printf("Calculating the correlation function...\n");
  firsti=ThisTask*Npart/NTasks;
  lasti=(ThisTask+1)*Npart/NTasks;
#ifdef ONEHALOTERM
  if (firsti!=0) {
    while (P[firsti].fofid==P[firsti-1].fofid) firsti++;
  } //if
  if (lasti!=Npart) {
    while (P[lasti-1].fofid==P[lasti].fofid) lasti++;
  } //if
#endif
  slabsize=BoxSize/NTasks;
  for (i=firsti; i<lasti; ++i) {
#ifndef ONEHALOTERM
    for (jj=i+1; (jj<Npart && P[jj].pos[0]-P[i].pos[0]<=ValidDistance) || (jj>=Npart && (BoxSize-P[i].pos[0])+P[jj%Npart].pos[0]<=ValidDistance); ++jj) {
#else
    for (jj=i+1; P[jj].fofid==P[i].fofid && ((jj<Npart && P[jj].pos[0]-P[i].pos[0]<=ValidDistance) || (jj>=Npart && (BoxSize-P[i].pos[0])+P[jj%Npart].pos[0]<=ValidDistance)); ++jj) {
#endif
      j=jj%Npart;
#ifndef SUBREGION
  #ifndef GRIDSHELLS
      dist_sq=SQR(myperiodic(P[i].pos[0]-P[j].pos[0]))+SQR(myperiodic(P[i].pos[1]-P[j].pos[1]))+SQR(myperiodic(P[i].pos[2]-P[j].pos[2]));
  #else
      dist_sq=SQR(myperiodic(gridconv(P[i].pos[0]-P[j].pos[0])))+SQR(myperiodic(gridconv(P[i].pos[1]-P[j].pos[1])))+SQR(myperiodic(gridconv(P[i].pos[2]-P[j].pos[2])));
  #endif
#else
      if (inside_sub(P[i].pos[0],P[i].pos[1],P[i].pos[2]) || inside_sub(P[j].pos[0],P[j].pos[1],P[j].pos[2])) {
  #ifndef GRIDSHELLS
      dist_sq=SQR(P[i].pos[0]-P[j].pos[0])+SQR(P[i].pos[1]-P[j].pos[1])+SQR(P[i].pos[2]-P[j].pos[2]);
  #else
      dist_sq=SQR(gridconv(P[i].pos[0]-P[j].pos[0]))+SQR(gridconv(P[i].pos[1]-P[j].pos[1]))+SQR(gridconv(P[i].pos[2]-P[j].pos[2]));
  #endif
#endif
      if (dist_sq<=SQR(ValidDistance) && dist_sq>=SQR(MinDistance)) {
#ifndef LOGBINS
	ibin=floor((sqrt(dist_sq)-MinDistance)*InvRange*nbins);
        DDbuf[ibin]+=P[i].mass*P[j].mass;
#else
	ibin=floor((0.5*log10(dist_sq)-LogMinDistance)*InvLogRange*nbins);
        DDbuf[ibin]+=P[i].mass*P[j].mass;
#endif
      } //if
#ifdef SUBREGION
      } //if
#endif
    } //for
  } //for
  MPI_Allreduce(DDbuf, DD, nbins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#ifndef SUBREGION
  for (i=0; i<Npart; ++i) MassTot+=P[i].mass;
#else
  for (i=0; i<Npart; ++i) {
    if (inside_sub(P[i].pos[0],P[i].pos[1],P[i].pos[2])) {
      MassTot+=P[i].mass;
    } //if
  } //for
#endif
#ifdef GRIDSHELLS
  int nn=ValidDistance/CellSize;
  shellcounts=mymalloc("shellcounts",nbins*sizeof(int));
  for (i=0; i<nbins; ++i) shellcounts[i]=0;
  for (i=-nn; i<=nn; ++i) {
    for (j=-nn; j<=nn; ++j) {
      for (k=-nn; k<=nn; ++k) {
	dist_sq=SQR(i*CellSize)+SQR(j*CellSize)+SQR(k*CellSize);
	if (dist_sq<=SQR(ValidDistance) && dist_sq>=SQR(MinDistance)) {
  #ifndef LOGBINS
	  ibin=floor((sqrt(dist_sq)-MinDistance)*InvRange*nbins);
  #else
	  ibin=floor((0.5*log10(dist_sq)-LogMinDistance)*InvLogRange*nbins);
  #endif
          shellcounts[ibin]++;
	} //if
      } //for
    } //for
  } //for
#endif
  for (i=0; i<nbins; ++i) {
#ifndef GRIDSHELLS
  #ifndef LOGBINS
    rlow=MinDistance+i*(MaxDistance-MinDistance)/nbins;
    rhigh=MinDistance+(i+1)*(MaxDistance-MinDistance)/nbins;
    RR=2.0*Pi/3.0*SQR(MassTot)/CUB(BoxSize)*(CUB(rhigh)-CUB(rlow));
  #else
    rlow=pow(10,3*(LogMinDistance+i*(LogMaxDistance-LogMinDistance)/nbins));
    rhigh=pow(10,3*(LogMinDistance+(i+1)*(LogMaxDistance-LogMinDistance)/nbins));
    RR=2.0*Pi/3.0*SQR(MassTot)/CUB(BoxSize)*(rhigh-rlow);
  #endif
#else
    RR=0.5*SQR(MassTot)/CUB(BoxSize)*shellcounts[i]*CUB(CellSize);
    //if (ThisTask==0) printf("bin %d:   %d    %f\n",i,shellcounts[i],RR);
#endif
    corr->bin[i]=DD[i]/RR-1.0;
  } //for
  if (ThisTask==0) {
#ifdef LOGBINS
    sprintf(logbit,"_log");
#else
    sprintf(logbit,"");
#endif
#ifdef GRIDSHELLS
    sprintf(gsbit,"_gs");
#else
    sprintf(gsbit,"");
#endif
#ifdef MASSCUT
    float masscutflt=MASSCUT;
    sprintf(massbit,"_mc%.2f",masscutflt);
#elif defined (MASSRANGE)
    sprintf(massbit,"_mr%d_%.2f-%.2f",MassType,MinMassR,MaxMassR);
#else
    sprintf(massbit,"");
#endif
#if defined (RANDOMFOFSAMPLE) || defined (RANDOMGALSAMPLE)
  #ifdef RANDOMFOFSAMPLE
    float randomflt=RANDOMFOFSAMPLE;
    sprintf(rndbit,"_rndfof%.1f",randomflt);
  #else
    float randomflt=RANDOMGALSAMPLE;
    sprintf(rndbit,"_rndgal%.1f",randomflt);
  #endif
#else
    sprintf(rndbit,"");
#endif
#if defined (ONLYCENTRALS) || defined (ONLYSATELLITES)
  #ifdef ONLYCENTRALS
    sprintf(onlybit,"_cen");
  #else
    sprintf(onlybit,"_sat");
  #endif
#else
    sprintf(onlybit,"");
#endif
#ifdef ONEHALOTERM
    sprintf(onebit,"_onehaloterm");
#else
    sprintf(onebit,"");
#endif
#ifdef NO_ORPHANS
    sprintf(orphanbit,"_no_orphans");
#else
    sprintf(orphanbit,"");
#endif
    sprintf(buf,"%s_correl%s%s%s%s%s%s%s.brute",fname,massbit,orphanbit,rndbit,onlybit,onebit,gsbit,logbit);
    fd=my_fopen(buf, "w");
    gsl_histogram_fprintf(fd, corr, "%e", "%e");
    fclose(fd);
    printf(" Writing results to %s\n", buf);
  } //if
#ifdef GRIDSHELLS
  myfree(shellcounts);
#endif
  myfree(DD);
  myfree(DDbuf);
  gsl_histogram_free(corr);
} //brute_correl
