#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "allvars.h"
#include "proto.h"
#include "statistics/stats_utils.h"
#include "mracs/particle_io.h"
#include "mracs/mra_grid.h"
#include "mracs/wavelet.h"
#include "statistics/bispectrum_fft.h"
#include "statistics/bispectrum_wavelet.h"
#include "statistics/three_pcf.h"

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
  MPI_Comm_size(MPI_COMM_WORLD, &NTask);
  int i;
  for(PTask = 0; NTask > (1 << PTask); PTask++);

  if(argc != 3)
  {
      if (ThisTask == 0)
      {
          printf("\n  usage: L-VelocityField-Halo-MXXL  <parameterfile>  <outputnum>\n");
          printf("  <parameterfile>    see readparmeterfile.c\n");
          printf("  <outputnum>        snapshot number\n\n");
      }
      exit(1);
    }

  read_parameter_file(argv[1]);

  /* Mode-based dispatch for new pipelines */
  if (strcmp(Mode, "paircount") == 0) {
      if (ThisTask == 0) printf("Mode: MRA pair-counting\n");

      Particle *particles = NULL;
      long long npart;
      MRAGrid mgrid;

      npart = read_millennium_galaxies(SimulationDir, &particles, 1);
      if (npart <= 0) {
          if (ThisTask == 0) printf("Failed to read particles from %s\n", SimulationDir);
          MPI_Finalize();
          return 1;
      }
      if (ThisTask == 0) printf("Read %lld particles\n", npart);

      mra_grid_build(particles, npart, 1ULL << 9, BoxSize, &mgrid);
      if (ThisTask == 0) printf("Built MRA grid: %llu cells\n", mgrid.grid_num);

      /* Convolve with Gaussian kernel */
      WaveletConfig wav_cfg;
      wav_cfg.resolution = 9;
      wav_cfg.base_type = 0;
      wav_cfg.phi = NULL;
      wav_cfg.phi_len = 0;
      mra_grid_convolve(&mgrid, &wav_cfg, KERNEL_GAUSSIAN, Smoothing, 0.0);
      if (ThisTask == 0) printf("Convolved with Gaussian kernel (R=%g)\n", Smoothing);

      mra_grid_write(&mgrid, "mra_density.dat");
      if (ThisTask == 0) printf("Wrote MRA density to mra_density.dat\n");

      free_particles(particles);
      mra_grid_free(&mgrid);

      if (ThisTask == 0) printf("Done (paircount mode)\n");
      MPI_Finalize();
      return 0;
  }

  if (strcmp(Mode, "bispec") == 0) {
      if (ThisTask == 0) printf("Mode: FFT bispectrum\n");
      /* Deferred to integration phase -- loads grid, does FFT, calls bispectrum_fft */
      if (ThisTask == 0) printf("bispec mode not yet integrated\n");
      MPI_Finalize();
      return 0;
  }

  if (strcmp(Mode, "bispec_wavelet") == 0) {
      if (ThisTask == 0) printf("Mode: wavelet bispectrum\n");

      Particle *particles = NULL;
      long long npart;
      npart = read_millennium_galaxies(SimulationDir, &particles, 1);
      if (npart <= 0) {
          if (ThisTask == 0) printf("Failed to read particles\n");
          MPI_Finalize();
          return 1;
      }

      MRAGrid mgrid;
      WaveletConfig wav_cfg;
      wav_cfg.resolution = 8;
      wav_cfg.base_type = 0;
      wav_cfg.phi = NULL;
      wav_cfg.phi_len = 0;

      mra_grid_build(particles, npart, 1ULL << 8, BoxSize, &mgrid);

      int nbins = 20;
      Result *bispec = (Result *) malloc(nbins * sizeof(Result));
      bispectrum_wavelet(&mgrid, &wav_cfg, KERNEL_GAUSSIAN, Smoothing, nbins, bispec);

      if (ThisTask == 0) {
          char fname[256];
          snprintf(fname, sizeof(fname), "%s_bispec_wavelet.dat", OutputDir);
          FILE *fp = fopen(fname, "w");
          if (fp) {
              int i;
              for (i = 0; i < nbins; i++) {
                  fprintf(fp, "%g %g %lld\n",
                          bispec[i].x, bispec[i].data, bispec[i].n);
              }
              fclose(fp);
              printf("Wrote wavelet bispectrum to %s\n", fname);
          }
      }

      free(bispec);
      free_particles(particles);
      mra_grid_free(&mgrid);
      if (ThisTask == 0) printf("Done (bispec_wavelet mode)\n");
      MPI_Finalize();
      return 0;
  }

  if (strcmp(Mode, "three_pcf") == 0) {
      if (ThisTask == 0) printf("Mode: 3-point correlation function\n");

      Particle *particles = NULL;
      long long npart;
      npart = read_millennium_galaxies(SimulationDir, &particles, 1);
      if (npart <= 0) {
          if (ThisTask == 0) printf("Failed to read particles\n");
          MPI_Finalize();
          return 1;
      }
      if (ThisTask == 0) printf("Read %lld particles\n", npart);

      WaveletConfig wav_cfg;
      wav_cfg.resolution = 8;
      wav_cfg.base_type = 0;
      wav_cfg.phi = NULL;
      wav_cfg.phi_len = 0;

      int nbins_r = 15;
      Result *tcf = (Result *) malloc(nbins_r * sizeof(Result));
      three_pcf(particles, npart, &wav_cfg, KERNEL_GAUSSIAN,
                Smoothing, BoxSize, nbins_r, tcf);

      if (ThisTask == 0) {
          char fname[256];
          snprintf(fname, sizeof(fname), "%s_three_pcf.dat", OutputDir);
          FILE *fp = fopen(fname, "w");
          if (fp) {
              int i;
              for (i = 0; i < nbins_r; i++) {
                  fprintf(fp, "%g %g %lld\n",
                          tcf[i].x, tcf[i].data, tcf[i].n);
              }
              fclose(fp);
              printf("Wrote 3PCF to %s\n", fname);
          }
      }

      free(tcf);
      free_particles(particles);
      if (ThisTask == 0) printf("Done (three_pcf mode)\n");
      MPI_Finalize();
      return 0;
  }

  SnapshotNum = atoi(argv[2]);

  HighMark_run = 0;

  mymalloc_init();

  init_grid();
 
#ifdef GEN_HALOFIELD
#ifndef DMP
  /* --- load in cata data --- */
  TotNumHalo = 0;
  load_catalogue();

  for (Dim = 0; Dim < 4; Dim ++)
  {
	memset(Grid, 0, maxfftsize * sizeof(fftw_real));
	density(PP, &Grid, LocNgroups, 1, Dim);
	save_grids(Dim);
  }

  myfree(PP);
#else
  load_snapshot();
#endif
#else
  /* --- load in mesh data --- */
  get_slabs_file();
#if defined (SMOOTH_HALOFIELD) && defined (SMOOTH_HALOVELFIELD)
  for (Dim = 0; Dim < 3; Dim ++)
  {
    if (ThisTask == 0)
    {
      if(Dim == 0) printf("Smoothing Vx Field \n");
      if(Dim == 1) printf("Smoothing Vy Field \n");
      if(Dim == 2) printf("Smoothing Vz Field \n");
    }
#endif
  memset(Grid, 0, maxfftsize * sizeof(fftw_real));
  load_grid_field(Dim);

  /* --- fft dens data --- */
  fft_Grid = (fftw_complex *) &Grid[0];
  if(ThisTask == 0)
      printf("fft forward \n");
  fftw_execute(fft_forward_plan);

  for (rep=0; rep<=0; rep++)
  {  
      if ( rep == 0 ){ 
          Smoothing = 0.0000001;
      } else {
          Smoothing = 1.25* pow(2, rep-1);
      }

      if (ThisTask == 0)
          printf("Smoothing :%g\n",Smoothing);
#ifdef SMOOTH_HALOFIELD
	  smooth();
	  save_grids(Dim);
	  myfree(sGrid);
#endif
#ifdef LINEAR_VELFIELD
	  for (Dim=0; Dim<3; Dim++)
	  {
        velocity_field(Dim);
		save_grids(Dim);
		myfree(Vel);
	  }
#endif

#if defined (SMOOTH_HALOFIELD) && defined (SMOOTH_HALOVELFIELD)
  }
#endif
  }
#endif

  clean_grid();

 
  MPI_Barrier(MPI_COMM_WORLD);
  if(ThisTask == 0)
      printf("Done\n"); 

  MPI_Finalize();		/* clean up & finalize MPI */

  return 0;
}




