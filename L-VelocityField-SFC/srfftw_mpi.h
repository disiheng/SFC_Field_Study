#ifndef SRFFTW_MPI_H
#define SRFFTW_MPI_H

/* Stub for single-precision FFTW MPI header.
 * Replace with real FFTW installation for full compilation. */
#define SEND_TYPE MPI_FLOAT

typedef float fftw_real;
typedef struct { fftw_real re, im; } fftw_complex;
typedef struct rfftwnd_mpi_plan_ { int dummy; } *rfftwnd_mpi_plan;

#endif
