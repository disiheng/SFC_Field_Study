EXEC   = L-VelHalo-MXXL

OBJS   = main.o mymalloc.o read_parameters.o allvars.o io.o smooth.o \
         velocity_field.o mxxl_io.o field_assign.o \
         statistics/stats_utils.o statistics/bispectrum_fft.o \
         statistics/bispectrum_wavelet.o statistics/three_pcf.o \
         mracs/particle_io.o mracs/sfc.o mracs/wavelet.o \
         mracs/kernel.o mracs/mra_grid.o

INCL   = allvars.h  proto.h  Makefile

#OPTIONS  += -DGEN_HALOFIELD
#OPTIONS  += -DCIC
#OPTIONS  += -DDMP
#OPTIONS  += -DSUBHALO

OPTIONS  += -DSMOOTH_HALOFIELD
#OPTIONS  += -DSMOOTH_HALOVELFIELD
#OPTIONS  += -DLINEAR_VELFIELD

#OPTIONS  += -DLOGDENS

#OPTIONS += -DDOUBLEPRECISION_FFTW
#OPTIONS += -DSWAP_ENDIAN
#OPTIONS += -DALPHA=0.8                # Exponent in weighting function

SYSTYPE="COSM"

CC       = mpicc        # sets the C-compiler (default)
OPTIMIZE = -Wall  -g   # optimization and warning flags (default)

MPICHLIB = -lmpich

ifeq ($(SYSTYPE),"COSM")
CC       =  mpicc   # sets the C-compiler
#CC       = /opt/parastation/mpi2-intel-mt/bin/mpicc
OPTIMIZE =  -g -O3 #-O3 -Wall -m64 -wd981 -wd2259 -wd1572 -ipoG
GSL_INCL = -I$(GSL_HOME)/include
GSL_LIBS = -L$(GSL_HOME)/lib -Xlinker -R -Xlinker $(GSL_HOME)/lib
FFTW_INCL= -I$(FFTW_HOME)/include
FFTW_LIBS= -L$(FFTW_HOME)/lib -Xlinker -R -Xlinker $(FFTW_HOME)/lib
MPICHLIB =
MPICHINC =
HDF5LIB  =
endif

# FFTW3: link fftw3_mpi + fftw3 (double precision)
FFTW_LIB = $(FFTW_LIBS) -lfftw3_mpi -lfftw3

LIBS   =  $(GSL_LIBS) -lgsl -lgslcblas $(FFTW_LIB) -lm -g

CFLAGS = $(OPTIMIZE) $(OPTIONS) $(GSL_INCL) $(FFTW_INCL) -I.

$(EXEC): $(OBJS)
	$(CC) $(OPTIONS) $(OPTIMIZE) $(OBJS) $(LIBS)   -o  $(EXEC)

$(OBJS): $(INCL)

statistics/%.o: statistics/%.c
	$(CC) $(CFLAGS) -c $< -o $@

mracs/%.o: mracs/%.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(EXEC)
