
EXEC   = L-VelHalo-MXXL

OBJS   = main.o mymalloc.o read_parameters.o allvars.o io.o smooth.o velocity_field.o mxxl_io.o 

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

#SYSTYPE="Ubuntu"
#SYSTYPE="OpenSuse"
#SYSTYPE="OpenSuse64"
#SYSTYPE="Genius"
#SYSTYPE="HLRB2"
#SYSTYPE="MPA"
#SYSTYPE="VIP"

#SYSTYPE="OpteronMPA-Gnu"
#SYSTYPE="OpteronMPA-Intel"

# module load mvapich2-1.2-sdr-intel/11.0
#SYSTYPE="OPA-Cluster64-Intel"

# module load mvapich2-1.2-sdr-gnu/4.1.2
#SYSTYPE="OPA-Cluster64-Gnu"

#SYSTYPE="JSC"
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



ifeq (SOFTDOUBLEDOUBLE,$(findstring SOFTDOUBLEDOUBLE,$(OPT)))
CC       =   mpiCC     # default C++ compiler
OPTIMIZE =   -g 
OPT     +=  -DX86FIX   # only needed for 32-bit intel/amd systems
endif

ifeq (VORONOI,$(findstring VORONOI,$(OPT)))
GMPLIB   =  -lgmp
endif

ifeq ($(SYSTYPE),"JSC")
CC       =  mpicc   # sets the C-compiler
#CC       = /opt/parastation/mpi2-intel-mt/bin/mpicc
OPTIMIZE =  -g -O3 -Wall -sox
ifeq (OPENMP,$(findstring OPENMP,$(OPT)))
OPTIMIZE +=  -openmp
endif
GSL_INCL = -I/lustre/jwork1/vsk021/libs/include
GSL_LIBS = -L/lustre/jwork1/vsk021/libs/lib -Xlinker -R -Xlinker /lustre/jwork1/vsk021/libs/lib
FFTW_INCL= -I/lustre/jwork1/vsk021/libs/include
FFTW_LIBS= -L/lustre/jwork1/vsk021/libs/lib -Xlinker -R -Xlinker /lustre/jwork1/vsk021/libs/lib
MPICHLIB = -L/opt/parastation/mpi2-intel/lib/ -Xlinker -R -Xlinker /opt/parastation/mpi2-intel/lib/ -lmpi
MPICHINC = -I/opt/parastation/mpi2-intel/include/
HDF5LIB  =
endif


ifeq ($(SYSTYPE),"Genius")
CC       =  mpixlc_r
OPTIMIZE =  -O3 -qarch=450 -qtune=450
GSL_INCL =  -I/u/vrs/BlueGene/Libs/include
GSL_LIBS =  -L/u/vrs/BlueGene/Libs/lib -static
FFTW_INCL=  -I/u/vrs/BlueGene/Libs/include
FFTW_LIBS=  -L/u/vrs/BlueGene/Libs/lib
MPICHLIB =
HDF5INCL =
HDF5LIB  =
OPT      +=  -DNOCALLSOFSYSTEM
endif



ifeq ($(SYSTYPE),"HLRB2")
CC       =  icc   # sets the C-compiler
OPTIMIZE =  -g -O3 -ipo  -mtune=itanium2 -mcpu=itanium2 -mp
GSL_INCL = -I/home/hlrb2/h1091/h1091aa/Libs/include
GSL_LIBS = -L/home/hlrb2/h1091/h1091aa/Libs/lib -Xlinker -R -Xlinker /home/hlrb2/h1091/h1091aa/Libs/lib
FFTW_INCL= -I/home/hlrb2/h1091/h1091aa/Libs/include
FFTW_LIBS= -L/home/hlrb2/h1091/h1091aa/Libs/lib -Xlinker -R -Xlinker /home/hlrb2/h1091/h1091aa/Libs/lib
MPICHLIB = -lmpi
HDF5INCL = 
HDF5LIB  = 
endif

ifeq ($(SYSTYPE),"Ubuntu")
CC       =  mpicc -Wall
OPTIMIZE =  
GSL_INCL =  
GSL_LIBS =  
FFTW_INCL=  -I/usr/include
FFTW_LIBS= -L/usr/lib
MPICHLIB = -L/usr/lib/mpi/gcc/openmpi/lib -Xlinker -R -Xlinker /usr/lib/mpi/gcc/openmpi/lib -lmpi 
HDF5INCL =
HDF5LIB  =
endif


ifeq ($(SYSTYPE),"OpenSuse")
CC       =  mpicc -Wall
OPTIMIZE =  
GSL_INCL =  
GSL_LIBS =  
FFTW_INCL=  
FFTW_LIBS=  
MPICHLIB = -L/usr/lib/mpi/gcc/openmpi/lib -Xlinker -R -Xlinker /usr/lib/mpi/gcc/openmpi/lib -lmpi 
HDF5INCL =
HDF5LIB  =
endif


ifeq ($(SYSTYPE),"OpenSuse64")
CC       =  mpicc
OPTIMIZE =  
GSL_INCL =  
GSL_LIBS =  
FFTW_INCL=  
FFTW_LIBS=  
MPICHLIB = -L/usr/lib64/mpi/gcc/openmpi/lib64 -Xlinker -R -Xlinker /usr/lib64/mpi/gcc/openmpi/lib64 -lmpi 
HDF5INCL =
HDF5LIB  =
endif




ifeq ($(SYSTYPE),"OPA-Cluster64-Intel")  # module load mvapich2-1.2-sdr-intel/11.0
CC       =   mpicc
ifeq (SOFTDOUBLEDOUBLE,$(findstring SOFTDOUBLEDOUBLE,$(OPT))) 
CC       =   mpicxx
endif
OPTIMIZE =  -O3 -g -m64 #-Wall 
ifeq (OPENMP,$(findstring OPENMP,$(OPT))) 
OPTIMIZE +=  -openmp
endif
GSL_INCL =  -I/afs/rzg/home/v/vrs/Libs/opteron64/mvapich2-1.2-sdr-intel-11.0/include
GSL_LIBS =  -L/afs/rzg/home/v/vrs/Libs/opteron64/mvapich2-1.2-sdr-intel-11.0/lib  -Xlinker -R -Xlinker /afs/rzg/home/v/vrs/Libs/opteron64/mvapich2-1.2-sdr-intel-11.0/lib
FFTW_INCL=  -I/afs/rzg/home/v/vrs/Libs/opteron64/mvapich2-1.2-sdr-intel-11.0/include
FFTW_LIBS=  -L/afs/rzg/home/v/vrs/Libs/opteron64/mvapich2-1.2-sdr-intel-11.0/lib
MPICHLIB =
HDF5INCL =  -I/afs/rzg/home/v/vrs/Libs/opteron64/mvapich2-1.2-sdr-intel-11.0/include
HDF5LIB  =  -L/afs/rzg/home/v/vrs/Libs/opteron64/mvapich2-1.2-sdr-intel-11.0/lib  -lhdf5 -lz 
ifeq (FFTW3,$(findstring FFTW3,$(OPT))) 
FFTW_INCL=  -I/afs/ipp-garching.mpg.de/home/r/reangulo/libs/mvapich2-1.2-sdr-gnu-11.0/include
FFTW_LIBS=  -L/afs/ipp-garching.mpg.de/home/r/reangulo/libs/mvapich2-1.2-sdr-gnu-11.0/lib
endif
OPT      +=  -DNOCALLSOFSYSTEM
endif


ifeq ($(SYSTYPE),"OPA-Cluster64-Gnu")  # module load mvapich2-1.2-sdr-gnu/4.1.2
CC       =   mpicc
ifeq (SOFTDOUBLEDOUBLE,$(findstring SOFTDOUBLEDOUBLE,$(OPT))) 
CC       =   mpicxx
endif
OPTIMIZE =  -O3 -g -Wall -m64
ifeq (OPENMP,$(findstring OPENMP,$(OPT))) 
OPTIMIZE +=  -fopenmp
endif
GSL_INCL =  -I/afs/rzg/home/v/vrs/Libs/opteron64/mvapich2-1.2-sdr-gnu-4.1.2/include
GSL_LIBS =  -L/afs/rzg/home/v/vrs/Libs/opteron64/mvapich2-1.2-sdr-gnu-4.1.2/lib  -Xlinker -R -Xlinker /afs/rzg/home/v/vrs/Libs/opteron64/mvapich2-1.2-sdr-gnu-4.1.2/lib
MPICHLIB =
HDF5INCL =  -I/afs/rzg/home/v/vrs/Libs/opteron64/mvapich2-1.2-sdr-gnu-4.1.2/include
HDF5LIB  =  -L/afs/rzg/home/v/vrs/Libs/opteron64/mvapich2-1.2-sdr-gnu-4.1.2/lib  -lhdf5 -lz 
OPT      +=  -DNOCALLSOFSYSTEM
FFTW_INCL=  -I/afs/rzg/home/v/vrs/Libs/opteron64/mvapich2-1.2-sdr-gnu-4.1.2/include
FFTW_LIBS=  -L/afs/rzg/home/v/vrs/Libs/opteron64/mvapich2-1.2-sdr-gnu-4.1.2/lib
ifeq (FFTW3,$(findstring FFTW3,$(OPT))) 
FFTW_INCL=  -I/afs/ipp-garching.mpg.de/home/r/reangulo/libs/mvapich2-1.2-sdr-gnu-11.0/include
FFTW_LIBS=  -L/afs/ipp-garching.mpg.de/home/r/reangulo/libs/mvapich2-1.2-sdr-gnu-11.0/lib
endif
endif


ifeq ($(SYSTYPE),"OpteronMPA-Intel")
CC       =  /usr/local/appl/openmpi-intel/bin/mpicc  -m64 # sets the C-compiler
ifeq (SOFTDOUBLEDOUBLE,$(findstring SOFTDOUBLEDOUBLE,$(OPT)))
CC       =  /usr/local/appl/openmpi-intel/bin/mpiCC  -m64
endif
OPTIMIZE =   -O3 -Wall
FFTW_INCL=  -I/afs/mpa/home/volker/Libs/include
FFTW_LIBS=  -L/afs/mpa/home/volker/Libs/lib -Xlinker -R -Xlinker /afs/mpa/home/volker/Libs/lib
MPICHLIB =
HDF5INCL =  -I/afs/mpa/home/volker/Libs/include
HDF5LIB  =  -L/afs/mpa/home/volker/Libs/lib -Xlinker -R -Xlinker /afs/mpa/home/volker/Libs/lib -lhdf5 -lz
endif

ifeq ($(SYSTYPE),"OpteronMPA-Gnu")
CC       =  /usr/local/appl/openmpi-gnu/bin/mpicc  -m64 # sets the C-compiler
ifeq (SOFTDOUBLEDOUBLE,$(findstring SOFTDOUBLEDOUBLE,$(OPT)))
CC       =  /usr/local/appl/openmpi-gnu/bin/mpiCC  -m64
endif
OPTIMIZE =   -O3 -Wall
GSL_INCL =  
GSL_LIBS =  
FFTW_INCL=  -I/afs/mpa/home/volker/Libs/include
FFTW_LIBS=  -L/afs/mpa/home/volker/Libs/lib -Xlinker -R -Xlinker /afs/mpa/home/volker/Libs/lib
MPICHLIB =
HDF5INCL =  -I/afs/mpa/home/volker/Libs/include
HDF5LIB  =  -L/afs/mpa/home/volker/Libs/lib -Xlinker -R -Xlinker /afs/mpa/home/volker/Libs/lib -lhdf5 -lz
endif


ifeq ($(SYSTYPE),"MPA")
CC       =  mpicc   # sets the C-compiler
ifeq (SOFTDOUBLEDOUBLE,$(findstring SOFTDOUBLEDOUBLE,$(OPT)))
CC       =  mpiCC
OPT     +=  -DX86FIX
endif
OPTIMIZE =   -g -Wall 
# GSL_INCL =  -I/usr/common/pdsoft/include
# GSL_LIBS =  -L/usr/common/pdsoft/lib
GSL_INCL =  -I/afs/mpa/home/volker/Libs/include
GSL_LIBS =  -L/afs/mpa/home/volker/Libs/lib
FFTW_INCL=  -I/afs/mpa/home/volker/Libs/include
FFTW_LIBS=  -L/afs/mpa/home/volker/Libs/lib -Xlinker -R -Xlinker /afs/mpa/home/volker/Libs/lib
MPICHLIB =
HDF5INCL =  -I/afs/mpa/home/volker/Libs/include
HDF5LIB  =  -L/afs/mpa/home/volker/Libs/lib -lhdf5 -lz 
endif


ifeq ($(SYSTYPE),"VIP")
CC       =  mpcc_r -g -qcpluscmt # -qflttrap=enable:zerodivide:nanq # sets the C-compiler
ifeq (FLTROUNDOFFREDUCTION,$(findstring FLTROUNDOFFREDUCTION,$(OPT)))
CC       =  mpcc_r  -qldbl128 -lC128_r  -qcpluscmt # sets the C-compiler
#            (this compiler has native support for 128bit long double, SOFTDOUBLEDOUBLE n\
ot needed)
endif
OPTIMIZE =  -q64 #-O3 -qarch=auto -qtune=auto -qinline # -qhot -qipa=level=2 -qfloat=hdflt
GSL_INCL = -I/afs/rzg/u/kud/include
GSL_LIBS = -L/afs/rzg/u/kud/lib -bstackpsize:64k -bdatapsize:64k -btextpsize:64k
FFTW_INCL= -I/afs/rzg/u/kud/include
FFTW_LIBS= -L/afs/rzg/u/kud/lib  -q64 # -qipa
MPICHLIB =
HDF5INCL =
HDF5LIB  =
OPT     +=  # -DPOWER6 -DALT_QSORT -DPAD_STRUCTURES
endif

ifeq (NOTYPEPREFIX_FFTW,$(findstring NOTYPEPREFIX_FFTW,$(OPT)))    # fftw installed with type prefix?
  FFTW_LIB = $(FFTW_LIBS) -lrfftw_mpi -lfftw_mpi -lrfftw -lfftw 
else
ifeq (DOUBLEPRECISION_FFTW,$(findstring DOUBLEPRECISION_FFTW,$(OPT)))
  FFTW_LIB = $(FFTW_LIBS) -ldrfftw_mpi -ldfftw_mpi -ldrfftw -ldfftw 
else
  FFTW_LIB = $(FFTW_LIBS) -lsrfftw_mpi -lsfftw_mpi -lsrfftw -lsfftw 
endif
endif


LIBS   =  $(GSL_LIBS) -lgsl -lgslcblas $(FFTW_LIB) -lm -g

CFLAGS = $(OPTIMIZE) $(OPTIONS) $(GSL_INCL) $(FFTW_INCL)

$(EXEC): $(OBJS) 
	$(CC) $(OPTIONS) $(OPTIMIZE) $(OBJS) $(LIBS)   -o  $(EXEC)  

$(OBJS): $(INCL) 


clean:
	rm -f $(OBJS) $(EXEC)









