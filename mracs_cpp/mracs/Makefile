#makefile for main, extr, phi, minitest and tinytest
#system type would be Linux or MacOS
SYSTEM	  = Linux

##############
#	Linux	#
##############
ifeq ($(SYSTEM), Linux)
#compiler&flags for Linux
CXX       = g++
OPTIMIZE  = -fopenmp -DMKL_lLP64 -m64 -std=c++20 -O3

OMP_PATH  = -L/usr/local/opt/libomp/lib
FFTW_INCL = -I/usr/local/include

#MRACS used Intel MKL lib for FFT
MKL_ROOT  = /opt/intel/oneapi/mkl/latest
MKL_PATH  = $(MKL_ROOT)/lib
MKL_INCL  = $(MKL_ROOT)/include

#link for Linux
MKL_LINK  = -L$(MKL_PATH) -Wl,--no-as-needed -lmkl_intel_ilp64 \
			-lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl
endif

##############
#	Macos	#
##############
ifeq ($(SYSTEM), MacOS)
#compiler&flags for MacOS
CXX       = clang++
OPTIMIZE  = -fopenmp -std=c++20 -O3 

OMP_PATH  = -L/usr/local/opt/libomp/lib
#LLVM_INCL = -I/usr/local/opt/llvm/include
FFTW_INCL = -I/usr/local/include

#MRACS used Intel MKL lib for FFT
MKL_ROOT  = /opt/intel/oneapi/mkl/latest
MKL_PATH  = $(MKL_ROOT)/lib
MKL_INCL  = $(MKL_ROOT)/include

#link for MacOS
MKL_LINK  = -L$(MKL_PATH) -Wl,-rpath,$(MKL_PATH) -lmkl_cdft_core -lmkl_intel_lp64 \
			-lmkl_tbb_thread -lmkl_core -ltbb -lstdc++ -lpthread -lm
endif


#MRACS source files 
SRC_DIR   = src
BUILD_DIR = build
SRCS      = csmain.cpp readin.cpp kernel.cpp
OBJS     := $(subst .cpp,.o,$(SRCS))
SRCS     := $(addprefix $(SRC_DIR)/, $(SRCS))
OBJS     := $(addprefix $(BUILD_DIR)/, $(OBJS))

#MRACS binary files
BINS     := mracs counting 2pc sigma_r sosta orcorr dipole

.PHONY: all 
.PRECIOUS: $(BUILD_DIR)/%.o
all: mracs counting 2pc sigma_r sosta orcorr dipole

%: $(OBJS) $(BUILD_DIR)/%.o
	$(CXX) $(OPTIMIZE) $(OMP_PATH) $(MKL_LINK) $^ -o $@

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(OPTIMIZE) $(FFTW_INCL) -c $^ -o $@ 


cleanall: clean
	rm -f $(BINS)

clean:
	rm -f $(BUILD_DIR)/*.o
