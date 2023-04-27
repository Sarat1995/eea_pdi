# Makefile

# MACHINE = debug
MACHINE = mac
# MACHINE = solidmatters
# MACHINE = ubuntu
# MACHINE = millipede-intel
# MACHINE = millipede


ifeq ($(MACHINE), debug)
CC = gcc
CC_FLAGS = -O0 -g -w
LAPACK = -framework veclib
FFTW = -lfftw3 -lm -L/opt/local/lib
FFTWI = -I/opt/local/include
CPPFLAGS = $(FFTWI)
endif

ifeq ($(MACHINE), mac)
CC = gcc
CC_FLAGS = -O3
LAPACK = -llapack
FFTW = -lfftw3 -lm -L/opt/local/lib
FFTWI = -I/opt/local/include
CPPFLAGS = $(FFTWI)
endif

ifeq ($(MACHINE), obsolete)
CC = gcc
CC_FLAGS = -O3
LAPACK = -framework Accelerate
FFTW = -lfftw3 -lm -L/opt/local/lib
FFTWI = -I/opt/local/include
CPPFLAGS = $(FFTWI)
endif

ifeq ($(MACHINE), solidmatters)
CC = icc
FFTWI = -I/software/apps/fftw/intel/3.3.3/include
FFTW = -lfftw3 -L/software/apps/fftw/intel/3.3.3/lib
CC_FLAGS = -B -O3 -lm
LAPACK = -llapack
LAPACKL = -L/usr/lib64/atlas
CPPFLAGS = $(FFTWI)
endif

ifeq ($(MACHINE), ubuntu)
CC = gcc-9
CC_FLAGS = -O3 -lm
LAPACK = -llapack
FFTW = -lfftw3
endif

ifeq ($(MACHINE), millipede-intel)
CC = icc
FFTWI = -I/cm/shared/apps/fftw/gcc/64/3.2.2/include
FFTW = -lfftw3 -L/cm/shared/apps/fftw/gcc/64/3.2.2/lib
CC_FLAGS = -B -O3 -lm
MKLDIR = /cm/shared/apps/intel-cluster-runtime/2.1/mkl/10.2.0.013/lib/em64t
LAPACK = -L$(MKLDIR) -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -Wl,-rpath=$(MKLDIR)
CPPFLAGS = $(FFTWI)
endif

ifeq ($(MACHINE), millipede)
CC = gcc
FFTWI = -I/cm/shared/apps/fftw/gcc/64/3.2.2/include
FFTW = -lfftw3 -L/cm/shared/apps/fftw/gcc/64/3.2.2/lib
CC_FLAGS = -B -O3 -lm
LAPACK = -llapack -lblas -lgcc
LAPACKL = -L/cm/shared/apps/lapack/gcc/64/3.2.1 -L/cm/shared/apps/blas/gcc/1/lib64
CPPFLAGS = $(FFTWI)
endif


all: ../bin/DirectRun ../bin/VibroNISE ../bin/Spectra ../bin/Pulse


../bin/DirectRun:			DirectRun.o ParmsHandleMod.o ToolsMod.o
	$(CC) -o ../bin/DirectRun	DirectRun.o ParmsHandleMod.o ToolsMod.o $(LAPACK) $(CC_FLAGS)

../bin/VibroNISE:			VibroNISE.o SubsMod.o ParmsHandleMod.o OverlapsMod.o BasisMod.o ToolsMod.o RandomMod.o
	$(CC) -o ../bin/VibroNISE	VibroNISE.o SubsMod.o ParmsHandleMod.o OverlapsMod.o BasisMod.o ToolsMod.o RandomMod.o $(LAPACK) $(CC_FLAGS)

../bin/Pulse:				Pulse.o ParmsHandleMod.o ToolsMod.o
	$(CC) -o ../bin/Pulse		Pulse.o ParmsHandleMod.o ToolsMod.o $(LAPACK) $(CC_FLAGS)

../bin/Spectra:			Spectra.o ParmsHandleMod.o ToolsMod.o
	$(CC) -o ../bin/Spectra	Spectra.o ParmsHandleMod.o ToolsMod.o $(LAPACK) $(FFTW) $(CC_FLAGS)

DirectRun.o:			DirectRun.c ParmsHandleMod.h ToolsMod.h ParmsMod.h GlobalsMod.h

VibroNISE.o:			VibroNISE.c VibroNISE.h SubsMod.h ParmsHandleMod.h OverlapsMod.h BasisMod.h ToolsMod.h RandomMod.h ParmsMod.h GlobalsMod.h

Pulse.o:			Pulse.c ParmsHandleMod.h ToolsMod.h ParmsMod.h GlobalsMod.h

Spectra.o:			Spectra.c ParmsHandleMod.h ToolsMod.h ParmsMod.h GlobalsMod.h

SubsMod.o:			SubsMod.c SubsMod.h BasisMod.h ToolsMod.h RandomMod.h ParmsMod.h GlobalsMod.h

ParmsHandleMod.o:		ParmsHandleMod.c ParmsHandleMod.h ToolsMod.h ParmsMod.h GlobalsMod.h

OverlapsMod.o:			OverlapsMod.c OverlapsMod.h GlobalsMod.h

BasisMod.o:			BasisMod.c BasisMod.h ToolsMod.h GlobalsMod.h

ToolsMod.o:			ToolsMod.c ToolsMod.h GlobalsMod.h

RandomMod.o:			RandomMod.c RandomMod.h

clean:
	rm -f *.o ../bin/VibroNISE ../bin/Spectra ../bin/Pulse ../bin/DirectRun
