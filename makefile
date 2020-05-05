#-------------------------------------------------------------------------------

HOST    := $(shell hostname)
OS      := $(shell uname)

HDF5DIR = /usr
FITSDIR = /usr

FC      = gfortran
FFLAGS  = -O3 -march=native -fopenmp -I$(HDF5DIR)/include
LD      = $(FC)
LDFLAGS = $(FFLAGS)
LIBS    =  -L$(HDF5DIR)/lib -lhdf5_fortran -lhdf5 -L$(FITSDIR)/lib -lcfitsio -lz

#-------------------------------------------------------------------------------

.SUFFIXES:
.SUFFIXES: .f90 .o

.f90.o:
	$(FC) -c $(FFLAGS) $*.f90

#-------------------------------------------------------------------------------

name = wtmodes

default: $(name).x

sources = $(name).f90 params.f90 fits.f90 hdf5.f90 wavelet.f90
objects = $(name).o   params.o   fits.o   hdf5.o   wavelet.o
files   = $(sources) makefile params.in

$(name).x: $(objects)
	$(LD) $(LDFLAGS) $(objects) -o $(name).x $(LIBS)

dist: $(files)
	tar czvf $(name).tar.gz $(files)

clean:
	rm -rf *.x *.o *.mod *.il

clean-all:
	rm -rf *.x *.o *.mod *.il *.dat *.bin *~

#-------------------------------------------------------------------------------

wtmodes.o  : wtmodes.f90 params.o fits.o hdf5.o wavelet.o
params.o   : params.f90
fits.o     : fits.f90 params.o
hdf5.o     : hdf5.f90 params.o
wavelet.o  : wavelet.f90 params.o fits.o hdf5.o
random.o   : random.f90

#-------------------------------------------------------------------------------
