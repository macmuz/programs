#
# Makefile for writing nc file code
# CVS:$Id: makefile,v 1.3 2000/04/19 22:05:59 pwjones Exp $
#

COMP = mpiifort

FLAGS = -O3 -assume byterecl -convert big_endian -fpic -mcmodel large -qopenmp
NETCDF_DIR = $(shell nf-config --prefix)
NETCDF_INC = -I$(NETCDF_DIR)/include
NETCDF_LIB = -L$(NETCDF_DIR)/lib
LIB  = $(NETCDF_LIB) -lnetcdf -lnetcdff
INCL = $(NETCDF_INC)
SRCDIR  =
EXEDIR  = .
OBJ1  = energy_ROMS.o
OBJ2  = energy_NEMO.o
	
all: \
	energy_ROMS\
	energy_NEMO

energy_ROMS: $(OBJ1)
	$(COMP) $(FLAGS) $(LIB) $(OBJ1) -o energy_ROMS

energy_ROMS.o energy_ROMS.mod: energy_ROMS.f90
	$(COMP) $(FLAGS) $(INCL) -c energy_ROMS.f90

energy_NEMO: $(OBJ2)
	$(COMP) $(FLAGS) $(LIB) $(OBJ2) -o energy_NEMO

energy_NEMO.o energy_NEMO.mod: energy_NEMO.f90
	$(COMP) $(FLAGS) $(INCL) -c energy_NEMO.f90

clean: 
	/bin/rm *.o
