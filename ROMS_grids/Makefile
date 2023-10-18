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
OBJ1  = rx1check.o
OBJ2  = ROMS_grid.o
	
all: \
	ROMS_grid\
	rx1check

rx1check: $(OBJ1)
	$(COMP) $(FLAGS) $(LIB) $(OBJ1) -o rx1check

rx1check.o rx1check.mod: rx1check.f90
	$(COMP) $(FLAGS) $(INCL) -c rx1check.f90

ROMS_grid: $(OBJ2)
	$(COMP) $(FLAGS) $(LIB) $(OBJ2) -o ROMS_grid

ROMS_grid.o ROMS_grid.mod: ROMS_grid.f90
	$(COMP) $(FLAGS) $(INCL) -c ROMS_grid.f90

clean: 
	/bin/rm *.o
