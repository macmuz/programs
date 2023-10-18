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
OBJ1  = new_grid.o
OBJ2  = HRDM_bathy.o
	
all: \
	new_grid\
	HRDM_bathy

new_grid: $(OBJ1)
	$(COMP) $(FLAGS) $(LIB) $(OBJ1) -o new_grid

new_grid.o new_grid.mod: new_grid.f90
	$(COMP) $(FLAGS) $(INCL) -c new_grid.f90

HRDM_bathy: $(OBJ2)
	$(COMP) $(FLAGS) $(LIB) $(OBJ2) -o HRDM_bathy

HRDM_bathy.o HRDM_bathy.mod: HRDM_bathy.f90
	$(COMP) $(FLAGS) $(INCL) -c HRDM_bathy.f90

clean: 
	/bin/rm *.o
