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
OBJ1  = new_flat.o 
OBJ2  = flat.o
	
all: \
	flat\
	new_flat

flat: $(OBJ2)
	$(COMP) $(FLAGS) $(LIB) $(OBJ2) -o flat

flat.o flat.mod: flat.f90
	$(COMP) $(FLAGS) $(INCL) -c flat.f90

new_flat: $(OBJ1)
	$(COMP) $(FLAGS) $(LIB) $(OBJ1) -o new_flat

new_flat.o new_flat.mod: new_flat.f90
	$(COMP) $(FLAGS) $(INCL) -c new_flat.f90

clean: 
	/bin/rm *.o
