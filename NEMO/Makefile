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
OBJ1  = nemo.o
OBJ2  = get_nemo.o
	
all: \
	nemo\
	get_nemo

nemo: $(OBJ1)
	$(COMP) $(FLAGS) $(LIB) $(OBJ1) -o nemo

nemo.o nemo.mod: nemo.f90
	$(COMP) $(FLAGS) $(INCL) -c nemo.f90

get_nemo: $(OBJ2)
	$(COMP) $(FLAGS) $(LIB) $(OBJ2) -o get_nemo

get_nemo.o get_nemo.mod: get_nemo.f90
	$(COMP) $(FLAGS) $(INCL) -c get_nemo.f90

clean: 
	/bin/rm *.o
