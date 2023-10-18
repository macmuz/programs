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
OBJ2  = prepare_bathy.o
OBJ3  = move_forcing.o
OBJ4  = meo_scheme.o
	
all: \
	prepare_bathy\
	rx1check\
	move_forcing\
	meo_scheme

rx1check: $(OBJ1)
	$(COMP) $(FLAGS) $(LIB) $(OBJ1) -o rx1check

rx1check.o rx1check.mod: rx1check.f90
	$(COMP) $(FLAGS) $(INCL) -c rx1check.f90

prepare_bathy: $(OBJ2)
	$(COMP) $(FLAGS) $(LIB) $(OBJ2) -o prepare_bathy

prepare_bathy.o prepare_bathy.mod: prepare_bathy.f90
	$(COMP) $(FLAGS) $(INCL) -c prepare_bathy.f90

move_forcing: $(OBJ3)
	$(COMP) $(FLAGS) $(LIB) $(OBJ3) -o move_forcing

move_forcing.o move_forcing.mod: move_forcing.f90
	$(COMP) $(FLAGS) $(INCL) -c move_forcing.f90

meo_scheme: $(OBJ4)
	$(COMP) $(FLAGS) $(LIB) $(OBJ4) -o meo_scheme

meo_scheme.o meo_scheme.mod: meo_scheme.f90
	$(COMP) $(FLAGS) $(INCL) -c meo_scheme.f90

clean: 
	/bin/rm *.o
