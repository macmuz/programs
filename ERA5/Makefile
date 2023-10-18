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
OBJ1  = forcing.o
OBJ2  = forcing2023.o
	
all: \
	forcing\
	forcing2023

forcing: $(OBJ1)
	$(COMP) $(FLAGS) $(LIB) $(OBJ1) -o forcing

forcing.o forcing.mod: forcing.f90
	$(COMP) $(FLAGS) $(INCL) -c forcing.f90

forcing2023: $(OBJ2)
	$(COMP) $(FLAGS) $(LIB) $(OBJ2) -o forcing2023

forcing2023.o forcing2023.mod: forcing2023.f90
	$(COMP) $(FLAGS) $(INCL) -c forcing2023.f90

clean: 
	/bin/rm *.o
