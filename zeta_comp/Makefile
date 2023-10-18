#
# Makefile for writing nc file code
# CVS:$Id: makefile,v 1.3 2000/04/19 22:05:59 pwjones Exp $
#

COMP = ifort

FLAGS = -O3 -assume byterecl -convert big_endian -fpic -mcmodel large -qopenmp
NETCDF_DIR = $(shell nf-config --prefix)
NETCDF_INC = -I$(NETCDF_DIR)/include
NETCDF_LIB = -L$(NETCDF_DIR)/lib
LIB  = $(NETCDF_LIB) -lnetcdf -lnetcdff
INCL = $(NETCDF_INC)
SRCDIR  =
EXEDIR  = .
OBJ1  = zeta.o
	
all: \
	zeta

zeta: $(OBJ1)
	$(COMP) $(FLAGS) $(LIB) $(OBJ1) -o zeta

zeta.o zeta.mod: zeta.f90
	$(COMP) $(FLAGS) $(INCL) -c zeta.f90

clean: 
	/bin/rm *.o
