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
OBJ1  = MO.o
OBJ2  = CT.o
	
all: \
	MO\
	CT

MO: $(OBJ1)
	$(COMP) $(FLAGS) $(LIB) $(OBJ1) -o MO

MO.o MO.mod: MO.f90
	$(COMP) $(FLAGS) $(INCL) -c MO.f90

CT: $(OBJ2)
	$(COMP) $(FLAGS) $(LIB) $(OBJ2) -o CT

CT.o CT.mod: CT.f90
	$(COMP) $(FLAGS) $(INCL) -c CT.f90

clean: 
	/bin/rm *.o
