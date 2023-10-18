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
OBJ  = hfilterbc.o
	
all: \
	hfilterbc

hfilterbc: $(OBJ)
	$(COMP) $(FLAGS) $(LIB) $(OBJ) -o hfilterbc

hfilterbc.o hfilterbc.mod: hfilterbc.f90
	$(COMP) $(FLAGS) $(INCL) -c hfilterbc.f90

clean: 
	/bin/rm *.o
