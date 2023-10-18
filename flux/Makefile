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
OBJ1  = flux.o
OBJ2  = flux_outflow.o
	
all: \
	flux\
	flux_outflow

flux: $(OBJ1)
	$(COMP) $(FLAGS) $(LIB) $(OBJ1) -o flux

flux.o flux.mod: flux.f90
	$(COMP) $(FLAGS) $(INCL) -c flux.f90

flux_outflow: $(OBJ2)
	$(COMP) $(FLAGS) $(LIB) $(OBJ2) -o flux_outflow

flux_outflow.o flux_outflow.mod: flux_outflow.f90
	$(COMP) $(FLAGS) $(INCL) -c flux_outflow.f90

clean: 
	/bin/rm *.o
