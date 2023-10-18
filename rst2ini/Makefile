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
OBJ1  = rst2ini.o
	
all: \
	rst2ini

rst2ini: $(OBJ1)
	$(COMP) $(FLAGS) $(LIB) $(OBJ1) -o rst2ini

rst2ini.o rst2ini.mod: rst2ini.f90
	$(COMP) $(FLAGS) $(INCL) -c rst2ini.f90

clean: 
	/bin/rm *.o
