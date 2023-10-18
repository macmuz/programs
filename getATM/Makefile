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
OBJ1  = getERA.o
OBJ2  = getUERRA.o
OBJ3  = getNEMO.o
OBJ4  = getWAVE.o
	
all: \
	getERA\
	getUERRA\
	getNEMO\
	getWAVE

getERA: $(OBJ1)
	$(COMP) $(FLAGS) $(LIB) $(OBJ1) -o getERA

getERA.o getERA.mod: getERA.f90
	$(COMP) $(FLAGS) $(INCL) -c getERA.f90

getUERRA: $(OBJ2)
	$(COMP) $(FLAGS) $(LIB) $(OBJ2) -o getUERRA

getUERRA.o getUERRA.mod: getUERRA.f90
	$(COMP) $(FLAGS) $(INCL) -c getUERRA.f90

getNEMO: $(OBJ3)
	$(COMP) $(FLAGS) $(LIB) $(OBJ3) -o getNEMO

getNEMO.o getNEMO.mod: getNEMO.f90
	$(COMP) $(FLAGS) $(INCL) -c getNEMO.f90

getWAVE: $(OBJ4)
	$(COMP) $(FLAGS) $(LIB) $(OBJ4) -o getWAVE

getWAVE.o getWAVE.mod: getWAVE.f90
	$(COMP) $(FLAGS) $(INCL) -c getWAVE.f90

clean: 
	/bin/rm *.o
