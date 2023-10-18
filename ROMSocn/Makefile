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
OBJ1  = ROMSocn.o
OBJ2 = OCNforcing.o
	
all: \
	ROMSocn\
	OCNforcing

ROMSocn: $(OBJ1)
	$(COMP) $(FLAGS) $(LIB) $(OBJ1) -o ROMSocn

ROMSocn.o ROMSocn.mod: ROMSocn.f90
	$(COMP) $(FLAGS) $(INCL) -c ROMSocn.f90

OCNforcing: $(OBJ2)
	$(COMP) $(FLAGS) $(LIB) $(OBJ2) -o OCNforcing

OCNforcing.o OCNforcing.mod: OCNforcing.f90
	$(COMP) $(FLAGS) $(INCL) -c OCNforcing.f90

clean: 
	/bin/rm *.o
