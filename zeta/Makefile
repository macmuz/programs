#
# Makefile for writing nc file code
# CVS:$Id: makefile,v 1.3 2000/04/19 22:05:59 pwjones Exp $
#

COMP = mpiifort
#COMP = gfortran

FLAGS = -O3 -assume byterecl -convert big_endian -fpic -mcmodel large -qopenmp
NETCDF_DIR = $(shell nf-config --prefix)
NETCDF_INC = -I$(NETCDF_DIR)/include
NETCDF_LIB = -L$(NETCDF_DIR)/lib
LIB  = $(NETCDF_LIB) -lnetcdf -lnetcdff
INCL = $(NETCDF_INC)
SRCDIR  =
EXEDIR  = .
OBJ1  = read_files.o
OBJ2  = repack.o
	
all: \
    read_files\
	repack

read_files: $(OBJ1)
	$(COMP) $(FLAGS) $(LIB) $(OBJ1) -o read_files

read_files.o read_files.mod: read_files.f90
	$(COMP) $(FLAGS) $(INCL) -c read_files.f90

repack: $(OBJ2)
	$(COMP) $(FLAGS) $(LIB) $(OBJ2) -o repack

repack.o repack.mod: repack.f90
	$(COMP) $(FLAGS) $(INCL) -c repack.f90

clean: 
	/bin/rm *.o
