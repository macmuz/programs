#
# Makefile for writing nc file code
# CVS:$Id: makefile,v 1.3 2000/04/19 22:05:59 pwjones Exp $
#

COMP = mpiifort

FLAGS = -O3 -assume byterecl -convert big_endian -fpic -mcmodel large -qopenmp
NETCDF_DIR = $(shell nf-config --prefix)
NETCDF_INC = -I$(NETCDF_DIR)/include
NETCDF_LIB = -L$(NETCDF_DIR)/lib
#LIB  = $(NETCDF_LIB) -lnetcdf -lnetcdff -L/apl/tryton/intel/2019.05/compilers_and_libraries/linux/mpi/intel64/lib
#INCL = $(NETCDF_INC) -I/apl/tryton/intel/2019.05/compilers_and_libraries/linux/mpi/intel64/include
LIB  = $(NETCDF_LIB) -lnetcdf -lnetcdff
INCL = $(NETCDF_INC)
SRCDIR  =
EXEDIR  = .
OBJ1  = mask.o
OBJ2  = interpolation.o
OBJ3  = coast.o
OBJ4  = linear_interp.o
	
all: \
	mask\
	interpolation\
	coast\
	linear_interp

mask: $(OBJ1)
	$(COMP) $(FLAGS) $(LIB) $(OBJ1) -o mask

mask.o mask.mod: mask.f90
	$(COMP) $(FLAGS) $(INCL) -c mask.f90

interpolation: $(OBJ2)
	$(COMP) $(FLAGS) $(LIB) $(OBJ2) -o interpolation

interpolation.o interpolation.mod: interpolation.f90
	$(COMP) $(FLAGS) $(INCL) -c interpolation.f90

coast: $(OBJ3)
	$(COMP) $(FLAGS) $(LIB) $(OBJ3) -o coast

coast.o coast.mod: coast.f90
	$(COMP) $(FLAGS) $(INCL) -c coast.f90

linear_interp: $(OBJ4)
	$(COMP) $(FLAGS) $(LIB) $(OBJ4) -o linear_interp

linear_interp.o linear_interp.mod: linear_interp.f90
	$(COMP) $(FLAGS) $(INCL) -c linear_interp.f90

clean: 
	/bin/rm *.o
