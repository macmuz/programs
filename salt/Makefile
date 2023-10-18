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
OBJ1  = salt_NEMO.o
OBJ2  = dens_test.o
OBJ3  = salt_ROMS.o
OBJ4  = new_salt_NEMO.o
	
all: \
	salt_NEMO\
	dens_test\
	salt_ROMS\
	new_salt_NEMO

salt_NEMO: $(OBJ1)
	$(COMP) $(FLAGS) $(LIB) $(OBJ1) -o salt_NEMO

salt_NEMO.o salt_NEMO.mod: salt_NEMO.f90
	$(COMP) $(FLAGS) $(INCL) -c salt_NEMO.f90

dens_test: $(OBJ2)
	$(COMP) $(FLAGS) $(LIB) $(OBJ2) -o dens_test

dens_test.o dens_test.mod: dens_test.f90
	$(COMP) $(FLAGS) $(INCL) -c dens_test.f90

salt_ROMS: $(OBJ3)
	$(COMP) $(FLAGS) $(LIB) $(OBJ3) -o salt_ROMS

salt_ROMS.o salt_ROMS.mod: salt_ROMS.f90
	$(COMP) $(FLAGS) $(INCL) -c salt_ROMS.f90

new_salt_NEMO: $(OBJ4)
	$(COMP) $(FLAGS) $(LIB) $(OBJ4) -o new_salt_NEMO

new_salt_NEMO.o new_salt_NEMO.mod: new_salt_NEMO.f90
	$(COMP) $(FLAGS) $(INCL) -c new_salt_NEMO.f90

clean: 
	/bin/rm *.o
