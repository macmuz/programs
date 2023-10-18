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
OBJ1  = cice_his.o
OBJ2  = roms_ave.o
OBJ3  = volume.o
OBJ4  = roms_temp.o
OBJ5  = nemo_temp.o
OBJ6  = nemo_ave.o
	
all: \
	cice_his\
	roms_ave\
	volume\
	roms_temp\
	nemo_temp\
	nemo_ave

cice_his: $(OBJ1)
	$(COMP) $(FLAGS) $(LIB) $(OBJ1) -o cice_his

cice_his.o cice_his.mod: cice_his.f90
	$(COMP) $(FLAGS) $(INCL) -c cice_his.f90

roms_ave: $(OBJ2)
	$(COMP) $(FLAGS) $(LIB) $(OBJ2) -o roms_ave

roms_ave.o roms_aves.mod: roms_ave.f90
	$(COMP) $(FLAGS) $(INCL) -c roms_ave.f90

volume: $(OBJ3)
	$(COMP) $(FLAGS) $(LIB) $(OBJ3) -o volume

volume.o volumes.mod: volume.f90
	$(COMP) $(FLAGS) $(INCL) -c volume.f90

roms_temp: $(OBJ4)
	$(COMP) $(FLAGS) $(LIB) $(OBJ4) -o roms_temp

roms_temp.o roms_temps.mod: roms_temp.f90
	$(COMP) $(FLAGS) $(INCL) -c roms_temp.f90

nemo_temp: $(OBJ5)
	$(COMP) $(FLAGS) $(LIB) $(OBJ5) -o nemo_temp

nemo_temp.o nemo_temps.mod: nemo_temp.f90
	$(COMP) $(FLAGS) $(INCL) -c nemo_temp.f90

nemo_ave: $(OBJ6)
	$(COMP) $(FLAGS) $(LIB) $(OBJ6) -o nemo_ave

nemo_ave.o nemo_aves.mod: nemo_ave.f90
	$(COMP) $(FLAGS) $(INCL) -c nemo_ave.f90

clean: 
	/bin/rm *.o
