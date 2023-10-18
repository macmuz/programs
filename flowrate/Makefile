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
OBJ1  = flowrate.o
OBJ2  = NEMO_sect.o
OBJ3  = sect.o
OBJ4  = NEMO_flow.o
OBJ5  = ROMS_flow.o
	
all: \
	flowrate\
	NEMO_sect\
	sect\
	NEMO_flow\
	ROMS_flow

flowrate: $(OBJ1)
	$(COMP) $(FLAGS) $(LIB) $(OBJ1) -o flowrate

flowrate.o flowrate.mod: flowrate.f90
	$(COMP) $(FLAGS) $(INCL) -c flowrate.f90

NEMO_sect: $(OBJ2)
	$(COMP) $(FLAGS) $(LIB) $(OBJ2) -o NEMO_sect

NEMO_sect.o NEMO_sect.mod: NEMO_sect.f90
	$(COMP) $(FLAGS) $(INCL) -c NEMO_sect.f90

sect: $(OBJ3)
	$(COMP) $(FLAGS) $(LIB) $(OBJ3) -o sect

sect.o sect.mod: sect.f90
	$(COMP) $(FLAGS) $(INCL) -c sect.f90

NEMO_flow: $(OBJ4)
	$(COMP) $(FLAGS) $(LIB) $(OBJ4) -o NEMO_flow

NEMO_flow.o NEMO_flow.mod: NEMO_flow.f90
	$(COMP) $(FLAGS) $(INCL) -c NEMO_flow.f90

ROMS_flow: $(OBJ5)
	$(COMP) $(FLAGS) $(LIB) $(OBJ5) -o ROMS_flow

ROMS_flow.o ROMS_flow.mod: ROMS_flow.f90
	$(COMP) $(FLAGS) $(INCL) -c ROMS_flow.f90

clean: 
	/bin/rm *.o
