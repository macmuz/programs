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
OBJ1  = flow_NEMO.o
OBJ2  = flow_ROMS.o
OBJ3  = new_flow_NEMO.o
	
all: \
	flow_NEMO\
	flow_ROMS\
	new_flow_NEMO

flow_NEMO: $(OBJ1)
	$(COMP) $(FLAGS) $(LIB) $(OBJ1) -o flow_NEMO

flow_NEMO.o flow_NEMO.mod: flow_NEMO.f90
	$(COMP) $(FLAGS) $(INCL) -c flow_NEMO.f90

flow_ROMS: $(OBJ2)
	$(COMP) $(FLAGS) $(LIB) $(OBJ2) -o flow_ROMS

flow_ROMS.o flow_ROMS.mod: flow_ROMS.f90
	$(COMP) $(FLAGS) $(INCL) -c flow_ROMS.f90

new_flow_NEMO: $(OBJ3)
	$(COMP) $(FLAGS) $(LIB) $(OBJ3) -o new_flow_NEMO

new_flow_NEMO.o new_flow_NEMO.mod: new_flow_NEMO.f90
	$(COMP) $(FLAGS) $(INCL) -c new_flow_NEMO.f90

clean: 
	/bin/rm *.o
