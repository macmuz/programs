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
OBJ2  = ROMS_bc.o
OBJ3 = ROMS_bc_spline.o
OBJ4 = ROMS_bc_obs_zeta.o
OBJ5 = robert.o
OBJ6 = ROMS_bc_bar.o
OBJ7 = new_ROMS_bc_bar.o
OBJ8 = old_ROMS_bc_bar.o
	
all: \
	ROMS_bc\
	ROMS_bc_spline\
	ROMS_bc_obs_zeta\
	robert\
	ROMS_bc_bar\
	new_ROMS_bc_bar\
	old_ROMS_bc_bar	

ROMS_bc: $(OBJ2)
	$(COMP) $(FLAGS) $(LIB) $(OBJ2) -o ROMS_bc

ROMS_bc.o ROMS_bc.mod: ROMS_bc.f90
	$(COMP) $(FLAGS) $(INCL) -c ROMS_bc.f90

ROMS_bc_spline: $(OBJ3)
	$(COMP) $(FLAGS) $(LIB) $(OBJ3) -o ROMS_bc_spline

ROMS_bc_spline.o ROMS_bc_spline.mod: ROMS_bc_spline.f90
	$(COMP) $(FLAGS) $(INCL) -c ROMS_bc_spline.f90

ROMS_bc_obs_zeta: $(OBJ4)
	$(COMP) $(FLAGS) $(LIB) $(OBJ4) -o ROMS_bc_obs_zeta

ROMS_bc_obs_zeta.o ROMS_bc_obs_zeta.mod: ROMS_bc_obs_zeta.f90
	$(COMP) $(FLAGS) $(INCL) -c ROMS_bc_obs_zeta.f90

robert: $(OBJ5)
	$(COMP) $(FLAGS) $(LIB) $(OBJ5) -o robert

robert.o robert.mod: robert.f90
	$(COMP) $(FLAGS) $(INCL) -c robert.f90

ROMS_bc_bar: $(OBJ6)
	$(COMP) $(FLAGS) $(LIB) $(OBJ6) -o ROMS_bc_bar

ROMS_bc_bar.o ROMS_bc_bar.mod: ROMS_bc_bar.f90
	$(COMP) $(FLAGS) $(INCL) -c ROMS_bc_bar.f90

old_ROMS_bc_bar: $(OBJ8)
	$(COMP) $(FLAGS) $(LIB) $(OBJ8) -o old_ROMS_bc_bar

old_ROMS_bc_bar.o old_ROMS_bc_bar.mod: old_ROMS_bc_bar.f90
	$(COMP) $(FLAGS) $(INCL) -c old_ROMS_bc_bar.f90

clean: 
	/bin/rm *.o
