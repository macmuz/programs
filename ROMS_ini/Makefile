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
OBJ1  = ROMS_rm_bubble.o
OBJ2  = ROMS_ini.o
OBJ3  = ROMS_ini_ebaltic.o
OBJ4  = create1992.o
OBJ5  = new_ROMS_ini.o
	
all: \
	ROMS_ini\
	ROMS_rm_bubble\
	ROMS_ini_ebaltic\
	create1992\
	new_ROMS_ini

ROMS_rm_bubble: $(OBJ1)
	$(COMP) $(FLAGS) $(LIB) $(OBJ1) -o ROMS_rm_bubble

ROMS_rm_bubble.o ROMS_rm_bubble.mod: ROMS_rm_bubble.f90
	$(COMP) $(FLAGS) $(INCL) -c ROMS_rm_bubble.f90

ROMS_ini: $(OBJ2)
	$(COMP) $(FLAGS) $(LIB) $(OBJ2) -o ROMS_ini

ROMS_ini.o ROMS_ini.mod: ROMS_ini.f90
	$(COMP) $(FLAGS) $(INCL) -c ROMS_ini.f90

ROMS_ini_ebaltic: $(OBJ3)
	$(COMP) $(FLAGS) $(LIB) $(OBJ3) -o ROMS_ini_ebaltic

ROMS_ini_ebaltic.o ROMS_ini_ebaltic.mod: ROMS_ini_ebaltic.f90
	$(COMP) $(FLAGS) $(INCL) -c ROMS_ini_ebaltic.f90

create1992: $(OBJ4)
	$(COMP) $(FLAGS) $(LIB) $(OBJ4) -o create1992

create1992.o create1992.mod: create1992.f90
	$(COMP) $(FLAGS) $(INCL) -c create1992.f90

new_ROMS_ini: $(OBJ5)
	$(COMP) $(FLAGS) $(LIB) $(OBJ5) -o new_ROMS_ini

new_ROMS_ini.o new_ROMS_ini.mod: new_ROMS_ini.f90
	$(COMP) $(FLAGS) $(INCL) -c new_ROMS_ini.f90

ROMS_ini_ebaltic.o ROMS_ini_ebaltic.mod: ROMS_ini_ebaltic.f90
clean: 
	/bin/rm *.o
