DIR_ROOT = ../../../..
include $(DIR_ROOT)/Makefile.cfg

SRC	=  \
	sub_bkw_prop.f90 \
	sub_bkw_action.f90

OBJ	= $(SRC:.f90=.o)

#SUBDIRS	=
FC_FLAGS	+=	-fpp
INCLUDE	+= -I..

%.o : %.f90
	$(FC) $(FC_FLAGS) $(INCLUDE)  $< -c -o $*.o 

%.o : %.f
	$(FC) $(FC_FLAGS) $(INCLUDE) $< -c -o $*.o 

all: build

build: $(OBJ)

clean:
	rm -rf ./*.o 
	rm -rf *.mod 
	rm -rf *~
	rm -rf *genmod.f90