DIR_ROOT = ../../../..
include $(DIR_ROOT)/Makefile.cfg

SRC	=  \
	sub_ptb_prop.f90 \
	sub_ptb_prop_0.f90 \
	sub_ptb_action.f90

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