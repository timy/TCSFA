DIR_ROOT = ../../..
include $(DIR_ROOT)/Makefile.cfg

SRC	= \
	sub_p0.f90 \
	sub_integrate.f90 \
	sub_traj_ptb_1.f90 \
	sub_cc.f90 \
	sub_action_0.f90

OBJ	= $(SRC:.f90=.o)

#DIR_PULSE	=	./pulse
FC_FLAGS	+=	-fpp
INCLUDE	+=

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