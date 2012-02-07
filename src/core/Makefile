DIR_ROOT = ../..
include $(DIR_ROOT)/Makefile.cfg

SRC	= \
	kepler_near_core.f90 \
	coulomb_tail.f90 \
	dbexp_int.f90 \
	local_minima.f90 \
	rk4.f90 \
	rk4_plot.f90 \
	crf.f90 \
	find_root.f90 \
	im_propagation.f90 \
	spline.f90 \
	action_bak.f90 \
	sub_integrate.f90 \
	sub_traj_ptb_1.f90 \
	sub_cc.f90 \
	propagate.f90

OBJ	= $(SRC:.f90=.o)

DIR_PULSE	=	./pulse
FC_FLAGS	+=	-fpp 
INCLUDE	+=	-I$(DIR_PULSE)



%.o : %.f90
	$(FC) $(FC_FLAGS) $(INCLUDE)  $< -c -o $*.o 

%.o : %.f
	$(FC) $(FC_FLAGS) $(INCLUDE) $< -c -o $*.o 

all: build

build: build_pulse $(OBJ)

build_pulse:
	cd $(DIR_PULSE) && $(MAKE) build

clean:
	cd $(DIR_PULSE) && $(MAKE) $@
	rm -rf ./*.o 
	rm -rf *.mod 
	rm -rf *~
	rm -rf *genmod.f90