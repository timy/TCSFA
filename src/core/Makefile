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
	pulse.f90 \
	pulse_ramp.f90 \
	pulse_const.f90 \
	pulse_sin2_old.f90 \
	im_propagation.f90 \
	spline.f90 \
	action_bak.f90 \
	sub_integrate.f90 \
	sub_traj_ptb_1.f90 \
	sub_cc.f90 \
	propagate.f90

FC_FLAGS	+= -fpp


%.o : %.f90
	$(FC) $(FC_FLAGS) $< -c -o $*.o 

%.o : %.f
	$(FC) $(FC_FLAGS) $< -c -o $*.o 

OBJ	= $(SRC:.f90=.o)


all: clean build

build: $(OBJ)

clean:
	rm -rf ./*.o 
	rm -rf *.mod 
	rm -rf *~
	rm -rf *genmod.f90