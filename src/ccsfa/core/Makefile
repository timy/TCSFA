DIR_ROOT = ../../..
include $(DIR_ROOT)/Makefile.cfg

SRC	= \
	dbexp_int.f90 \
	Kepler.f90 \
	local_minima.f90 \
	rk4_re.f90 \
	crf.f90 \
	find_root.f90 \
	pulse.f90 \
	im_propagation.f90 \
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