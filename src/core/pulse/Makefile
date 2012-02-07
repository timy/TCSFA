DIR_ROOT = ../../..
include $(DIR_ROOT)/Makefile.cfg

SRC	= \
	pulse.f90 \
	pulse_ramp.f90 \
#	pulse_sin2_old.f90

FC_FLAGS	+= -fpp


%.o : %.f90
	$(FC) $(FC_FLAGS) $< -c -o $*.o 

%.o : %.f
	$(FC) $(FC_FLAGS) $< -c -o $*.o 

OBJ	= $(SRC:.f90=.o)


all: build

build: $(OBJ)

clean:
	rm -rf ./*.o 
	rm -rf *.mod 
	rm -rf *~
	rm -rf *genmod.f90