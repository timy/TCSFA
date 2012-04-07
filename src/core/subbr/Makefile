DIR_ROOT = ../../..
include $(DIR_ROOT)/Makefile.cfg

SRC	=  \
	sub_p0.f90 \
	sub_rk4.f90 \
	sub_integrate.f90 \
	sub_traj_0.f90 \

OBJ	= $(SRC:.f90=.o)

SUBDIRS		=	./ptb ./bkw ./bkw_eps
FC_FLAGS	+=	-fpp
INCLUDE	+=

%.o : %.f90
	$(FC) $(FC_FLAGS) $(INCLUDE)  $< -c -o $*.o 

%.o : %.f
	$(FC) $(FC_FLAGS) $(INCLUDE) $< -c -o $*.o 

all: build

build: $(OBJ) build_subdir

clean: clean_subdir
	rm -rf ./*.o 
	rm -rf *.mod 
	rm -rf *~
	rm -rf *genmod.f90

build_subdir: $(SUBDIRS)
	for subdir in $(SUBDIRS); do (cd $${subdir}; $(MAKE) build); done

clean_subdir: $(SUBDIRS)
	for subdir in $(SUBDIRS); do (cd $${subdir}; $(MAKE) clean); done