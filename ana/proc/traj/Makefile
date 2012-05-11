DIR_ENTRY	=	main
DIR_ROOT	=	../../..
include $(DIR_ROOT)/Makefile.cfg

PARALLEL	=	0

SUBDIRS = $(DIR_BASE)/$(DIR_MAIN)/$(DIR_ENTRY) 

include $(DIR_ROOT)/Makefile.cfg

SRC			=	main.f90

OUT			=	$(EXE_NAME)
DIR_CORE	= 	$(DIR_ROOT)/src/core
DIR_PULSE	=	$(DIR_CORE)/pulse
DIR_SUBBR	=	$(DIR_CORE)/subbr
DIR_SUBBR_PTB	=	$(DIR_SUBBR)/ptb
DIR_SUBBR_BKW	=	$(DIR_SUBBR)/bkw
DIR_SUBBR_BKW_EPS	=	$(DIR_SUBBR)/bkw_eps

FC_FLAGS	+= 	-fpp
INCLUDE	+= 	-I. -I$(DIR_CORE) -I$(DIR_PULSE) -I$(DIR_SUBBR)
LIBS		+=
LK_FLAGS	+= 	#-m64

%.o : %.f90
	$(FC) $(FC_FLAGS) $(INCLUDE) $< -c -o $@

%.o : %.f
	$(FC) $(FC_FLAGS) $(INCLUDE) $< -c -o $@ 

OBJ_ENTRY	=	$(SRC:.f90=.o)
OBJ_CORE	+=	$(DIR_CORE)/*.o $(DIR_SUBBR)/*.o $(DIR_PULSE)/*.o \
	$(DIR_SUBBR_PTB)/*.o $(DIR_SUBBR_BKW)/*.o $(DIR_SUBBR_BKW_EPS)/*.o

all: build

build: head build_core build_entry

head:
	cd $(DIR_ROOT)/src && python config.py traj

build_core:
	cd $(DIR_CORE) && $(MAKE) build

build_entry: $(OBJ_ENTRY)
	$(FC) $(LK_FLAGS) $(OBJ_ENTRY) $(OBJ_CORE) $(LIBS) -o $(OUT)	

clean:
	cd $(DIR_CORE) && $(MAKE) $@
	$(RM) $(DIR_ROOT)/src/include/*.h
	$(RM) ./*.o
	$(RM) *.mod 
	$(RM) *~ 
	$(RM) *genmod.f90
	$(RM) ../../dat/traj/*.dat
