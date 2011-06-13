DIR_ROOT = ../../..
include $(DIR_ROOT)/Makefile.cfg

SRC		= 	data_proc.f90 test.f90

OUT		=	$(DIR_ROOT)/$(EXE_NAME)
DIR_CORE	= 	../core

FC_FLAGS	+= 	-fpp
INCLUDE		+= 	-I. -I$(DIR_CORE)
LIBS		+= 	-L$(DIR_CORE)
LK_FLAGS	+= 	#-m64


%.o : %.f90
	$(FC) $(FC_FLAGS) $(INCLUDE) $< -c -o $@

%.o : %.f
	$(FC) $(FC_FLAGS) $(INCLUDE) $< -c -o $@ 

OBJ_ENTRY	= $(SRC:.f90=.o)
OBJ_CORE	= $(DIR_CORE)/*.o

SUBDIRS = $(DIR_CORE)

all: build

build: head build_core build_entry

head:
	cd .. && python config.py

build_core:
	cd $(DIR_CORE) && $(MAKE) build

build_entry: $(OBJ_ENTRY)
	$(FC) $(LK_FLAGS) $(OBJ_ENTRY) $(OBJ_CORE) $(LIBS) -o $(OUT)

clean:
	cd $(DIR_CORE) && $(MAKE) $@
	$(RM) ./*.o
	$(RM) *.mod 
	$(RM) *~ 

cleanall: clean
	$(RM) ../include/*.h