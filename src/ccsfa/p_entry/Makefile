DIR_ROOT = ../../..
include $(DIR_ROOT)/Makefile.cfg

SRC		= 	console_slave.f90 console_master.f90

OUT		=	$(DIR_ROOT)/$(EXE_NAME)
DIR_CORE	= 	../core
DIR_MMFF	= 	$(DIR_ROOT)/src/mmff

FC_FLAGS	+= 	-fpp
INCLUDE		+= 	-I. -I$(DIR_CORE) -I$(DIR_MMFF)
LIBS		+= 	-L$(DIR_CORE) -L$(DIR_ROOT)/lib
LK_FLAGS	+= 	-m64 -lmmff -Wl,-rpath,./lib

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