DIR_ROOT = ../..
include $(DIR_ROOT)/Makefile.cfg

MODULE_NAME	= libmmff
OUTPUT		= $(DIR_ROOT)/lib/$(MODULE_NAME).so
OUTPUT_STATIC	= $(MODULE_NAME).a

SRC_LIBRARY 	= \
	write_log.f90 \
	mmff_util.f90 \
	mmff_data.f90 \
	mmff_base_dbl.f90 \
	mmff_base_dcp.f90 \
	mmff_base_int.f90 \
	mmff_workload.f90 \
	mmff_core_dbl.f90 \
	mmff_core_dcp.f90 \
	mmff_core_int.f90 \
	mmff_toolkit.f90 \
	mmff_main.f90

FORTRAN_FLAGS 	:= $(FORTRAN_FLAGS) -fPIC -I./src
LD_FLAGS	:= $(LD_FLAGS) -fPIC

%.o : %.f90
	$(FC) $(FORTRAN_FLAGS) $< -c -o $*.o 

%.o : %.f
	$(FC) $(FORTRAN_FLAGS) $< -c -o $*.o 

OBJ_LIBRARY	= $(SRC_LIBRARY:.f90=.o)

all: lib
lib: $(OUTPUT)

$(OUTPUT): $(OBJ_LIBRARY) 
	$(LD) $(OBJ_LIBRARY) $(LD_FLAGS) $(LD_SHARED) $(LIB_DIR) $(LIBS) -o $@ 

cleanall:
	$(RM) *.o 
	$(RM) *~ 
	$(RM) *.mod 
	$(RM) $(OUTPUT)