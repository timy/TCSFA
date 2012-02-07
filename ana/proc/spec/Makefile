FC	=	/usr/local/Packages/mvapich2-tcp-intel-1.2rc2/bin/mpif90	
F_OPT	=	-O3 -fpp
L_OPT	=
CORE_DIR = 	../../../src/core/
SRC	=	read_info.f90 read_traj.f90 calc_w.f90 main_single.f90
OBJ	=	$(SRC:.f90=.o)
TARGET	=	spec

%.o : %.f90
	$(FC) $(F_OPT) -c $< -o $*.o

all : $(TARGET)

$(TARGET) : $(OBJ)
	$(FC) $(L_OPT) $(OBJ) -o $@

clean :
	rm -f $(OBJ)
	rm -f $(TARGET)
	rm -f *.mod
	rm -f *genmod.f90
	rm -f *~

run :
#	qsub run.sh
	rm -f $(OBJ)
	qsub job.sh -q opt.q