FC	=	ifort
F_OPT	=	-warn -O3
L_OPT	=
SRC	=	read_info.f90 read_traj.f90 selc.f90
OBJ	=	$(SRC:.f90=.o)
TARGET	=	selc

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