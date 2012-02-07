FC	=	ifort
F_OPT	=	-warn -O3 -fpp
L_OPT	=
SRC	=	count_lines.f90 read_selc.f90 stat.f90
OBJ	=	$(SRC:.f90=.o)
TARGET	=	stat

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
	rm -f ../../dat/hist_*.dat