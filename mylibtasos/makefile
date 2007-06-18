.SUFFIXES: .f90 .c .F90 
SHELL = /bin/sh

###############################
# PGI Compilers               # 
###############################
#FFLAGS = -O2 -Mextend -Mpreprocess -Munroll -Mvect -Mdepchk -module ./Modules
#FFLAGS = -fast -tpp7 -Bstatic
#FFLAGS = -O2 -Mextend -Minfo -Mneginfo=loop -Mpreprocess -Mbounds -Mchkptr -Mchkstk -Mdclchk -Bstatic -module ./Modules
#CFLAGS = -fast -Bstatic -tpp7
#FC = pgf90
#CC = pgcc
###############################
# GNU Compilers               #
###############################
FC = gfortran
FFLAGS  = -O3 -fno-automatic 

OBJECTS = \
functions.o\
global.o\
cubicauxiliary.o\
initialize.o \
cubicinterpolation.o\
grid.o\
forces.o\
pbc.o\

SUFFIX = .f90 

TARGET = libtasos.a

$(TARGET): $(OBJECTS) 
	ar cr libtasos.a $(OBJECTS) 
	ranlib libtasos.a

.f90.o:
	$(FC) $(FFLAGS) -c $*$(SUFFIX)

clean :
	rm -f *.o $(TARGET) *.mod 

$(OBJECTS): global.f90

%.o: %.mod
