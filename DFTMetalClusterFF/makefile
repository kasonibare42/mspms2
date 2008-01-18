.SUFFIXES: .f
SHELL = /bin/sh

FC = g77
FFLAGS  = -O3 

OBJECTS = \
ffieldCu.o\

SUFFIX = .f

TARGET = DFTMetalClusterFF.a

$(TARGET): $(OBJECTS) 
	ar cr DFTMetalClusterFF.a $(OBJECTS) 
	ranlib DFTMetalClusterFF.a

.f.o:
	$(FC) $(FFLAGS) -c $*$(SUFFIX)

clean :
	rm -f *.o $(TARGET) *.mod 

$(OBJECTS): ffieldCu.f

%.o: %.mod

