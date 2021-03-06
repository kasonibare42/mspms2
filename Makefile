# Project: mspms2

CC = gcc

OBJ = mspms2.o random.o erfrc.o rafrc.o nvtnh.o sffrc.o nptnh.o
BIN = mspms2.x
LIBS = mylibtasos/libtasos.a /usr/local/lib/libgfortran.a
#LIBS = mylibtasos/libtasos.a /usr/lib/gcc/i386-redhat-linux/4.1.2/libgfortran.a
CFLAGS = -std=c99
RM = rm -f

.PHONY: all all-before all-after clean clean-custom

all: mspms2.x

clean: clean-custom 
	${RM} $(OBJ) $(BIN)

$(BIN): $(OBJ)
	$(CC) -static $(OBJ) -o mspms2.x $(LIBS) $(CFLAGS) -O2 -lgfortran -lgsl -lgslcblas -lm

mspms2.o: mspms2.c
	$(CC) -c mspms2.c $(CFLAGS)

random.o: random.c 
	$(CC) -c random.c $(CFLAGS)

erfrc.o: erfrc.c
	$(CC) -c erfrc.c $(CFLAGS)

rafrc.o: rafrc.c
	$(CC) -c rafrc.c $(CFLAGS)

nvtnh.o: nvtnh.c
	$(CC) -c nvtnh.c $(CFLAGS)

sffrc.o: sffrc.c
	$(CC) -c sffrc.c $(CFLAGS)
 
nptnh.o: nptnh.c
	$(CC) -c nptnh.c $(CFLAGS)

