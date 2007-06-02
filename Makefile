# Project mspms-2

CC = gcc

OBJ = mspms-2.o random.o erfrc.o rafrc.o nvtnh.o
BIN = mspms-2.x
RM = rm -f

.PHONY: all all-before all-after clean clean-custom

all: mspms-2.x

clean: clean-custom 
	${RM} $(OBJ) $(BIN)

$(BIN): $(OBJ)
	$(CC) $(OBJ) -o mspms-2.x -lm -std=c99 

mspms-2.o: mspms-2.c
	$(CC) -c mspms-2.c -std=c99

random.o: random.c 
	$(CC) -c random.c -std=c99

erfrc.o: erfrc.c
	$(CC) -c erfrc.c -std=c99

rafrc.o: rafrc.c
	$(CC) -c rafrc.c -std=c99

nvtnh.o: nvtnh.c
	$(CC) -c nvtnh.c -std=c99
