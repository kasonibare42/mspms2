# Project mspms-2

CC = gcc

OBJ = mspms-2.o random.o
BIN = mspms-2.x

all: mspms-2.x

$(BIN): $(OBJ)
	$(CC) $(OBJ) -o mspms-2.x -lm

mspms-2.o: mspms-2.c
	$(CC) -c mspms-2.c

random.o: random.c 
	$(CC) -c random.c

