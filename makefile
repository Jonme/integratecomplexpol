# makefile for testing operations on complex numbers with a struct of mpfi structs

CC=gcc
CFLAGS= -g -Wall

default: all 

all: complex_func.o Simple_polynomial.o
	$(CC) $(CFLAGS) -o testing complex_func.o Simple_polynomial.o \
	-lm -lmpfr -lgmp -lmpfi

one_z: complex_func.o integrate_1_over_z.o
	$(CC) $(CFLAGS) -o one_over_z complex_func.o integrate_1_over_z.o \
	-lm -lmpfr -lgmp -lmpfi

complex_func.o: complex_func.c complex_func.h
	$(CC) $(CFLAGS) -c complex_func.c
	
Simple_polynomial.o: Simple_polynomial.c complex_func.h
	$(CC) $(CFLAGS) -c Simple_polynomial.c
	
clean:
	$(RM) testing one_over_z *.o *~
