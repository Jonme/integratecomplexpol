# makefile for testing operations on complex numbers with a struct of mpfi structs
# Use make one_z to compile code for integrating 1/z 
#default: all 
#all: complex_func.o Simple_polynomial.o
	#	$(CC) $(CFLAGS) -o testing complex_func.o Simple_polynomial.o $(LIBS)
#one_z: complex_func.o integrate_1_over_z.o
	#	$(CC) $(CFLAGS) -o one_over_z complex_func.o integrate_1_over_z.o $(LIBS)
#complex_func.o: complex_func.c complex_func.h
	#	$(CC) $(CFLAGS) -c complex_func.c
#Simple_polynomial.o: Simple_polynomial.c complex_func.h
	#	$(CC) $(CFLAGS) -c Simple_polynomial.c

CC=gcc
CFLAGS= -g -Wall
LIBS= -lm -lmpfr -lgmp -lmpfi
OBJ= complex_func.o integrate_one_divided_by_z.o testing_functions.o integrate_twoz_div_by_sqrz_plus1.o
DEP=complex_func.h cintegrals.h

# Creates all .o files from all .c files
%.o: %.c $(DEP)
	$(CC) $(CFLAGS) -c -o $@ $< 

# Links and compiles into an executable.
one_z: $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

.PHONY: clean
	
clean:
	$(RM) testing one_z *.o *~



