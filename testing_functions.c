/*
Function calculating the line integral of f(z) = 1 / z on a closed square 
integrand has known pole at the origin.


Author: @Jonme.
*/
#include <stdio.h>
#include "mpfi.h"
#include "mpfi_io.h"
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include <mpfr.h>
#include "complex_func.h"
#include "cintegrals.h"
#define DEBUG 0
#define SUM 1

int main(int argc, char *argv[])
{
	if (argc != 5)
	{
		printf("Correct inputs are: number of pts, precision,\
		integral start, end\n");
		 return -1;
	}
	// Declaring variables.
	cmpfi sum1, sum2;
	int N, p;
	double start, end;
	
	// Collecting inputs and assigning to variables
	N = atoi(argv[1]);
	p = atoi(argv[2]);
	start = atof(argv[3]);
	end = atof(argv[4]);
	
	// Initializing summation variable
	cInit(sum1);
	cInit(sum2);
	
	// Setting precision
	mpfr_prec_t prec;
	prec = p;
	mpfr_set_default_prec(prec);
	
	// Calling integral function
	cint_one_div_z(sum1, start, end, N);
	cint_twoz_div_sqrz_plus1(sum2, start, end, N);
	
	// Printing results
	printf("This is the result of integrating 1/z: \n");
	cPrint(sum1);
	
	
	printf("This is the result of integrating 2z/(z^2 + 1): \n");
	cPrint(sum2);
	
	// Clearing
	cClear(sum1);
	cClear(sum2);
}

