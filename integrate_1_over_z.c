/*
Program testing line integrals of simple polynomials in the complex plane,
using Riemann sums. Integrating around square counter-clockwise. 
f(z) = 1 / z
Function has known pole at the origin.

----- Use make one_z to compile code for integrating 1/z -------

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

#define DEBUG 0
#define SUM 1

int main(int argc, char *argv[])
{
	// Check for correct number of inputs..
	if (argc != 3)
	{
		printf("This program takes two inputs, number of points, precision \n");
		return -1;
	}
	// Declare variables and constants.
	cmpfi z1_sum, z2_sum, z3_sum, z4_sum, z1_tmp, z2_tmp, z3_tmp, z4_tmp;
	cmpfi sum, z, dz, dz_tmp, val1, val1c, val4, valN, val0, test_div, test_mul;
	mpfi_t interval, calculations;
	int N, p, left, right, length;
	
	// Initialize variables and constants.
	mpfi_init(interval);
	mpfi_init(calculations);
	
	cInit(test_div);
	cInit(test_mul);
	cInit(z1_sum);
	cInit(z2_sum);
	cInit(z3_sum);
	cInit(z4_sum);
	cInit(z1_tmp);
	cInit(z2_tmp);
	cInit(z3_tmp);
	cInit(z4_tmp);
	cInit(sum);
	cInit(z);
	cInit(dz);
	cInit(dz_tmp);
	cInit(val1);
	cInit(val1c);
	cInit(val4);
	cInit(valN);
	cInit(val0);

	// Collect inputs
	N = atoi(argv[1]); // Number of steps as unsigned int.
	p = atoi(argv[2]); // Computational precision as unsigned int.
	
	// Assign values to constants
	left = -1; // Integral limit
	right = 1;	// Integral limit
	length = right - left; // Interval size
	cSet_ui(val1, 1, 0); // Value 1
	cSet_ui(val1c, 0, 1); // Value i
	cSet_ui(valN, N, 0); // Value N from input
	mpfi_set_ui(interval, length); // Interval size
	mpfi_set_ui(calculations, N); // Number of calculations
 

	
	// Precision
	mpfr_prec_t prec;
	prec = p; // Set prec to p from input
	mpfr_set_default_prec(prec);

	// Set sums to zero
	cSet_ui(z1_sum, 0, 0);
	cSet_ui(z2_sum, 0, 0);
	cSet_ui(z3_sum, 0, 0);
	cSet_ui(z4_sum, 0, 0);
		
	// Calculate and set step size, dz
	mpfi_div(dz->real, interval, calculations);
	mpfi_set_ui(dz->imag, 0);
	
	// First step to create an upper and lower Riemann sum
	cSet_si(z, left, left); // Setting z
	cSet_si(dz_tmp, left, 0); // Creating first temporary step
	cAdd(dz_tmp, dz_tmp, dz); // This is Re(z) + dz

	// Testing multiplication and division
	#if DEBUG
	cDiv(test_div, val1, z);
	printf("Testing division on complex numbers \n");
	cPrint(test_div);
	cMul(test_mul, z, dz);
	printf("Testing multiplication on complex numbers \n");
	cPrint(test_mul);
	#endif
	
	// Setting start of inteval to z and end of interval to z + dz
	mpfi_union(z->real, z->real, dz_tmp->real);

	

	#if DEBUG
	// Print
	cPrint(z);
	printf("Before loops\n");
	#endif
	
	// Integral 1
	for (int i = 0; i < N; i++)
	{
		cDiv(z1_tmp, val1, z); // This is our integrand 1/z
		cMul(z1_tmp, z1_tmp, dz); // Multiplying with step size dz
		cAdd(z1_sum, z1_sum, z1_tmp); // Summing the results
		cAdd(z, z, dz); // Advancing z to next interval
	}
	#if SUM	
	// Print
	printf("Sum 1: \n");
	cPrint(z1_sum);
	#endif
	
	// Integral 2
	// Calculate and set step size
	mpfi_div(dz->imag, interval, calculations);
	mpfi_set_ui(dz->real, 0);
	
	// First step
	cSet_si(z, right, left);
	cSet_si(dz_tmp, right, left);
	cAdd(dz_tmp, dz_tmp, dz);

	// Setting start and end of first interval.
	mpfi_union(z->imag, z->imag, dz_tmp->imag);
	for (int i = 0; i < N; i++)
	{
		cDiv(z2_tmp, val1, z); // This is our integrand 1/z
		// Debugging imaginary part of the integral
		#if DEBUG
		if (i == 0)
		{
		cPrint(z2_tmp);
		}
		#endif
		cMul(z2_tmp, z2_tmp, dz); // Multiplying with step size dz
		// Debugging imaginary part of the integral
		#if DEBUG
		if (i == 0)
		{
		cPrint(z2_tmp);
		}
		#endif
		cAdd(z2_sum, z2_sum, z2_tmp); // Summing the results
		cAdd(z, z, dz); // Advanzing z to next interval
	}
	#if SUM
	// Print
	printf("Sum 2: \n");
	cPrint(z2_sum);
	#endif
	
	// Integral 3
	// Calculate and set step size
	mpfi_div(dz->real, interval, calculations);
	mpfi_set_ui(dz->imag, 0);
	
	// Change the sign of dz to flip the bounds of the integral
	cChange_sign(dz, 3); // Inputs cmpfi struct and flag what sign to change 3 for all
	
	// First step
	cSet_si(z, right, right);
	cSet_si(dz_tmp, right, 0);
	cAdd(dz_tmp, dz_tmp, dz);
	
	// Setting start and end of first interval
	mpfi_union(z->real, z->real, dz_tmp->real);
	for (int i = 0; i < N; i++)
	{
		cDiv(z3_tmp, val1, z); // This is our integrand 1/z
		cMul(z3_tmp, z3_tmp, dz); // Multiplying with step size dz
		cAdd(z3_sum, z3_sum, z3_tmp); // Summing the results
		cAdd(z, z, dz); // Advancing z to next interval
	}
	#if SUM
	// Print
	printf("Sum 3: \n");
	cPrint(z3_sum);
	#endif

	// Integral 4
	// Calculate and set step size
	mpfi_div(dz->imag, interval, calculations);
	mpfi_set_ui(dz->real, 0);
	
	// Change the sign of dz to flip the bounds of the integral
	cChange_sign(dz, 3);
	
	// First step
	cSet_si(z, left, right);
	cSet_si(dz_tmp, 0, right);
	cAdd(dz_tmp, dz_tmp, dz);

	// Setting start and end of first interval
	mpfi_union(z->imag, z->imag, dz_tmp->imag);
	for (int i = 0; i < N; i++)
	{
		cDiv(z4_tmp, val1, z); // This is our integrand 1/z
		cMul(z4_tmp, z4_tmp, dz); // Multiplying with step size dz
		cAdd(z4_sum, z4_sum, z4_tmp); // Summing the results
		cAdd(z, z, dz); // Advancing z to next interval
	}
	#if SUM
	// Print
	printf("Sum 4: \n");
	cPrint(z4_sum);
	#endif
	
	// Adding sums together
	cAdd(sum, z1_sum, z2_sum);
	cAdd(sum, sum, z3_sum);
	cAdd(sum, sum, z4_sum);
	
	// Printing total sum
	printf("Total: \n");
	cPrint(sum);

	// Clearing used variables and constants
	mpfi_clear(interval);
	mpfi_clear(calculations);
	cClear(z1_sum);
	cClear(z2_sum);
	cClear(z3_sum);
	cClear(z4_sum);
	cClear(z1_tmp);
	cClear(z2_tmp);
	cClear(z3_tmp);
	cClear(z4_tmp);
	cClear(sum);
	cClear(z);
	cClear(dz);
	cClear(dz_tmp);
	cClear(val1c);
	cClear(val1);
	cClear(val4);
	cClear(valN);
	cClear(val0);	
	
}
