/*
Program testing line integrals of simple polynomials in the complex plane,
using Riemann sums. Integrating around squares counter-clockwise. 
f(z) = z² + 1.
Function has known zeros in +i and -i.
Integrating ( f' / f ) * dz.
f'/f = (2*z) / (z² + 1) 
Author: @JM.
*/
#include <stdio.h>
#include "mpfi.h"
#include "mpfi_io.h"
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include <mpfr.h>
#include "complex_func.h"
int main(int argc, char *argv[])
{
	// Check for correct number of inputs..
	if (argc != 3)
	{
		printf("This program takes two inputs, number of points, precision \n");
		return -1;
	}
	
	cmpfi z1_sum, z2_sum, z3_sum, z4_sum, z1_tmp, z2_tmp, z3_tmp, z4_tmp;
	cmpfi sum, z, dz, dz_tmp, val1, val4, valN, val0;
	mpfi_t interval, calculations;
	
	mpfi_init(interval);
	mpfi_init(calculations);
	
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
	cInit(val4);
	cInit(valN);
	cInit(val0);

	// Declare inputs and constants	
	int N, p, left, right, length;
	left = -1;
	right = 1;
	length = 2;
	// Collect inputs
	N = atoi(argv[1]); // Number of steps as unsigned int.
	p = atoi(argv[2]); // Computational precision as unsigned int.
	
	// Constants
	cSet_ui(val1, 1, 0);
	
	// Precision
	mpfr_prec_t prec;
	prec = p;
	mpfr_set_default_prec(prec);
	
	// Step size
	mpfi_set_ui(interval, length);
	mpfi_set_ui(calculations, N);
	cSet_ui(valN, N, 0);
	mpfi_div(dz->real, interval, calculations);
	mpfi_set_d(dz->imag, 0.0);
	
	
	// First step
	cSet_si(z, left, left);
	cSet_si(dz_tmp, left, 0);
	cAdd(dz_tmp, dz_tmp, dz);
	//mpfi_put(z->real, dz_tmp->real);
	//cUnion(z, z, dz_tmp);
	mpfi_union(z->real, z->real, dz_tmp->real);
	cSet_ui(z1_sum, 0, 0);
	cSet_ui(z2_sum, 0, 0);
	cSet_ui(z3_sum, 0, 0);
	cSet_ui(z4_sum, 0, 0);
	cPrint(z);
	printf("Before loops\n");
	// Integral 1
	for (int i = 0; i < N; i++)
	{
		cDiv(z1_tmp, val1, z);
		cMul(z1_tmp, z1_tmp, dz);
		cAdd(z1_sum, z1_sum, z1_tmp);
		cAdd(z, z, dz);
		//cPrint(z);
	}	
	printf("Sum 1: \n");
	cPrint(z1_sum);
	// Integral 2
	mpfi_div(dz->imag, interval, calculations);
	mpfi_set_d(dz->real, 0.0);	
	cSet_si(z, right, left);
	cSet_si(dz_tmp, 0, left);
	cAdd(dz_tmp, dz_tmp, dz);
	//mpfi_put(z->imag, dz_tmp->imag);
	mpfi_union(z->imag, z->imag, dz_tmp->imag);
	
	for (int i = 0; i < N; i++)
	{
		cDiv(z2_tmp, val1, z);
		cMul(z2_tmp, z2_tmp, dz);
		cAdd(z2_sum, z2_sum, z2_tmp);
		cAdd(z, z, dz);
		//cPrint(z);
	}
	printf("Sum 2: \n");
	cPrint(z2_sum);
	
	
	// Integral 3
	mpfi_div(dz->real, interval, calculations);
	mpfi_set_d(dz->imag, 0.0);
	cChange_sign(dz, 1);
	cSet_si(z, right, right);
	cSet_si(dz_tmp, right, 0);
	cAdd(dz_tmp, dz_tmp, dz);
	//mpfi_put(z->real, dz_tmp->real);
	//cUnion(z, z, dz_tmp);
	mpfi_union(z->real, z->real, dz_tmp->real);
	for (int i = 0; i < N; i++)
	{
		cDiv(z3_tmp, val1, z);
		cMul(z3_tmp, z3_tmp, dz);
		cAdd(z3_sum, z3_sum, z3_tmp);
		cAdd(z, z, dz);
		//cPrint(z);
	}
	printf("Sum 3: \n");
	cPrint(z3_sum);
	// Integral 4
	mpfi_div(dz->imag, interval, calculations);
	mpfi_set_d(dz->real, 0.0);
	cChange_sign(dz, 2);
	cSet_si(z, left, right);
	cSet_si(dz_tmp, 0, right);
	cAdd(dz_tmp, dz_tmp, dz);
	//mpfi_put(z->imag, dz_tmp->imag);
	//cUnion(z, z, dz_tmp);
	mpfi_union(z->imag, z->imag, dz_tmp->imag);
	for (int i = 0; i < N; i++)
	{
		cDiv(z4_tmp, val1, z);
		cMul(z4_tmp, z4_tmp, dz);
		cAdd(z4_sum, z4_sum, z4_tmp);
		cAdd(z, z, dz);
		//cPrint(z);
	}
	printf("Sum 4: \n");
	cPrint(z4_sum);

	
	cAdd(sum, z1_sum, z2_sum);
	cAdd(sum, sum, z3_sum);
	cAdd(sum, sum, z4_sum);
	
	printf("Total: \n");
	cPrint(sum);

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
	cClear(val1);
	cClear(val4);
	cClear(valN);
	cClear(val0);	
	
}
