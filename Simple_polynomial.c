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
	// Declare variables
	mpfi_t zi, left, right, first;	
	// Integral 1
	cmpfi z1, tmp_z1, tmp_z2, dz1, z1_sum, z2_sum, z3_sum, z4_sum, cval2;
	cmpfi z1_sum_tmp, z2_sum_tmp, z3_sum_tmp, z4_sum_tmp, sum_total, cval1;
	
	
	// Declare constants
	mpfi_t val1, val2, ilength, valN, dz;
	int length;	
	
	// Initialize declarations
	mpfi_init(val1);
	mpfi_init(val2);
	mpfi_init(ilength);
	mpfi_init(valN);
	mpfi_init(dz);
	mpfi_init(zi);
	mpfi_init(first);
	mpfi_init(left);
	mpfi_init(right);
	cInit(z1);
	cInit(dz1);
	cInit(z1_sum);
	cInit(z1_sum_tmp);
	cInit(tmp_z1);
	cInit(tmp_z2);
	cInit(cval1);
	cInit(cval2);
	cInit(z2_sum);
	cInit(z2_sum_tmp);
	cInit(z3_sum);
	cInit(z3_sum_tmp);
	cInit(z4_sum);
	cInit(z4_sum_tmp);
	cInit(sum_total);
	
	
	// Declare inputs
	int N, p;
	long int startval;
	// Collect inputs
	N = atoi(argv[1]); // Number of steps as unsigned int.
	p = atoi(argv[2]); // Computational precision as unsigned int.
	mpfi_set_ui(valN, N);
	

	
	// Set precision
	mpfr_prec_t prec;
	prec = p;
	mpfr_set_default_prec(prec);
	
	
	// Give values to constants
	mpfi_set_ui(val1, 1);
	mpfi_set_ui(val2, 2);
	cSet_ui(cval1, 1, 0);
	cSet_ui(cval2, 2, 0);

	startval = -2;
	
	// Create intervals
	mpfi_set_si(left, startval);
	mpfi_set_si(right, 2);
	length = 4;
	mpfi_set_ui(ilength, length);
	
	
	// Create step size from inputs and interval size
	mpfi_div(dz, ilength, valN);
	
	
	
	// Set starting positions and values for integral 1
	mpfi_set(zi, left);
	mpfi_add(first, zi, dz);
	mpfi_put(zi, first);

	mpfi_set_si(z1->imag, -2);
	mpfi_set(z1->real, zi);
	mpfi_set(dz1->real, dz);
	mpfi_set_ui(dz1->imag, 0);
	//cSet(cval2, val2, val2);
	mpfi_set_d(z1_sum->real, 0.0);
	mpfi_set_d(z1_sum->imag, 0.0);
	// Loop for integral 1
	for (int i = 0; i < N; i++)
	{
		// integrate 2*z / (z² + 1)
		// parametrization 1: z(t) = t -2i, t in [-2, 2]
		// Do ( 2*(t -2i) / ((t - 2i)² +1) ) * dt
		/*
		cMul(tmp_z1, z1, cval2);
		cMul(tmp_z2, z1, z1);
		mpfi_add_ui(tmp_z2->real, tmp_z2->real, 1);
		cDiv(z1_sum_tmp, tmp_z1, tmp_z2);
		cMul(z1_sum_tmp, z1_sum_tmp, dz1);
		cAdd(z1_sum, z1_sum, z1_sum_tmp);
		cAdd(z1, z1, dz1);
		*/
		cDiv(z1_sum_tmp, cval1, z1);
		cMul(z1_sum_tmp, z1_sum_tmp, dz1);
		cAdd(z1_sum, z1_sum, z1_sum_tmp);
		cAdd(z1, z1, dz1);
	}

	mpfr_t real_left, real_right, imaginary_left, imaginary_right;
	mpfr_init(real_left);
	mpfr_init(imaginary_left);
	mpfi_get_left(real_left, z1_sum->real);
	mpfi_get_left(imaginary_left, z1_sum->imag);
	mpfr_printf("z1 real: %.15Rf \n", real_left);
	mpfr_printf("z1 imaginary: %.15Rf \n", imaginary_left);


	// Set starting positions and values for integral 2
	mpfi_set(zi, left);
	mpfi_add(first, zi, dz);
	mpfi_put(zi, first);

	mpfi_set(z1->imag, left);
	mpfi_set_ui(z1->real, 2);
	mpfi_set_ui(dz1->real, 0);
	mpfi_set(dz1->imag, dz);
	mpfi_set_d(z2_sum->real, 0.0);
	mpfi_set_d(z2_sum->imag, 0.0);
	
	// Loop for integral 2
	for (int i = 0; i < N; i++)
	{
		// integrate 2*z / (z² + 1)
		// parametrization 2: z(t) = 2 + ti, t in [-2, 2]
		// Do ( 2*(2 + ti) / ((2 + ti)² +1) ) * dz
		/*
		cMul(tmp_z1, z1, cval2);
		cMul(tmp_z2, z1, z1);
		mpfi_add_ui(tmp_z2->real, tmp_z2->real, 1);
		cDiv(z2_sum_tmp, tmp_z1, tmp_z2);
		cMul(z2_sum_tmp, z2_sum_tmp, dz1);
		cAdd(z2_sum, z2_sum, z2_sum_tmp);
		cAdd(z1, z1, dz1);
		*/
		cDiv(z2_sum_tmp, cval1, z1);
		cMul(z2_sum_tmp, z2_sum_tmp, dz1);
		cAdd(z2_sum, z2_sum, z2_sum_tmp);
		cAdd(z1, z1, dz1);
	}
	
	
	mpfi_get_left(real_left, z2_sum->real);
	mpfi_get_left(imaginary_left, z2_sum->imag);
	mpfr_printf("z2 real: %.15Rf \n", real_left);
	mpfr_printf("z2 imaginary: %.15Rf \n", imaginary_left);
	
	// Set starting values for integral 3
	mpfi_set(zi, right);
	mpfi_sub(first, zi, dz);
	mpfi_put(first, zi);
	
	mpfi_set_ui(z1->imag, 2);
	mpfi_set(z1->real, right);
	mpfi_set(dz1->real, dz);
	mpfi_set_ui(dz1->imag, 0);
	cChange_sign(dz1, 1);
	mpfi_set_d(z3_sum->real, 0.0);
	mpfi_set_d(z3_sum->imag, 0.0);
	
	
	// Loop for integral 3
	for (int i = 0; i < N; i++)
	{
		// Integrate 2*z / (z² + 1)
		// Parametrization 3: t +2i, t in [2, -2]
		// Do (2*(t + 2i) / ((t + 2i)² + 1) dz
		/*
		cMul(tmp_z1, z1, cval2);
		cMul(tmp_z2, z1, z1);
		mpfi_add_ui(tmp_z2->real, tmp_z2->real, 1);
		cDiv(z3_sum_tmp, tmp_z1, tmp_z2);
		cMul(z3_sum_tmp, z3_sum_tmp, dz1);
		cAdd(z3_sum, z3_sum, z3_sum_tmp);
		cAdd(z1, z1, dz1);
		*/
		cDiv(z3_sum_tmp, cval1, z1);
		cMul(z3_sum_tmp, z3_sum_tmp, dz1);
		cAdd(z3_sum, z3_sum, z3_sum_tmp);
		cAdd(z1, z1, dz1);
	}

	mpfi_get_left(real_left, z3_sum->real);
	mpfi_get_left(imaginary_left, z3_sum->imag);
	mpfr_printf("z3 real: %.15Rf \n", real_left);
	mpfr_printf("z3 imaginary: %.15Rf \n", imaginary_left);

	
	// Set starting positions and values for integral 4
	mpfi_set(zi, right);
	mpfi_sub(first, zi, dz);
	mpfi_put(first, zi);
	
	mpfi_set(z1->imag, right);
	mpfi_set_si(z1->real, -2);
	mpfi_set_ui(dz1->real, 0);
	mpfi_set(dz1->imag, dz);
	cChange_sign(dz1, 2);
	mpfi_set_d(z4_sum->real, 0.0);
	mpfi_set_d(z4_sum->imag, 0.0);	
	
	
	// Loop for integral 4
	for (int i = 0; i < N; i++)
	{
		// integrate 2*z / (z² + 1)
		// parametrization 2: z(t) = 2 + ti, t in [2, -2]
		// Do ( 2*(2 + ti) / ((2 + ti)² +1) ) * dz
		/*
		cMul(tmp_z1, z1, cval2);
		cMul(tmp_z2, z1, z1);
		mpfi_add_ui(tmp_z2->real, tmp_z2->real, 1);
		cDiv(z4_sum_tmp, tmp_z1, tmp_z2);
		cMul(z4_sum_tmp, z4_sum_tmp, dz1);
		cAdd(z4_sum, z4_sum, z4_sum_tmp);
		cSub(z1, z1, dz1);
		*/
		
		cDiv(z4_sum_tmp, cval1, z1);
		cMul(z4_sum_tmp, z4_sum_tmp, dz1);
		cAdd(z4_sum, z4_sum, z4_sum_tmp);
		cAdd(z1, z1, dz1);
	}

	mpfi_get_left(real_left, z4_sum->real);
	mpfi_get_left(imaginary_left, z4_sum->imag);
	mpfr_printf("z4 real: %.15Rf \n", real_left);
	mpfr_printf("z4 imaginary: %.15Rf \n", imaginary_left);
	
	cAdd(sum_total, z1_sum, z2_sum);
	cAdd(sum_total, sum_total, z3_sum);
	cAdd(sum_total, sum_total, z4_sum);
	
	cPrint(sum_total);
	mpfr_init(real_right);
	mpfr_init(imaginary_right);
	mpfi_get_left(real_left, sum_total->real);
	mpfi_get_right(real_right, sum_total->real);
	mpfi_get_left(imaginary_left, sum_total->imag);
	mpfi_get_right(imaginary_right, sum_total->imag);
	mpfr_printf("total real left: %.15Rf \n", real_left);
	mpfr_printf("total real right: %.15Rf \n", real_right);
	mpfr_printf("total imaginary left: %.15Rf \n", imaginary_left);
	mpfr_printf("total imaginary right: %.15Rf \n", imaginary_right);
	
	
	cClear(sum_total);
	cClear(z1);
	cClear(dz1);
	cClear(z1_sum);
	cClear(z1_sum_tmp);
	cClear(tmp_z1);
	cClear(tmp_z2);
	cClear(cval1);
	cClear(cval2);
	cClear(z2_sum);
	cClear(z2_sum_tmp);
	cClear(z3_sum);
	cClear(z3_sum_tmp);
	cClear(z4_sum);
	cClear(z4_sum_tmp);
		

	mpfi_clear(left);
	mpfi_clear(right);
	mpfi_clear(val1);
	mpfi_clear(val2);
	mpfi_clear(ilength);
	mpfi_clear(valN);
	mpfi_clear(dz);
	mpfi_clear(zi);	
	mpfi_clear(first);
	
	mpfr_clear(real_left);
	mpfr_clear(real_right);
	mpfr_clear(imaginary_left);
	mpfr_clear(imaginary_right);
	
	return 1;
	
}
