/* 
Function integrating f'(z)/f(z) where f(z) = z^2 + 1 in a counter-clockwise fashion around a square
of chosen limits. z is complex.
f(z) has known zeroes at +- i.

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
void cint_twoz_div_sqrz_plus1(cmpfi sum, double start, double end, int Nr)
{
	// Declaring internal variables and constants
	cmpfi z1, z2, z3, z4; // Intervals for the integrals
	cmpfi f1_tmp, f2_tmp, f3_tmp, f4_tmp; // Temporary function values
	cmpfi f1_prime, f2_prime, f3_prime, f4_prime; // Temp function derivatives 
	cmpfi z1_aux, z2_aux, z3_aux, z4_aux; // Temporary to create intervals.
	cmpfi z1_sum, z2_sum, z3_sum, z4_sum; // Sums for the integrals
	cmpfi z1_tmp_sum, z2_tmp_sum, z3_tmp_sum, z4_tmp_sum; // Temporary sums
	cmpfi val1, val2, dz1, dz3, idz2, idz4; // Constants
	mpfi_t step, N, interval_length;
	double length;
	// Initialize structs
	cInit(z1);
	cInit(z2);
	cInit(z3);
	cInit(z4);
	cInit(f1_tmp);
	cInit(f2_tmp);
	cInit(f3_tmp);
	cInit(f4_tmp);
	cInit(f1_prime);
	cInit(f2_prime);
	cInit(f3_prime);
	cInit(f4_prime);
	cInit(z1_aux);
	cInit(z2_aux);
	cInit(z3_aux);
	cInit(z4_aux);
	cInit(z1_sum);
	cInit(z2_sum);
	cInit(z3_sum);
	cInit(z4_sum);
	cInit(z1_tmp_sum);
	cInit(z2_tmp_sum);
	cInit(z3_tmp_sum);
	cInit(z4_tmp_sum);
	cInit(val1);
	cInit(val2);
	cInit(dz1);
	cInit(dz3);
	cInit(idz2);
	cInit(idz4);
	mpfi_init(step);
	mpfi_init(N);
	mpfi_init(interval_length);
	// Set constants 
	length = end - start;
	mpfi_set_ui(N, Nr);
	mpfi_set_d(interval_length, length);
	cSet_ui(val1, 1, 0);
	cSet_ui(val2, 2, 0);
	// Zeroing sums
	cSet_ui(z1_sum, 0, 0);
	cSet_ui(z2_sum, 0, 0);
	cSet_ui(z3_sum, 0, 0);
	cSet_ui(z4_sum, 0, 0);
	// Calculate and set step sizes
	mpfi_div(step, interval_length, N); // Calculating step size
	mpfi_set(dz1->real, step); // Assigning Re(pos_dz) = step size
	mpfi_set_ui(dz1->imag, 0); // Assigning Im(pos_dz) = 0
	mpfi_set(dz3->real, step); // Assigning Re(neg_dz) = step size
	mpfi_set_ui(dz3->imag, 0); // Assigning Im(neg_dz) = 0
	cChange_sign(dz3, 3); // Changing the sign of neg_dz
	mpfi_set(idz2->imag, step); // Assigning Im(pos_idz) = step size
	mpfi_set_ui(idz2->real, 0); // Assigning Re(pos_idz) = 0
	mpfi_set(idz4->imag, step); // Assigning Im(neg_idz) = step size
	mpfi_set_ui(idz4->real, 0); // Assigning Re(neg_idz) = 0
	cChange_sign(idz4, 3); // Chaning the sign of neg_idz
	// Set intervals
	cSet_d(z1, start, start);
	cSet_d(z2, end, start);
	cSet_d(z3, end, end);
	cSet_d(z4, start, end);
	// Take first steps and create a union between step 0 and step 1
	cSet_c(z1_aux, z1);
	cSet_c(z2_aux, z2);
	cSet_c(z3_aux, z3);
	cSet_c(z4_aux, z4);
	cAdd(z1_aux, z1_aux, dz1);
	cAdd(z2_aux, z2_aux, idz2);
	cAdd(z3_aux, z3_aux, dz3);
	cAdd(z4_aux, z4_aux, idz4);
	cUnion(z1, z1, z1_aux);
	cUnion(z2, z2, z2_aux);
	cUnion(z3, z3, z3_aux);
	cUnion(z4, z4, z4_aux);
	// Calculate integrals
	for (int i = 0; i < Nr; i++)
	{
		// Creating integrands for our four integrals
		// This is f'(z) = 2z
		cMul(f1_prime, val2, z1);
		cMul(f2_prime, val2, z2);
		cMul(f3_prime, val2, z3);
		cMul(f4_prime, val2, z4);
		
		// This is f(z) = z^2 + 1
		//cSqr(f1_tmp, z1);
		//cSqr(f2_tmp, z2);
		//cSqr(f3_tmp, z3);
		//cSqr(f4_tmp, z4);
		cMul(f1_tmp, z1, z1);
		cMul(f2_tmp, z2, z2);
		cMul(f3_tmp, z3, z3);
		cMul(f4_tmp, z4, z4);
		cAdd(f1_tmp, f1_tmp, val1);
		cAdd(f2_tmp, f2_tmp, val1);
		cAdd(f3_tmp, f3_tmp, val1);
		cAdd(f4_tmp, f4_tmp, val1);
		
		// This is f'(z)/f(z)
		cDiv(z1_tmp_sum, f1_prime, f1_tmp);
		cDiv(z2_tmp_sum, f2_prime, f2_tmp);
		cDiv(z3_tmp_sum, f3_prime, f3_tmp);
		cDiv(z4_tmp_sum, f4_prime, f4_tmp);
		
		// Multiplying integrands with corresponding step size
		cMul(z1_tmp_sum, z1_tmp_sum, dz1);
		cMul(z2_tmp_sum, z2_tmp_sum, idz2);
		cMul(z3_tmp_sum, z3_tmp_sum, dz3);
		cMul(z4_tmp_sum, z4_tmp_sum, idz4);
		
		// Adding result to corresponding sum
		cAdd(z1_sum, z1_sum, z1_tmp_sum);
		cAdd(z2_sum, z2_sum, z2_tmp_sum);
		cAdd(z3_sum, z3_sum, z3_tmp_sum);
		cAdd(z4_sum, z4_sum, z4_tmp_sum);
		
		// Advancing interval
		cAdd(z1, z1, dz1);
		cAdd(z2, z2, idz2);
		cAdd(z3, z3, dz3);
		cAdd(z4, z4, idz4);
	}

	// Summing partial results
	cAdd(sum, z1_sum, z2_sum);
	cAdd(sum, sum, z3_sum);
	cAdd(sum, sum, z4_sum);

	// Clear structs
	cClear(z1);
	cClear(z2);
	cClear(z3);
	cClear(z4);
	cClear(f1_tmp);
	cClear(f2_tmp);
	cClear(f3_tmp);
	cClear(f4_tmp);
	cClear(f1_prime);
	cClear(f2_prime);
	cClear(f3_prime);
	cClear(f4_prime);
	cClear(z1_aux);
	cClear(z2_aux);
	cClear(z3_aux);
	cClear(z4_aux);
	cClear(z1_sum);
	cClear(z2_sum);
	cClear(z3_sum);
	cClear(z4_sum);
	cClear(z1_tmp_sum);
	cClear(z2_tmp_sum);
	cClear(z3_tmp_sum);
	cClear(z4_tmp_sum);
	cClear(val1);
	cClear(val2);
	cClear(dz1);
	cClear(dz3);
	cClear(idz2);
	cClear(idz4);
	mpfi_clear(step);
	mpfi_clear(N);
	mpfi_clear(interval_length);
}
