/* 
Function integrating f(x) = cos(x)/(1+x^2),
x is real.
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
#include "math.h"

void int_cosx_div_1_plus_sqrx(mpfi_t sum, double start, double end, int Nr)
{
	// Declaring variables and constants
	mpfi_t cosx, divisor, val1, valN, x, dx, x_tmp, interval, fx;
	mpfr_t diam;
	double length;
	// Initializing variables and constants
	mpfi_init(cosx);
	mpfi_init(divisor);
	mpfi_init(val1);
	mpfi_init(valN);
	mpfi_init(x);
	mpfi_init(dx);
	mpfi_init(x_tmp);
	mpfi_init(interval);
	mpfi_init(fx);
	mpfr_init(diam);

	// Assigning values	
	length = end - start; // Interval length
	mpfi_set_d(interval, length);
	mpfi_set_si(valN, Nr); // Number of steps
	mpfi_set_d(x, start); // Start
	mpfi_set_d(x_tmp, start); // Temp start pos to calc first interval
	mpfi_set_ui(val1, 1);
	mpfi_set_ui(sum, 0);
	
	// Creating interval and step size
	mpfi_div(dx, interval, valN);
	mpfi_add(x_tmp, x_tmp, dx);
	mpfi_union(x, x, x_tmp);
	
	for (int i = 0; i < Nr; i++)
	{
		// Calculating integral with Riemann sum.
		mpfi_cos(cosx, x);
		mpfi_sqr(x_tmp, x);
		mpfi_add(divisor, val1, x_tmp);
		mpfi_div(fx, cosx, divisor);
		mpfi_diam_abs(diam, x);
		mpfi_set_fr(dx, diam);
		mpfi_mul(fx, fx, dx);
		mpfi_add(sum, sum, fx);
		mpfi_add(x, x, dx);
	}
	
	// Getting ends of sum interval and placing in mpfr structs to print.
	mpfr_t left, right;
	mpfr_init(left);
	mpfr_init(right);
	mpfi_get_left(left, sum);
	mpfi_get_right(right, sum);
	printf("The integral of cos(x)/(1+x^2) from -1 to 1 with %d points is\n", Nr);
	mpfr_printf("= [ %.10Rf, ", left ); mpfr_printf("%.10Rf ] \n", right );
	
	// Clearing variables.
	mpfr_clear(left);
	mpfr_clear(right);
	mpfr_clear(diam);
	mpfi_clear(cosx);
	mpfi_clear(divisor);
	mpfi_clear(val1);
	mpfi_clear(valN);
	mpfi_clear(x);
	mpfi_clear(dx);
	mpfi_clear(x_tmp);
	mpfi_clear(interval);
	mpfi_clear(fx);
}

