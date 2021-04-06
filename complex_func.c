/*
 Functions for operations on complex numbers with mpfi library.


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


// Initialize complext struct.
void cInit(cmpfi z)
{
	mpfi_init( z->real );
	mpfi_init( z->imag );
}

// Clear
void cClear(cmpfi z)
{
	mpfi_clear( z->real );
	mpfi_clear( z->imag ); 
}


// Fill struct
void cSet(cmpfi z, mpfi_t a, mpfi_t b)
{
	mpfi_set( z->real, a);
	mpfi_set( z->imag, b);
}


// Function for complex addition
void cAdd(cmpfi z, cmpfi w, cmpfi u)
{
	// complex format a+bi, c + di
	// Do a + c
	// Do bi + di
	mpfi_add( z->real, w->real, u->real);
	mpfi_add( z->imag, w->imag, u->imag);
}

// Function for complex subtraction

void cSub(cmpfi z, cmpfi w, cmpfi u)
{
	// complex format a+bi, c + di
	// Do a - c
	// Do bi - di
	mpfi_sub(z->real, w->real, u->real);
	mpfi_sub(z->imag, w->imag, u->imag);
	
	
}


// Function for complex multiplication
void cMul(cmpfi z, cmpfi w, cmpfi u)
{
	// Initializing temporary variables.
	mpfi_t temp_real, temp_imag;
	mpfi_init(temp_real);
	mpfi_init(temp_imag);

	// Real product
	mpfi_mul( z->real, w->real, u->real);
	mpfi_mul(temp_real, w->imag, u->imag);
	mpfi_sub(z->real, z->real, temp_real);
	
	// Imaginary product
	mpfi_mul(temp_imag, w->real, u->imag);
	mpfi_mul( z->imag, w->imag, u->real);
	mpfi_add( z->imag, z->imag, temp_imag);
	
	// Clearing temporary variables.
	mpfi_clear(temp_real);
	mpfi_clear(temp_imag);
	
	
	// complex format a+bi, c + di
	// Do a * c + bi * di
	// Do a * di + c * bi
}

// Function for complex division

void cDiv(cmpfi z, cmpfi w, cmpfi u)
{
	// Initializing temporary variables.
	mpfi_t temp_real, temp_imag, temp_divisor1, temp_divisor2;
	mpfi_init(temp_real);
	mpfi_init(temp_imag);
	mpfi_init(temp_divisor1);
	mpfi_init(temp_divisor2);
	
	// Real dividend
	mpfi_mul(z->real, w->real, u->real);
	mpfi_mul(temp_real, w->imag, u->imag);
	mpfi_add(z->real, z->real, temp_real);
	
	// Divisor
	mpfi_mul(temp_divisor1, u->real, u->real);
	mpfi_mul(temp_divisor2, u->imag, u->imag);	
	mpfi_add(temp_divisor1, temp_divisor1, temp_divisor2);
	
	// Imaginary dividend
	mpfi_mul(z->imag, w->imag, u->real);
	mpfi_mul(temp_imag, w->real, u->imag);
	mpfi_sub(z->imag, z->imag, temp_imag);
	
	// Placing result in appropriate struct
	mpfi_div(z->real, z->real, temp_divisor1);
	mpfi_div(z->imag, z->imag, temp_divisor1);
	
	// Clearing temporary variables.
	mpfi_clear(temp_real);
	mpfi_clear(temp_imag);
	mpfi_clear(temp_divisor1);
	mpfi_clear(temp_divisor2);
	// complex format a+bi, c + di
	// Do a * c + bi * di / (c * c + d * d)
	// Do bi * c - a * di / (c * c + d * d)	
}


///////////////////////////////// ADD TO .h FILE ////////////////////////////////

// Set unsigned int, a real part, b imaginary part
void cSet_ui(cmpfi z, unsigned int a, unsigned int b)
{
	mpfi_set_ui(z->real, a);
	mpfi_set_ui(z->imag, b);
}	

// Set signed int, a real part, b imaginary part
void cSet_si(cmpfi z, int a, int b)
{
	mpfi_set_si(z->real, a);
	mpfi_set_si(z->imag, b);
}

// Set double, a real part, b imaginary part
void cSet_d(cmpfi z, long double a, long double b)
{
	mpfi_set_d(z->real, a);
	mpfi_set_d(z->imag, b);
}

// Set real part
void cSet_real(cmpfi z, mpfi_t a)
{
	mpfi_set(z->real, a);
}

// Set imaginary part
void cSet_imag(cmpfi z, mpfi_t b)
{
	mpfi_set(z->imag, b);
}


// Print left and right bound of real and imaginary intervals.
void cPrint(cmpfi z)
{
	mpfr_t real_l, real_r, imag_l, imag_r;
	mpfr_init(real_l);
	mpfr_init(real_r);
	mpfr_init(imag_l);
	mpfr_init(imag_r);
	
	// Get sides and place int temp variables
	mpfi_get_left(real_l, z->real);
	mpfi_get_right(real_r, z->real);
	mpfi_get_left(imag_l, z->imag);
	mpfi_get_right(imag_r, z->imag);
	
	// Print format Re(z) = [ a, b]
	mpfr_printf("Re(z) = [ %.10Rf,", real_l); mpfr_printf(" %.10Rf ] \n", real_r);
	mpfr_printf("Im(z) = [ %.10Rf,", imag_l); mpfr_printf(" %.10Rf ] \n", imag_r);
	
	mpfr_clear(real_l);
	mpfr_clear(real_r);
	mpfr_clear(imag_l);
	mpfr_clear(imag_r);
}


// Create a union between w, and u and place it in z
void cUnion(cmpfi z, cmpfi w, cmpfi u)
{
	mpfi_union(z->real, w->real, u->real);
	mpfi_union(z->imag, w->imag, u->imag);
}

// Change sign

void cChange_sign(cmpfi z, int flag)
{
mpfi_t real_flip, imag_flip;

	if (flag == 1) // Change sign of real part, flag 1
	{
		mpfi_init(real_flip);
		mpfi_set_si(real_flip, 0);
		mpfi_sub(z->real, real_flip, z->real);
		mpfi_clear(real_flip);
	}
	else if (flag == 2) // Change sign of imaginary part, flag 2
	{
		mpfi_init(imag_flip);
		mpfi_set_si(imag_flip, 0);
		mpfi_sub(z->imag, imag_flip, z->imag);
		mpfi_clear(imag_flip);
	}
	else if (flag == 3) // Change signs of both, flag 3
	{
		mpfi_init(real_flip);
		mpfi_set_si(real_flip, 0);
		mpfi_sub(z->real, real_flip, z->real);
		mpfi_init(imag_flip);
		mpfi_set_si(imag_flip, 0);
		mpfi_sub(z->imag, imag_flip, z->imag);
		mpfi_clear(real_flip);
		mpfi_clear(imag_flip);
	}
	else // Wrong input error.
	{
		printf("function input (complex mpfi struct, int flag)\n");
		printf("flag == 1 flips sign of real part,\n");
		printf("flag == 2 flips sign of imaginary part,\n");
		printf("flag == 3 flips sign of both.");
		printf("Input error, sign flip function.");
		exit(-1);
	}
}


// Get diameter of intervals in w and place in z
void cDiam(cmpfi z, cmpfi w)
{
	mpfr_t temp_real, temp_imag;
	mpfr_init(temp_real);
	mpfr_init(temp_imag);
	mpfi_diam_abs(temp_real, w->real);
	mpfi_diam_abs(temp_imag, w->imag);
	mpfi_set_fr(z->real, temp_real);
	mpfi_set_fr(z->imag, temp_imag);
	mpfr_clear(temp_real);
	mpfr_clear(temp_imag);
}





























