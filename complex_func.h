/*
 .h file defining functions and structs.

Author: @JM.
*/
#include <stdio.h>
#include "mpfi.h"
#include "mpfi_io.h"
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include <mpfr.h>





/*
	Creating a struct of mpfi structs handling complex operations
	
*/

typedef struct complex_mpfi{
  mpfi_t real;
  mpfi_t imag;
}complex_mpfi;

typedef complex_mpfi cmpfi[1];
typedef complex_mpfi *cmpfi_ptr;
typedef const complex_mpfi cmpfi_srcptr;

void cInit(cmpfi);
void cClear(cmpfi);

// Assign operations
void cSet(cmpfi, mpfi_t, mpfi_t);

void cSet_ui(cmpfi, unsigned int, unsigned int);

void cSet_si(cmpfi, int, int);

void cSet_d(cmpfi, long double, long double);

void cSet_real(cmpfi, mpfi_t);

void cSet_imag(cmpfi, mpfi_t);

// Print function(s)
void cPrint(cmpfi);


// operations on complex numbers
void cAdd(cmpfi, cmpfi, cmpfi);

void cSub(cmpfi, cmpfi, cmpfi);

void cMul(cmpfi, cmpfi, cmpfi);

void cDiv(cmpfi, cmpfi, cmpfi);

void cChange_sign(cmpfi, int);

// Set operations
void cUnion(cmpfi, cmpfi, cmpfi);
void cDiam(cmpfi, cmpfi);


