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

// Initialization and clear
// @input: cmpfi struct
void cInit(cmpfi);
void cClear(cmpfi);

// Assign operations
// @input: cmpfi struct, mpfi struct for real value, mpfi_struct for imag value
void cSet(cmpfi, mpfi_t, mpfi_t);
// @input: cmpfi struct, unsigned int for real value, unsigned int for imag value
void cSet_ui(cmpfi, unsigned int, unsigned int);
// @input: cmpfi struct, int for real value, int for imag value
void cSet_si(cmpfi, int, int);
// @input: cmpfi struct, double for real value, double for imag value
void cSet_d(cmpfi, long double, long double);
// @input: cmpfi struct, mpfi struct for real value
void cSet_real(cmpfi, mpfi_t);
// @input: cmpfi struct, mpfi struct for imag value
void cSet_imag(cmpfi, mpfi_t);
// @input: cmpfi struct, cmpfi struct
void cSet_c(cmpfi, cmpfi);

// Print function(s)
// @input: cmpfi struct
void cPrint(cmpfi);


// operations on complex numbers
// @input: cmpfi struct for sum, cmpfi struct, cmpfi struct
void cAdd(cmpfi, cmpfi, cmpfi);
// @input: cmpfi struct for difference, cmpfi struct, cmpfi struct
void cSub(cmpfi, cmpfi, cmpfi);
// @inout: cmpfi struct for product, cmpfi struct, cmpfi struct
void cMul(cmpfi, cmpfi, cmpfi);
// @input: cmpfi struct for result, cmpfi struct
void cSqr(cmpfi, cmpfi);
// @input: cmpfi struct for ratio, cmpfi struct, cmpfi struct
void cDiv(cmpfi, cmpfi, cmpfi);
// @input: cmpfi struct, int flag for which sign to change.
void cChange_sign(cmpfi, int);

// Set operations
// @input: cmpfi struct for result, cmpfi struct, cmpfi struct
void cUnion(cmpfi, cmpfi, cmpfi);
// @input: cmpfi struct for result, cmpfi struct, cmpfi struct
void cDiam(cmpfi, cmpfi);


