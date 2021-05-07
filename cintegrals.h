/*
.h file defining functions calculating closed curve integrals on the complex plane.
Author: @JM.
*/
#include <stdio.h>
#include "mpfi.h"
#include "mpfi_io.h"
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include <mpfr.h>
//#include "complex_func.h"

// Integrals on simple polynomials

// This is 1/z
// @inputs: result struct for sum, starting pos, ending pos, number of calcs
void cint_one_div_z(cmpfi, double, double, int);

// This is 2z/(z^2 + 1)
// @inputs: result struct for sum, starting pos, ending pos, number of calcs
void cint_twoz_div_sqrz_plus1(cmpfi, double, double, int);

// This is cos(x)/(1+x^2), x in R ========= This is MPFI ONLY  ==========
// @inputs: result struct for sum, starting pos, ending pos, number of calcs
void int_cosx_div_1_plus_sqrx(mpfi_t, double, double, int);
