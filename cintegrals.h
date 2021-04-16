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
void cint_one_div_z(cmpfi, double, double, int);

// This is 2z/(z^2 + 1)
void cint_twoz_div_sqrz_plus1(cmpfi, double, double, int);
