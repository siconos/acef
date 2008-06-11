/**********
Copyright 1991 Regents of the University of California.  All rights reserved.
Author:	1987 Kartikeya Mayaram, U. C. Berkeley CAD Group
**********/

#include "ngspice.h"
#include "accuracy.h"

/*
 * this function computes the bernoulli function 
 * f(x) = x / ( exp (x) - 1 ) and its derivatives. the function
 * f(-x) alongwith its derivative is also computed.
 */

/*
 * input    delta-psi 
 * outputs  f(x), df(x)/dx, f(-x), and df(-x)/dx
 */

void bernoulli (double x, double *pfx, double *pDfxDx, double *pfMx, 
                double *pDfMxDx, BOOLEAN derivAlso)
{

}
