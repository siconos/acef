/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 1985 Wayne A. Christopher, U. C. Berkeley CAD Group
**********/

/*
 * Routines to do complex mathematical functions. These routines require
 * the -lm libraries. We sacrifice a lot of space to be able
 * to avoid having to do a seperate call for every vector element,
 * but it pays off in time savings.  These routines should never
 * allow FPE's to happen.
 *
 * Complex functions are called as follows:
 *  cx_something(data, type, length, &newlength, &newtype),
 *  and return a char * that is cast to complex or double.
 */

#include <ngspice.h>
#include <cpdefs.h>
#include <ftedefs.h>
#include <dvec.h>

#include "cmath.h"
#include "cmath2.h"

static double *
d_tan(double *dd, int length)
{
    double *d;
    int i;

    d = alloc_d(length);
    for (i = 0; i < length; i++) {
	rcheck(cos(degtorad(dd[i])) != 0, "tan");
	d[i] = sin(degtorad(dd[i])) / cos(degtorad(dd[i]));
    }
    return d;
}

static complex *
c_tan(complex *cc, int length)
{
    complex *c;
    int i;

    c = alloc_c(length);
    for (i = 0; i < length; i++) {
	double u, v;

	rcheck(cos(degtorad(realpart(&cc[i]))) *
	       cosh(degtorad(imagpart(&cc[i]))), "tan");
	rcheck(sin(degtorad(realpart(&cc[i]))) *
	       sinh(degtorad(imagpart(&cc[i]))), "tan");
	u = degtorad(realpart(&cc[i]));
	v = degtorad(imagpart(&cc[i]));
        /* The Lattice C compiler won't take multi-line macros, and
         * CMS won't take >80 column lines....
         */
#define xx1 sin(u) * cosh(v)
#define xx2 cos(u) * sinh(v)
#define xx3 cos(u) * cosh(v)
#define xx4 sin(u) * sinh(v)
        cdiv(xx1, xx2, xx3, xx4, realpart(&c[i]), imagpart(&c[i]));
    }
    return c;
}

void *
cx_tan(void *data, short int type, int length, int *newlength, short int *newtype, ...)
{
    *newlength = length;
    if (type == VF_REAL) {
        *newtype = VF_REAL;
	return (void *) d_tan((double *) data, length);
    } else {
        *newtype = VF_COMPLEX;
        return (void *) c_tan((complex *) data, length);
    }
}



void *
cx_atan(void *data, short int type, int length, int *newlength, short int *newtype, ...)
{
    double *d;

    d = alloc_d(length);
    *newtype = VF_REAL;
    *newlength = length;
    if (type == VF_COMPLEX) {
	complex *cc = (complex *) data;
	int i;

        for (i = 0; i < length; i++)
            d[i] = radtodeg(atan(realpart(&cc[i])));
    } else {
	double *dd = (double *) data;
	int i;

        for (i = 0; i < length; i++)
            d[i] = radtodeg(atan(dd[i]));
    }
    return ((void *) d);
}


double
cx_max_local(void *data, short int type, int length)
{
    double largest = 0.0;

    if (type == VF_COMPLEX) {
	complex *cc = (complex *) data;
	int i;

        for (i = 0; i < length; i++)
            if (cmag(&cc[i]) > largest)
                largest = cmag(&cc[i]);
    } else {
	double *dd = (double *) data;
	int i;

        for (i = 0; i < length; i++)
            if (FTEcabs(dd[i]) > largest)
                largest = FTEcabs(dd[i]);
    }
    return largest;
}

/* Normalize the data so that the magnitude of the greatest value is 1. */

void *
cx_norm(void *data, short int type, int length, int *newlength, short int *newtype, ...)
{
    double largest = 0.0;

    largest = cx_max_local(data, type, length);
    if (largest == 0.0) {
	fprintf(cp_err, "Error: can't normalize a 0 vector\n");
	return (NULL);
    }

    *newlength = length;
    if (type == VF_COMPLEX) {
	complex *c;
	complex *cc = (complex *) data;
	int i;

        c = alloc_c(length);
        *newtype = VF_COMPLEX;

        for (i = 0; i < length; i++) {
            realpart(&c[i]) = realpart(&cc[i]) / largest;
            imagpart(&c[i]) = imagpart(&cc[i]) / largest;
        }
        return ((void *) c);
    } else {
	double *d;
	double *dd = (double *) data;
	int i;

        d = alloc_d(length);
        *newtype = VF_REAL;

        for (i = 0; i < length; i++)
            d[i] = dd[i] / largest;
        return ((void *) d);
    }
}

void *
cx_uminus(void *data, short int type, int length, int *newlength, short int *newtype, ...)
{
    *newlength = length;
    if (type == VF_COMPLEX) {
	complex *c;
	complex *cc = (complex *) data;
	int i;

        c = alloc_c(length);
        *newtype = VF_COMPLEX;
        for (i = 0; i < length; i++) {
            realpart(&c[i]) = - realpart(&cc[i]);
            imagpart(&c[i]) = - imagpart(&cc[i]);
        }
        return ((void *) c);
    } else {
	double *d;
	double *dd = (double *) data;
	int i;

        d = alloc_d(length);
        *newtype = VF_REAL;
        for (i = 0; i < length; i++)
            d[i] = - dd[i];
        return ((void *) d);
    }
}

void *
cx_rnd(void *data, short int type, int length, int *newlength, short int *newtype, ...)
{
    *newlength = length;
    if (type == VF_COMPLEX) {
	complex *c;
	complex *cc = (complex *) data;
	int i;

        c = alloc_c(length);
        *newtype = VF_COMPLEX;
        for (i = 0; i < length; i++) {
	    int j, k;

            j = floor(realpart(&cc[i]));
            k = floor(imagpart(&cc[i]));
            realpart(&c[i]) = j ? random() % j : 0;
            imagpart(&c[i]) = k ? random() % k : 0;
        }
        return ((void *) c);
    } else {
	double *d;
	double *dd = (double *) data;
	int i;

        d = alloc_d(length);
        *newtype = VF_REAL;
        for (i = 0; i < length; i++) {
	    int j;

            j = floor(dd[i]);
            d[i] = j ? random() % j : 0;
        }
        return ((void *) d);
    }
}

/* Compute the mean of a vector. */

void *
cx_mean(void *data, short int type, int length, int *newlength, short int *newtype, ...)
{
    *newlength = 1;
    rcheck(length > 0, "mean");
    if (type == VF_REAL) {
	double *d;
	double *dd = (double *) data;
	int i;

        d = alloc_d(1);
        *newtype = VF_REAL;
        for (i = 0; i < length; i++)
            *d += dd[i];
        *d /= length;
        return ((void *) d);
    } else {
	complex *c;
	complex *cc = (complex *) data;
	int i;

        c = alloc_c(1);
        *newtype = VF_COMPLEX;
        for (i = 0; i < length; i++) {
            realpart(c) += realpart(cc + i);
            imagpart(c) += imagpart(cc + i);
        }
        realpart(c) /= length;
        imagpart(c) /= length;
        return ((void *) c);
    }
}


void *
cx_length(void *data, short int type, int length, int *newlength, short int *newtype, ...)
{
    double *d;

    *newlength = 1;
    *newtype = VF_REAL;
    d = alloc_d(1);
    *d = length;
    return ((void *) d);
}

/* Return a vector from 0 to the magnitude of the argument. Length of the
 * argument is irrelevent.
 */


void *
cx_vector(void *data, short int type, int length, int *newlength, short int *newtype, ...)
{
    complex *cc = (complex *) data;
    double *dd = (double *) data;
    int i, len;
    double *d;

    if (type == VF_REAL)
        len = FTEcabs(*dd);
    else
        len = cmag(cc);
    if (len == 0)
        len = 1;
    d = alloc_d(len);
    *newlength = len;
    *newtype = VF_REAL;
    for (i = 0; i < len; i++)
        d[i] = i;
    return ((void *) d);
}

/* Create a vector of the given length composed of all ones. */


void *
cx_unitvec(void *data, short int type, int length, int *newlength, short int *newtype, ...)
{
    complex *cc = (complex *) data;
    double *dd = (double *) data;
    int i, len;
    double *d;

    if (type == VF_REAL)
        len = FTEcabs(*dd);
    else
        len = cmag(cc);
    if (len == 0)
        len = 1;
    d = alloc_d(len);
    *newlength = len;
    *newtype = VF_REAL;
    for (i = 0; i < len; i++)
        d[i] = 1;
    return ((void *) d);
}

/* Calling methods for these functions are:
 *  cx_something(data1, data2, datatype1, datatype2, length)
 *
 * The length of the two data vectors is always the same, and is the length
 * of the result. The result type is complex iff one of the args is
 * complex.
 */

void *
cx_plus(void *data1, void *data2, short int datatype1, short int datatype2, int length, ...)
{
    double *dd1 = (double *) data1;
    double *dd2 = (double *) data2;
    double *d;
    complex *cc1 = (complex *) data1;
    complex *cc2 = (complex *) data2;
    complex *c, c1, c2;
    int i;

    if ((datatype1 == VF_REAL) && (datatype2 == VF_REAL)) {
        d = alloc_d(length);
        for (i = 0; i < length; i++)
            d[i] = dd1[i] + dd2[i];
        return ((void *) d);
    } else {
        c = alloc_c(length);
        for (i = 0; i < length; i++) {
            if (datatype1 == VF_REAL) {
                realpart(&c1) = dd1[i];
                imagpart(&c1) = 0.0;
            } else {
                realpart(&c1) = realpart(&cc1[i]);
                imagpart(&c1) = imagpart(&cc1[i]);
            }
            if (datatype2 == VF_REAL) {
                realpart(&c2) = dd2[i];
                imagpart(&c2) = 0.0;
            } else {
                realpart(&c2) = realpart(&cc2[i]);
                imagpart(&c2) = imagpart(&cc2[i]);
            }
            realpart(&c[i]) = realpart(&c1) + realpart(&c2);
            imagpart(&c[i]) = imagpart(&c1) + imagpart(&c2);
        }
        return ((void *) c);
    }
}

void *
cx_minus(void *data1, void *data2, short int datatype1, short int datatype2, int length, ...)
{
    double *dd1 = (double *) data1;
    double *dd2 = (double *) data2;
    double *d;
    complex *cc1 = (complex *) data1;
    complex *cc2 = (complex *) data2;
    complex *c, c1, c2;
    int i;

    if ((datatype1 == VF_REAL) && (datatype2 == VF_REAL)) {
        d = alloc_d(length);
        for (i = 0; i < length; i++)
            d[i] = dd1[i] - dd2[i];
        return ((void *) d);
    } else {
        c = alloc_c(length);
        for (i = 0; i < length; i++) {
            if (datatype1 == VF_REAL) {
                realpart(&c1) = dd1[i];
                imagpart(&c1) = 0.0;
            } else {
                realpart(&c1) = realpart(&cc1[i]);
                imagpart(&c1) = imagpart(&cc1[i]);
            }
            if (datatype2 == VF_REAL) {
                realpart(&c2) = dd2[i];
                imagpart(&c2) = 0.0;
            } else {
                realpart(&c2) = realpart(&cc2[i]);
                imagpart(&c2) = imagpart(&cc2[i]);
            }
            realpart(&c[i]) = realpart(&c1) - realpart(&c2);
            imagpart(&c[i]) = imagpart(&c1) - imagpart(&c2);
        }
        return ((void *) c);
    }
}

void *
cx_times(void *data1, void *data2, short int datatype1, short int datatype2, int length, ...)
{
    double *dd1 = (double *) data1;
    double *dd2 = (double *) data2;
    double *d;
    complex *cc1 = (complex *) data1;
    complex *cc2 = (complex *) data2;
    complex *c, c1, c2;
    int i;

    if ((datatype1 == VF_REAL) && (datatype2 == VF_REAL)) {
        d = alloc_d(length);
        for (i = 0; i < length; i++)
            d[i] = dd1[i] * dd2[i];
        return ((void *) d);
    } else {
        c = alloc_c(length);
        for (i = 0; i < length; i++) {
            if (datatype1 == VF_REAL) {
                realpart(&c1) = dd1[i];
                imagpart(&c1) = 0.0;
            } else {
                realpart(&c1) = realpart(&cc1[i]);
                imagpart(&c1) = imagpart(&cc1[i]);
            }
            if (datatype2 == VF_REAL) {
                realpart(&c2) = dd2[i];
                imagpart(&c2) = 0.0;
            } else {
                realpart(&c2) = realpart(&cc2[i]);
                imagpart(&c2) = imagpart(&cc2[i]);
            }
            realpart(&c[i]) = realpart(&c1) * realpart(&c2)
                - imagpart(&c1) * imagpart(&c2);
            imagpart(&c[i]) = imagpart(&c1) * realpart(&c2)
                + realpart(&c1) * imagpart(&c2);
        }
        return ((void *) c);
    }
}

void *
cx_mod(void *data1, void *data2, short int datatype1, short int datatype2, int length, ...)
{
    double *dd1 = (double *) data1;
    double *dd2 = (double *) data2;
    double *d;
    complex *cc1 = (complex *) data1;
    complex *cc2 = (complex *) data2;
    complex *c, c1, c2;
    int i, r1, r2, i1, i2, r3, i3;

    if ((datatype1 == VF_REAL) && (datatype2 == VF_REAL)) {
        d = alloc_d(length);
        for (i = 0; i < length; i++) {
            r1 = floor(FTEcabs(dd1[i]));
            rcheck(r1 > 0, "mod");
            r2 = floor(FTEcabs(dd2[i]));
            rcheck(r2 > 0, "mod");
            r3 = r1 % r2;
            d[i] = (double) r3;
        }
        return ((void *) d);
    } else {
        c = alloc_c(length);
        for (i = 0; i < length; i++) {
            if (datatype1 == VF_REAL) {
                realpart(&c1) = dd1[i];
                imagpart(&c1) = 0.0;
            } else {
                realpart(&c1) = realpart(&cc1[i]);
                imagpart(&c1) = imagpart(&cc1[i]);
            }
            if (datatype2 == VF_REAL) {
                realpart(&c2) = dd2[i];
                imagpart(&c2) = 0.0;
            } else {
                realpart(&c2) = realpart(&cc2[i]);
                imagpart(&c2) = imagpart(&cc2[i]);
            }
            r1 = floor(FTEcabs(realpart(&c1)));
            rcheck(r1 > 0, "mod");
            r2 = floor(FTEcabs(realpart(&c2)));
            rcheck(r2 > 0, "mod");
            i1 = floor(FTEcabs(imagpart(&c1)));
            rcheck(i1 > 0, "mod");
            i2 = floor(FTEcabs(imagpart(&c2)));
            rcheck(i2 > 0, "mod");
            r3 = r1 % r2;
            i3 = i1 % i2;
            realpart(&c[i]) = (double) r3;
            imagpart(&c[i]) = (double) i3;
        }
        return ((void *) c);
    }
}


/* Routoure JM : Compute the max of a vector. */

void *
cx_max(void *data, short int type, int length, int *newlength, short int *newtype, ...)
{
    *newlength = 1;
    /* test if length >0 et affiche un message d'erreur */
    rcheck(length > 0, "mean");
    if (type == VF_REAL) {
      double largest=0.0;
      double *d;
      double *dd = (double *) data;
      int i;
      
      d = alloc_d(1);
      *newtype = VF_REAL;
      largest=dd[0];
      for (i = 1; i < length; i++)
        if (dd[i]>largest) largest=dd[i];
      *d=largest;
      return ((void *) d);
    } else { 
      double largest_real=0.0;
      double largest_complex=0.0;
      complex *c;
      complex *cc = (complex *) data;
      int i;
      
      c = alloc_c(1);
      *newtype = VF_COMPLEX;
      largest_real=realpart(cc);
      largest_complex=imagpart(cc);
      for (i = 0; i < length; i++) {
        if (realpart(cc + i)>largest_real) largest_real=realpart(cc + i);
        if (imagpart(cc + i)>largest_complex) largest_complex=imagpart(cc + i);
        }
        realpart(c) = largest_real;
        imagpart(c) = largest_complex;
        return ((void *) c);
    }
}
/* Routoure JM : Compute the min of a vector. */

void *
cx_min(void *data, short int type, int length, int *newlength, short int *newtype, ...)
{
    *newlength = 1;
    /* test if length >0 et affiche un message d'erreur */
    rcheck(length > 0, "mean");
    if (type == VF_REAL) {
      double smallest;
      double *d;
      double *dd = (double *) data;
      int i;
      
      d = alloc_d(1);
      *newtype = VF_REAL;
      smallest=dd[0];
      for (i = 1; i < length; i++)
        if (dd[i]<smallest) smallest=dd[i];
      *d=smallest;
      return ((void *) d);
    } else { 
      double smallest_real;
      double smallest_complex;
      complex *c;
      complex *cc = (complex *) data;
      int i;
      
      c = alloc_c(1);
      *newtype = VF_COMPLEX;
      smallest_real=realpart(cc);
      smallest_complex=imagpart(cc);
      for (i = 1; i < length; i++) {
        if (realpart(cc + i)<smallest_real) smallest_real=realpart(cc + i);
        if (imagpart(cc + i)<smallest_complex) smallest_complex=imagpart(cc + i);
        }
        realpart(c) = smallest_real;
        imagpart(c) = smallest_complex;
        return ((void *) c);
    }
}


/* Routoure JM : Compute the differential  of a vector. */

void *
cx_d(void *data, short int type, int length, int *newlength, short int *newtype, ...)
{
    *newlength = length;
    /* test if length >0 et affiche un message d'erreur */
    rcheck(length > 0, "deriv");
    if (type == VF_REAL) {
      double *d;
      double *dd = (double *) data;
      int i;
      
      d = alloc_d(length);
      *newtype = VF_REAL;
      d[0]=dd[1]-dd[0];
      d[length-1]=dd[length-1]-dd[length-2];
      for (i = 1; i < length-1; i++)
        d[i]=dd[i+1]-dd[i-1];

      return ((void *) d);
    } else { 

      complex *c;
      complex *cc = (complex *) data;
      int i;
      
      c = alloc_c(length);
      *newtype = VF_COMPLEX;
      realpart(c)=realpart(cc+1)-realpart(cc);
      imagpart(c)=imagpart(cc+1)-imagpart(cc);
      realpart(c+length-1)=realpart(cc+length-1)-realpart(cc+length-2);
      imagpart(c+length-1)=imagpart(cc+length-1)-imagpart(cc+length-2);
      
      
      for (i = 1; i < (length-1); i++) {
        realpart(c+i)=realpart(cc+i+1)-realpart(cc+i-1);
        imagpart(c+i)=imagpart(cc+i+1)-imagpart(cc+i-1);

        }
        return ((void *) c);
    }
}
