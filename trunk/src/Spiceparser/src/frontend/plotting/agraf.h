/*************
 * Header file for agraf.c
 * 1999 E. Rouat
 ************/

#ifndef _AGRAF_H
#define _AGRAF_H

#include <dvec.h>
#include <bool.h>
#include <plot.h>

void ft_agraf(double *xlims, double *ylims, struct dvec *xscale,
	      struct plot *plot, struct dvec *vecs,
	      double xdel, double ydel, bool xlog, bool ylog, 
	      bool nointerp);

#endif
