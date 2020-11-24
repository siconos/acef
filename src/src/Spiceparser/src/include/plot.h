#ifndef _PLOT_H
#define _PLOT_H

#include "wordlist.h"
#include "bool.h"
#include "dvec.h"

/* The information for a particular set of vectors that come from one
 * plot.  */
struct plot {
    char *pl_title;		/* The title card. */
    char *pl_date;		/* Date. */
    char *pl_name;		/* The plot name. */
    char *pl_typename;		/* Tran1, op2, etc. */
    struct dvec *pl_dvecs;	/* The data vectors in this plot. */
    struct dvec *pl_scale;	/* The "scale" for the rest of the vectors. */
    struct plot *pl_next;	/* List of plots. */
    wordlist *pl_commands;	/* Commands to execute for this plot. */
    struct variable *pl_env;	/* The 'environment' for this plot. */
    char *pl_ccom;		/* The ccom struct for this plot. */
    bool pl_written;		/* Some or all of the vecs have been saved. */
    int pl_ndims;		/* Number of dimensions */
} ;


#endif
