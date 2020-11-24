/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
**********/


#include "ngspice.h"
#include "cktdefs.h"
#include "complex.h"
#include "sperror.h"
#include "compdefs.h"



int
COMPpzLoad(GENmodel *inModel, CKTcircuit *ckt, SPcomplex *s)
        /* actually load the current resistance value into the 
         * sparse matrix previously provided 
         */
{
    COMPmodel *model = (COMPmodel *)inModel;
    COMPinstance *here;
    double m;

    return(OK);
}
