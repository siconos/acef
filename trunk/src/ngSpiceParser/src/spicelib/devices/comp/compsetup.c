/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 1985 Thomas L. Quarles
Modified: Apr 2000 Paolo Nenzi
**********/

#include "ngspice.h"
#include "smpdefs.h"
#include "compdefs.h"
#include "sperror.h"


int 
COMPsetup(SMPmatrix *matrix, GENmodel *inModel, CKTcircuit*ckt, int *state)
        /* load the resistor structure with those pointers needed later 
         * for fast matrix loading 
         */
{

    return(OK);
}
