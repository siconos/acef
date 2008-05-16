/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 1985 Thomas L. Quarles
Modified: Apr 2000 - Paolo Nenzi
**********/

#include "ngspice.h"
#include "cktdefs.h"
#include "compdefs.h"
#include "sperror.h"


/* actually load the current resistance value into the sparse matrix
 * previously provided */
int
COMPload(GENmodel *inModel, CKTcircuit *ckt)
{
    COMPmodel *model = (COMPmodel *)inModel;
    double m;
    double difference;
    double factor;

    return(OK);
}


/* actually load the current resistance value into the sparse matrix
 * previously provided */
int
COMPacload(GENmodel *inModel, CKTcircuit *ckt)
{
    COMPmodel *model = (COMPmodel *)inModel;
    double m;
    double difference;
    double factor;

    return(OK);
}
