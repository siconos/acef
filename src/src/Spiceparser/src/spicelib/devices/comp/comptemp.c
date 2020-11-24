/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
**********/

#include "ngspice.h"
#include "cktdefs.h"
#include "compdefs.h"
#include "sperror.h"

int
COMPtemp(GENmodel *inModel, CKTcircuit *ckt)
        /* perform the temperature update to the resistors
         * calculate the conductance as a function of the
         * given nominal and current temperatures - the
         * resistance given in the struct is the nominal
         * temperature resistance
         */
{
    COMPmodel *model =  (COMPmodel *)inModel;
    COMPinstance *here;
    double factor;
    double difference;


    return(OK);
}
