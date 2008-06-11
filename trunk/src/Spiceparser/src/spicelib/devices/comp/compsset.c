/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
**********/

#include "ngspice.h"
#include "smpdefs.h"
#include "compdefs.h"
#include "sperror.h"
#include "cktdefs.h"


int
COMPsSetup(SENstruct *info, GENmodel *inModel)
        /* loop through all the devices and 
         * assign parameter #s to design parameters 
         */
{
    COMPmodel *model = (COMPmodel *)inModel;
    COMPinstance *here;

    /*  loop through all the resistor models */
    for( ; model != NULL; model = model->COMPnextModel ) {

        /* loop through all the instances of the model */
        for (here = model->COMPinstances; here != NULL ;
            here=here->COMPnextInstance) {
	    
	    if (here->COMPowner != ARCHme) continue;

            if(here->COMPsenParmNo){
                here->COMPsenParmNo = ++(info->SENparms);
            }
        }
    }
    return(OK);
}
