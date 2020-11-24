/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
**********/
/*
 */

#include "ngspice.h"
#include "compdefs.h"
#include "sperror.h"


int
COMPdelete(GENmodel *inModel, IFuid name, GENinstance **inst)
{
    COMPmodel *model = (COMPmodel *)inModel;
    COMPinstance **fast = (COMPinstance **)inst;
    COMPinstance **prev = NULL;
    COMPinstance *here;

    for( ; model ; model = model->COMPnextModel) {
        prev = &(model->COMPinstances);
        for(here = *prev; here ; here = *prev) {
            if(here->COMPname == name || (fast && here==*fast) ) {
                *prev= here->COMPnextInstance;
                FREE(here);
                return(OK);
            }
            prev = &(here->COMPnextInstance);
        }
    }
    return(E_NODEV);
}
