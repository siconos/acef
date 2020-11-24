/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
**********/

#include "ngspice.h"
#include "compdefs.h"
#include "sperror.h"


int
COMPmDelete(GENmodel **inModel, IFuid modname, GENmodel *kill)
{
    COMPmodel **model = (COMPmodel **)inModel;
    COMPmodel *modfast = (COMPmodel *)kill;
    COMPinstance *here;
    COMPinstance *prev = NULL;
    COMPmodel **oldmod;
    oldmod = model;
    for( ; *model ; model = &((*model)->COMPnextModel)) {
        if( (*model)->COMPmodName == modname || 
                (modfast && *model == modfast) ) goto delgot;
        oldmod = model;
    }
    return(E_NOMOD);

delgot:
    *oldmod = (*model)->COMPnextModel; /* cut deleted device out of list */
    for(here = (*model)->COMPinstances ; here ; here = here->COMPnextInstance) {
        if(prev) FREE(prev);
        prev = here;
    }
    if(prev) FREE(prev);
    FREE(*model);
    return(OK);

}
