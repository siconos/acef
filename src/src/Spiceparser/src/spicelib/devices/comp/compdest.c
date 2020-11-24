/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
**********/
/*
 */

#include "ngspice.h"
#include "compdefs.h"


void
COMPdestroy(GENmodel **inModel)
{
    COMPmodel **model = (COMPmodel **)inModel;
    COMPinstance *here;
    COMPinstance *prev = NULL;
    COMPmodel *mod = *model;
    COMPmodel *oldmod = NULL;

    for( ; mod ; mod = mod->COMPnextModel) {
        if(oldmod) FREE(oldmod);
        oldmod = mod;
        prev = (COMPinstance *)NULL;
        for(here = mod->COMPinstances ; here ; here = here->COMPnextInstance) {
            if(prev) FREE(prev);
            prev = here;
        }
        if(prev) FREE(prev);
    }
    if(oldmod) FREE(oldmod);
    *model = NULL;
}
