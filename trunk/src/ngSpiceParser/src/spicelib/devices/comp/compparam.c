/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
**********/

#include "ngspice.h"
#include "const.h"
#include "ifsim.h"
#include "compdefs.h"
#include "sperror.h"


int
COMPparam(int param, IFvalue *value, GENinstance *inst, IFvalue *select)
{
    COMPinstance *here = (COMPinstance *)inst;
    switch(param) {
        case COMP_VMIN:
            here->COMPVmoins = value->rValue ;
            break;
	case COMP_VMAX:
            here->COMPVplus = value->rValue;
            break;   
        case COMP_VEPSILON:
            here->COMPEpsilon = value->rValue;
            break;
        default:
            return(E_BADPARM);
    }
    return(OK);
}
