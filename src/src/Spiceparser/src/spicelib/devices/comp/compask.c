/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
**********/

#include "ngspice.h"
#include "const.h"
#include "compdefs.h"
#include "ifsim.h"
#include "cktdefs.h"
#include "sperror.h"


/* TODO : there are "double" value compared with 0 (eg: vm == 0)
 *        Need to substitute this check with a suitable eps.
 *        PN 2003
 */ 

int
COMPask(CKTcircuit *ckt, GENinstance *inst, int which, IFvalue *value, 
       IFvalue *select)
{
    COMPinstance *fast = (COMPinstance *)inst;
    double vr;
    double vi;
    double sr;
    double si;
    double vm;
    static char *msg = "Current and power not available for ac analysis";
    
    switch(which) {
        case COMP_VMIN:
            value->rValue = fast->COMPVmoins;
            return(OK);
	case COMP_VMAX:
            value->rValue = fast->COMPVplus;
            return(OK);    
        case COMP_VEPSILON:
            value->rValue = fast->COMPEpsilon;
            return(OK);
        case COMP_VOFFSET:
            value->rValue = fast->COMPOffset;
            return(OK);
        default:
            return(E_BADPARM);
    }
    /* NOTREACHED */
}
