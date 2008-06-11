/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
**********/


#include "ngspice.h"
#include "const.h"
#include "cktdefs.h"
#include "ifsim.h"
#include "compdefs.h"
#include "sperror.h"
#include "devdefs.h"


int 
COMPmodAsk(CKTcircuit *ckt, GENmodel *inModel, int which, IFvalue *value)
{

    COMPmodel *model = (COMPmodel *)inModel;
    return(E_BADPARM);

}

