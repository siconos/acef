/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
**********/

#include "ngspice.h"
#include "const.h"
#include "ifsim.h"
#include "compdefs.h"
#include "sperror.h"


int
COMPmParam(int param, IFvalue *value, GENmodel *inModel)
{
    COMPmodel *model = (COMPmodel *)inModel;
    return(E_BADPARM);
}
