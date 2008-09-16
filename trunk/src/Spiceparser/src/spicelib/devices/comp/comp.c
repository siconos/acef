/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
**********/

#include "ngspice.h"
#include "compdefs.h"
#include "devdefs.h"

#include "ifsim.h"

IFparm COMPpTable[] = { /* parameters */ 
 IOPP(  "vmin", 	COMP_VMIN, 	     IF_REAL,    "Vmin"),
 IOPP(  "vmax", 	COMP_VMAX, 	     IF_REAL,    "Vmax"),
 IOPP(  "vepsilon", 	COMP_VEPSILON, 	     IF_REAL,    "Vepsilon"),
 IOPP(  "voffset", 	COMP_VOFFSET, 	     IF_REAL,    "Voffset")
};

IFparm COMPmPTable[] = {
};

char *COMPnames[] = {
    "Z+",
    "Z-",
    "Zout",
    "Zoffset"
};

int	COMPnSize = NUMELEMS(COMPnames);
int	COMPpTSize = NUMELEMS(COMPpTable);
int	COMPmPTSize = NUMELEMS(COMPmPTable);
int	COMPiSize = sizeof(COMPinstance);
int	COMPmSize = sizeof(COMPmodel);
