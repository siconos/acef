/* Configuration file for ng-spice */
#include <config.h>

/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
**********/

#include "ngspice.h"


#define CONFIG


#include "devdefs.h"
#include "noisedef.h"
#include "suffix.h"


/* XXX Should be -1 ? There is always an extra null element at the end ? */
static char * specSigList[] = {
    "time"
};

static IFparm nodeParms[] = {
    IP( "nodeset",PARM_NS ,IF_REAL,"suggested initial voltage"),
    IP( "ic",PARM_IC ,IF_REAL,"initial voltage"),
    IP( "type",PARM_NODETYPE ,IF_INTEGER,"output type of equation")
};

IFsimulator SIMinfo = {
    "ngspice",			/* name */
    "Circuit level simulation program",	/* more about me */
    Spice_Version,		/* version */

    CKTinit,			/* newCircuit function */
    CKTdestroy,			/* deleteCircuit function */

    CKTnewNode,			/* newNode function */
    CKTground,			/* groundNode function */
    CKTbindNode,		/* bindNode function */
    CKTfndNode,			/* findNode function */
    CKTinst2Node,		/* instToNode function */
    CKTsetNodPm,		/* setNodeParm function */
    CKTaskNodQst,		/* askNodeQuest function */
    CKTdltNod,			/* deleteNode function */

    CKTcrtElt,			/* newInstance function */
    CKTparam,			/* setInstanceParm function */
    CKTask,			/* askInstanceQuest function */
    CKTfndDev,			/* findInstance funciton */
    CKTdltInst,			/* deleteInstance function */

    CKTmodCrt,			/* newModel function */
    CKTmodParam,		/* setModelParm function */
    CKTmodAsk,			/* askModelQuest function */
    CKTfndMod,			/* findModel function */
    CKTdltMod,			/* deleteModel function */

    CKTnewTask,			/* newTask function */
    CKTnewAnal,			/* newAnalysis function */
    CKTsetAnalPm,		/* setAnalysisParm function */
    CKTaskAnalQ,		/* askAnalysisQuest function */
    CKTfndAnal,			/* findAnalysis function */
    CKTfndTask,			/* findTask function */
    CKTdelTask,			/* deleteTask function */

    CKTdoJob,			/* doAnalyses function */
    CKTtrouble,			/* non-convergence message function */

    0,				/* Initialized in SIMinit() */
    NULL,			/* Initialized in SIMinit() */
    0,				/* Initialized in SIMinit() */
    NULL,			/* Initialized in SIMinit() */

    sizeof(nodeParms)/sizeof(IFparm),
    nodeParms,

    sizeof(specSigList)/sizeof(char *),
    specSigList,
};
