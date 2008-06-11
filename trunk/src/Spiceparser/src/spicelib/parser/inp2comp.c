/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 1988 Thomas L. Quarles
Modified: Paolo Nenzi 2000
Remarks:  This code is based on a version written by Serban Popescu which
          accepted an optional parameter ac. I have adapted his code
	  to conform to INP standard. (PN)
**********/

#include "ngspice.h"
#include <stdio.h>
#include "ifsim.h"
#include "inpdefs.h"
#include "inpmacs.h"
#include "fteext.h"
#include "inp.h"

/* undefine to add tracing to this file */
/*#define TRACE*/

void INP2COMP(void *ckt, INPtables * tab, card * current)
{
/* parse a resistor card */
/* Rname <node> <node> [<val>][<mname>][w=<val>][l=<val>][ac=<val>] */

    int mytype;			/* the type we determine resistors are */
    int type = 0;		/* the type the model says it is */
    char *line;			/* the part of the current line left to parse */
    char *saveline;		/* ... just in case we need to go back... */
    char *name;			/* the resistor's name */
    char *model;		/* the name of the resistor's model */
    char *nname1;		/* the first node's name */
    char *nname2;		/* the second node's name */
    char *nname3;		/* the third node's name */
    void *node1;		/* the first node's node pointer */
    void *node2;		/* the second node's node pointer */
    void *node3;		/* the second node's node pointer */
    double val1;			/* temp to held resistance */
    double val2;			/* temp to held resistance */
    double val3;			/* temp to held resistance */
    int error;			/* error code temporary */
    int error1;			/* secondary error code temporary */
    INPmodel *thismodel;	/* pointer to model structure describing our model */
    void *mdfast = NULL;	/* pointer to the actual model */
    void *fast;			/* pointer to the actual instance */
    IFvalue ptemp;		/* a value structure to package resistance into */
    int waslead;		/* flag to indicate that funny unlabeled number was found */
    double leadval;		/* actual value of unlabeled number */
    IFuid uid;			/* uid for default model */

    char *s,*t;/* Temporary buffer and pointer for translation */
    int i;
	
#ifdef TRACE
    printf("In INP2COMP()\n");
    printf("  Current line: %s\n",current->line);
#endif

    mytype = INPtypelook("Comparator");
    if (mytype < 0) {
	LITERR("Device type Comparator not supported by this binary\n");
	return;
    }
    line = current->line;
    INPgetTok(&line, &name, 1);			/* Cname */
    INPinsert(&name, tab);
    INPgetNetTok(&line, &nname1, 1);		/* <node> */
    INPtermInsert(ckt, &nname1, tab, &node1);
    INPgetNetTok(&line, &nname2, 1);		/* <node> */
    INPtermInsert(ckt, &nname2, tab, &node2);
    INPgetNetTok(&line, &nname3, 1);		/* <node> */
    INPtermInsert(ckt, &nname3, tab, &node3);
    /*val1 = INPevaluate(&line, &error1, 1);	
    val2 = INPevaluate(&line, &error1, 1);	
    val3 = INPevaluate(&line, &error1, 1);*/	
    /* either not a number -> model, or
     * follows a number, so must be a model name
     * -> MUST be a model name (or null)
     */
    
    saveline = line;		/* save then old pointer */
    

#ifdef TRACE
    printf("(Translated) Res line: %s\n",current->line);
#endif

    INPgetTok(&line, &model, 1);

    if (*model) {
	/* token isn't null */
	if (INPlookMod(model)) {
	    /* If this is a valid model connect it */
	    INPinsert(&model, tab);
	    thismodel = (INPmodel *) NULL;
	    current->error = INPgetMod(ckt, model, &thismodel, tab);
	    if (thismodel != NULL) {
		if (mytype != thismodel->INPmodType) {
		    LITERR("incorrect model type for comparator");
		    return;
		}
		mdfast = thismodel->INPmodfast;
		type = thismodel->INPmodType;
	    }
	} else {
	    tfree(model);
	    /* It is not a model */
	    line = saveline;	/* go back */
	    type = mytype;
	    if (!tab->defCOMPmod) {	/* create default R model */
		IFnewUid(ckt, &uid, (IFuid) NULL, "COMP", UID_MODEL,
			 (void **) NULL);
		IFC(newModel, (ckt, type, &(tab->defCOMPmod), uid));
	    }
	    mdfast = tab->defCOMPmod;
	}
	IFC(newInstance, (ckt, mdfast, &fast, name));
    } else {
	tfree(model);
	/* The token is null and a default model will be created */
	type = mytype;
	if (!tab->defCOMPmod) {
	    /* create default R model */
	    IFnewUid(ckt, &uid, (IFuid) NULL, "COMP", UID_MODEL,
		     (void **) NULL);
	    IFC(newModel, (ckt, type, &(tab->defCOMPmod), uid));
	}
	IFC(newInstance, (ckt, tab->defCOMPmod, &fast, name));
    }

    if (error1 == 0) {		/* got a resistance above */
	ptemp.rValue = val1;
	GCA(INPpName, ("comparator", &ptemp, ckt, type, fast))
    }

    IFC(bindNode, (ckt, fast, 1, node1));
    IFC(bindNode, (ckt, fast, 2, node2));
    IFC(bindNode, (ckt, fast, 3, node3));
    PARSECALL((&line, ckt, type, fast, &leadval, &waslead, tab));
    if (waslead) {
	ptemp.rValue = leadval;
	GCA(INPpName, ("comparator", &ptemp, ckt, type, fast));
    }

    return;
}
