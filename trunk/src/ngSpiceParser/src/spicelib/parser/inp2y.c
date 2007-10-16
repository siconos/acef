/**********
Copyright 1992 Regents of the University of California.  All rights
reserved.
Author: 1992 Charles Hough
**********/

#include "ngspice.h"
#include <stdio.h>
#include "ifsim.h"
#include "inpdefs.h"
#include "inpmacs.h"
#include "fteext.h"
#include "inp.h"

void
INP2Y(ckt,tab,current)
    void *ckt;
    INPtables *tab;
    card *current;

{
/* parse a txl card */
/* Yxxxx node1 node2 name */

int mytype; /* the type to determine txl */
int mytype2; /* the type to determine cpl */
int type;   /* the type the model says it is */
char *line; /* the part of the current line left to parse */
char *name; /* the resistor's name */
char *buf; /* temporary buffer for parsing */
char *model;    /* the name of the resistor's model */
char *nname1;   /* the first node's name */
char *nname2;   /* the second node's name */
char rname1[10], rname2[10], rname3[10];
char cname1[10], cname2[10], cname3[10], cname4[10];
char *internal1, *internal2;
char *ground1, *ground2;
void *node1; /* the first node's node pointer */
void *node2; /* the second node's node pointer */
void *gnode1, *gnode2, *inode1, *inode2;
int error;      /* error code temporary */
int error1;     /* secondary error code temporary */
INPmodel *thismodel;    /* pointer to model structure describing our model */
void *mdfast;    /* pointer to the actual model */
void *fast;  /* pointer to the actual instance */
void *mdfast2, *mdfast3, *mdfast4, *mdfast5, *mdfast6;
void *fast2, *fast3, *fast4, *fast5, *fast6;
IFuid uid;      /* uid for default model */
GENinstance *txl;
IFvalue ptemp;  /* a value structure to package into */
double lval=0, rval=0, cval=0, lenval=0;
int lenvalgiven = 0;

    mytype = INPtypelook("TransLine");
    mytype2 = INPtypelook("CplLines");

    if(mytype < 0 ) {
        LITERR("Device type TransLine not supported by this binary\n")
        return;
    }
    line = current->line;
    INPgetTok(&line,&name,1);
    INPinsert(&name,tab);
    INPgetNetTok(&line,&nname1,1);
    INPtermInsert(ckt,&nname1,tab,&node1);
    INPgetNetTok(&line,&ground1,1);
    INPtermInsert(ckt,&ground1,tab,&gnode1);
    INPgetNetTok(&line,&nname2,1);
    INPtermInsert(ckt,&nname2,tab,&node2);
    INPgetNetTok(&line,&ground2,1);
    INPtermInsert(ckt,&ground2,tab,&gnode2);
	
    INPgetTok(&line,&model,1);
	if (strcmp(model, "len") == 0) {
		lenval = INPevaluate(&line,&error1,1);
		lenvalgiven = 1;
		INPgetTok(&line,&model,1);
	}
	if(*model) { /* token isn't null */
		INPinsert(&model,tab);
		thismodel = (INPmodel *)NULL;
		current->error = INPgetMod(ckt,model,&thismodel,tab);
		if(thismodel != NULL) {
			if (thismodel->INPmodType == mytype2) {
				INP2P(ckt,tab,current);
				return;
			}
			else if (mytype != thismodel->INPmodType) { 
				LITERR("incorrect model type")
				return;
			}
			line = thismodel->INPmodLine->line;
			INPgetTok(&line,&buf,1);  /* throw out .model */
			INPgetTok(&line,&buf,1);  /* throw out model name */
			INPgetTok(&line,&buf,1);  /* throw out txl */
			INPgetTok(&line,&buf,1);
			while (*line != '\0') {
				if (*buf == 'R' || *buf == 'r') {
					INPgetTok(&line,&buf,1);
					rval = atof(buf);
				}
				if ((strcmp(buf,"L") == 0)  || (strcmp(buf,"l") == 0)) {
					INPgetTok(&line,&buf,1);
					lval = atof(buf);
				}
				if ((strcmp(buf,"C") == 0)  || (strcmp(buf,"c") == 0)) {
					INPgetTok(&line,&buf,1);
					cval = atof(buf);
				}
				if (lenvalgiven == 0) {
					if (strcmp(buf,"length")== 0) {
						INPgetTok(&line,&buf,1);
						lenval = atof(buf);
					}
				}
				INPgetTok(&line,&buf,1);
			}
			if (lenval && rval && lval && rval/lval > 1.6e10) {
				/* use 3-pi model */
				rval = 3.0 / (rval * lenval);
				cval = cval * lenval / 6.0;

				type = INPtypelook("Resistor");

				/* resistor between node1 and internal1 */
				internal1 = (char *) malloc (10);
				strcpy(internal1, "txlnd1");   
				INPtermInsert(ckt, &internal1, tab, &inode1);
				if(!tab->defRmod) {
					/* create default R model */
					IFnewUid(ckt,&uid,(IFuid)NULL,"R",UID_MODEL,(void**)NULL);
					IFC(newModel, (ckt,type,&(tab->defRmod),uid))
				}
				mdfast = tab->defRmod;
				strcpy(rname1, "txlres1");
				IFC(newInstance,(ckt,mdfast,&fast,rname1))
				IFC(bindNode,(ckt,fast,1,node1))
				IFC(bindNode,(ckt,fast,2,inode1))
				ptemp.rValue = rval;
				GCA(INPpName,("resistance",&ptemp,ckt,type,fast))

				/* resistor between internal1 and internal2 */
				internal2 = (char *) malloc (10);
				strcpy(internal2, "txlnd2");
				INPtermInsert(ckt, &internal2, tab, &inode2);
				strcpy(rname2, "txlres2");
				mdfast2 = tab->defRmod;
				IFC(newInstance,(ckt,mdfast2,&fast2,rname2))
				IFC(bindNode,(ckt,fast2,1,inode1))
				IFC(bindNode,(ckt,fast2,2,inode2))
				ptemp.rValue = rval;
				GCA(INPpName,("resistance",&ptemp,ckt,type,fast2))

				/* resistor between internal2 and node2 */
				strcpy(rname3, "txlres3");
				mdfast3 = tab->defRmod;
				IFC(newInstance,(ckt,mdfast3,&fast3,rname3))
				IFC(bindNode,(ckt,fast3,1,inode2))
				IFC(bindNode,(ckt,fast3,2,node2))
				ptemp.rValue = rval;
				GCA(INPpName,("resistance",&ptemp,ckt,type,fast3))

				/* capacitor on node1 */
				type = INPtypelook("Capacitor");
				if(!tab->defCmod) {
					IFnewUid(ckt,&uid,(IFuid)NULL,"C",UID_MODEL,(void**)NULL);
					IFC(newModel,(ckt,type,&(tab->defCmod),uid))
				}
				mdfast4 = tab->defCmod;
				strcpy(cname1, "txlcap1");
				IFC(newInstance,(ckt,mdfast4,&fast4,cname1))
				IFC(bindNode,(ckt,fast4,1,node1))
				IFC(bindNode,(ckt,fast4,2,gnode1))
				ptemp.rValue = cval;
				GCA(INPpName,("capacitance",&ptemp,ckt,type,fast4))

				/* capacitor on internal1 */
				strcpy(cname2, "txlcap2");
				mdfast4 = tab->defCmod;
				IFC(newInstance,(ckt,mdfast4,&fast4,cname2))
				IFC(bindNode,(ckt,fast4,1,inode1))
				IFC(bindNode,(ckt,fast4,2,gnode1))
				ptemp.rValue = cval * 2;
				GCA(INPpName,("capacitance",&ptemp,ckt,type,fast4))

				/* capacitor on internal2 */
				strcpy(cname3, "txlcap3");
				mdfast5 = tab->defCmod;
				IFC(newInstance,(ckt,mdfast5,&fast5,cname3))
				IFC(bindNode,(ckt,fast5,1,inode2))
				IFC(bindNode,(ckt,fast5,2,gnode1))
				ptemp.rValue = cval * 2;
				GCA(INPpName,("capacitance",&ptemp,ckt,type,fast5))

				/* capacitor on node2 */
				strcpy(cname4, "txlcap4");
				mdfast6 = tab->defCmod;
				IFC(newInstance,(ckt,mdfast6,&fast6,cname4))
				IFC(bindNode,(ckt,fast6,1,node2))
				IFC(bindNode,(ckt,fast6,2,gnode1))
				ptemp.rValue = cval;
				GCA(INPpName,("capacitance",&ptemp,ckt,type,fast6))
				return;

			}

			/* use regular txl model */
			mdfast = thismodel->INPmodfast;
			type = thismodel->INPmodType;
		} else {
			type = mytype;
			if(!tab->defYmod) {
				/* create default Y model */
				IFnewUid(ckt,&uid,(IFuid)NULL,"Y",UID_MODEL,(void**)NULL);
				IFC(newModel, (ckt,type,&(tab->defYmod),uid))
			}
			mdfast = tab->defYmod;
    	}
    	IFC(newInstance,(ckt,mdfast,&fast,name))
	} else  {
		LITERR("model name is not found")
		return;
	}
	
	if (error1 == 0 && lenvalgiven) {
		ptemp.rValue = lenval;
		GCA(INPpName,("length",&ptemp,ckt,type,fast))
	}

    IFC(bindNode,(ckt,fast,1,node1))
    IFC(bindNode,(ckt,fast,2,node2))

	txl = (GENinstance *)fast;

    return;
}
