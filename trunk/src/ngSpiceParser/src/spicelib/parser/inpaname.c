/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 1985 Thomas L. Quarles
**********/

    /*
     * INPaName()
     *
     *  Take a parameter by Name and ask for the specified value 
     * *dev is -1 if type unknown, otherwise, device type
     * **fast is a device, and will be set if possible.
     */

#include "ngspice.h"
#include <stdio.h>
#include "cpdefs.h"
#include "fteext.h"
#include "ifsim.h"
#include "iferrmsg.h"
#include "inp.h"

int
INPaName(char *parm, IFvalue * val, void *ckt, int *dev, char *devnam,
	 void **fast, IFsimulator * sim, int *dataType, IFvalue * selector)
		    /* the name of the parameter to set */
		    /* the parameter union containing the value to set */
		    /* the circuit this device is a member of */
		    /* the device type code to the device being parsed */
		    /* the name of the device */
		     /* direct pointer to device being parsed */
			    /* the simulator data structure */
		    /* the datatype of the returned value structure */
			/* data sub-selector for questions */
{
    int error;			/* int to store evaluate error return codes in */
    int i;

    /* find the instance - don't know about model, so use null there,
     * otherwise pass on as much info as we have about the device
     * (name, type, direct pointer) - the type and direct pointer
     * WILL be set on return unless error is not OK
     */
    error = (*(sim->findInstance)) (ckt, dev, fast, devnam, (void *) NULL,
				    (char *) NULL);
    if (error)
	return (error);

    /* now find the parameter - hunt through the parameter tables for
     * this device type and look for a name match of an 'ask'able
     * parameter.
     */
    for (i = 0; i < (*(*(sim->devices)[*dev]).numInstanceParms); i++) {
	if (strcmp(parm,
		   ((*(sim->devices)[*dev]).instanceParms[i].keyword)) == 0
	    && (((*(sim->devices)[*dev]).instanceParms[i].dataType) &
		IF_ASK)) {
	    /* found it, so we ask the question using the device info we got
	     * above and put the results in the IFvalue structure our caller
	     * gave us originally
	     */
	    error = (*(sim->askInstanceQuest)) (ckt, *fast,
						(*(sim->devices)[*dev]).
						instanceParms[i].id, val,
						selector);
	    if (dataType)
		*dataType =
		    (*(sim->devices)[*dev]).instanceParms[i].dataType;
	    return (error);
	}
    }
    return (E_BADPARM);
}
