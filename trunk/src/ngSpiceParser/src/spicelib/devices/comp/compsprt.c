/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
**********/

/* Pretty print the sensitivity info for all 
 * the resistors in the circuit.
 */

#include "ngspice.h"
#include "cktdefs.h"
#include "compdefs.h"
#include "sperror.h"


void
COMPsPrint(GENmodel *inModel, CKTcircuit *ckt)
{
    COMPmodel *model = (COMPmodel *)inModel;
    COMPinstance *here;
    printf("COMPARATORS-----------------\n");

    /*  loop through all the comparator models */
    for( ; model != NULL; model = model->COMPnextModel ) {

        printf("Model name:%s\n",model->COMPmodName);

        /* loop through all the comparator of the model */
        for (here = model->COMPinstances; here != NULL ;
                here=here->COMPnextInstance) {
	    
	    if (here->COMPowner != ARCHme) continue;

            printf("    Instance name:%s\n",here->COMPname);
            printf("      Positive, negative, out nodes name: %s, %s, %s\n",
            CKTnodName(ckt,here->COMPposNode),CKTnodName(ckt,here->COMPnegNode),CKTnodName(ckt,here->COMPoutNode));
            printf("      Positive, negative nodes: %d, %d, %d\n",
		   here->COMPposNode,here->COMPnegNode,here->COMPoutNode);
	            
	    printf("  Vmin: %f \n",here->COMPVmoins);
	    printf("  Vmax: %f \n",here->COMPVplus);
	    printf("  Epsilon: %f \n",here->COMPEpsilon);

        }
    }
}
