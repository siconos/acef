/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 2008 O. Bonnefon
**********/

#ifndef COMP
#define COMP

#include "ifsim.h"
#include "cktdefs.h"
#include "gendefs.h"
#include "complex.h"
#include "noisedef.h"

        /* definitions used to describe resistors */


/* information used to describe a single instance */

typedef struct sCOMPinstance {
    struct sCOMPmodel *COMPmodPtr;            /* backpointer to model */
    struct sCOMPinstance *COMPnextInstance;   /* pointer to next instance of 
                                             * current model*/

    IFuid COMPname;      /* pointer to character string naming this instance */
    int COMPowner;       /* number of owner process */
    int COMPstate;       /* not used but needed for sructure consistency */
    int COMPposNode;     /* number of positive node of comparator */
    int COMPnegNode;     /* number of negative node of comparator */

  int COMPoutNode;     /* number of out node of comparator */
  double COMPVplus;
  double COMPVmoins;
  double COMPEpsilon;
  double COMPOffset;
  int    COMPsenParmNo; 


} COMPinstance ;


/* per model data */

typedef struct sCOMPmodel {       /* model structure for a resistor */
    int COMPmodType; /* type index of this device type */
    struct sCOMPmodel *COMPnextModel; /* pointer to next possible model in 
                                     * linked list */
    COMPinstance * COMPinstances; /* pointer to list of instances that have this
                                 * model */
    IFuid COMPmodName;       /* pointer to character string naming this model */

} COMPmodel;

#define COMP_VMIN 1
#define COMP_VMAX 2
#define COMP_VEPSILON 3
#define COMP_VOFFSET 4

/* model questions */

#include "compext.h"

#endif /*COMP*/
