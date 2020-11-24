/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: Jeffrey M. Hsu
Modified 1999 Emmanuel Rouat
**********/

#include <config.h>
#include "ngspice.h"

#ifndef X_DISPLAY_MISSING

#include "cpstd.h"
#include "hlpdefs.h"

void newtopic(Widget w, caddr_t client_data, caddr_t call_data), delete(Widget w, caddr_t client_data, caddr_t call_data), quit(Widget w, caddr_t client_data, caddr_t call_data); 

/* Create a new window... */
bool
hlp_xdisplay(topic *top)
{
    return (TRUE);

}

void
newtopic(Widget w, caddr_t client_data, caddr_t call_data)
{
        return;
}

void
delete(Widget w, caddr_t client_data, caddr_t call_data)
{
}

void
quit(Widget w, caddr_t client_data, caddr_t call_data)
{
}

void
hlp_xkillwin(topic *top)
{
    return;
}

/* rip out font changes and write at end of buffer */

#endif  /* X_DISPLAY_MISSING */
