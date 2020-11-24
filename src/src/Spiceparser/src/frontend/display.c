/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
$Id: display.c,v 1.13 2005/05/31 00:12:18 sjborley Exp $
**********/


#include <ngspice.h>
#include <graph.h>
#include <ftedev.h>
#include <fteinput.h>
#include <cpdefs.h>     /* for VT_STRING */
#include <ftedefs.h>        /* for mylog() */

#include "display.h"
#include "variable.h"
#include "error.h"

/* static declarations */
static void gen_DatatoScreen(GRAPH *graph, double x, double y, int *screenx, int *screeny);
static int gen_Input(REQUEST *request, RESPONSE *response);
static int nop(void);
static int nodev(void);

#ifndef X_DISPLAY_MISSING
#include "plotting/x11.h"
#endif

#ifdef HAS_WINDOWS	/* Graphic-IO under MS Windows */
#include "wdisp/windisp.h"
#include "wdisp/winprint.h"
#endif

#include "plotting/plot5.h"
#include "postsc.h"
#include "hpgl.h"


DISPDEVICE device[] = {

    {"error", 0, 0, 0, 0, 0, 0, nop, nop,
    nop, nop,
    nop, nop, nop, nop, nop,
    nop, nop, nop,
    nop, nop, nop, gen_Input,
    (void*)nop,},


};

DISPDEVICE *dispdev = device + NUMELEMS(device) - 1;

#define XtNumber(arr)       (sizeof(arr) / sizeof(arr[0]))


DISPDEVICE *FindDev(char *name)
{
    int i;

    for (i=0; i < XtNumber(device); i++) {
      if (!strcmp(name, device[i].name)) {
        return(&device[i]);
      }
    }
    sprintf(ErrorMessage, "Can't find device %s.", name);
    internalerror(ErrorMessage);
    return(&device[0]);

}

void
DevInit(void)
{
#ifndef X_DISPLAY_MISSING
    char buf[128];   /* va: used with NOT X_DISPLAY_MISSING only */
#endif /* X_DISPLAY_MISSING */

/* note: do better determination */

/*
    dumb tradition that got passed on from gi_interface
    to do compile time determination
*/

    dispdev = NULL;

#ifndef X_DISPLAY_MISSING
    /* determine display type */
    if (getenv("DISPLAY") || cp_getvar("display", VT_STRING, buf)) {
       dispdev = FindDev("X11");
    }
#endif


#ifdef HAS_WINDOWS
	 if (!dispdev) {
      dispdev = FindDev("Windows");
    }
#endif

    if (!dispdev) {
	externalerror(
	 "no graphics interface; please check compiling instructions");
	dispdev = FindDev("error");
    } else if ((*(dispdev->Init))()) {
      fprintf(cp_err,
        "Warning: can't initialize display device for graphics.\n");
      dispdev = FindDev("error");
    }

}

/* NewViewport is responsible for filling in graph->viewport */
int
NewViewport(GRAPH *pgraph)
{

    return (*(dispdev->NewViewport))(pgraph);

}

void DevClose(void)
{

    (*(dispdev->Close))();

}

void DevClear(void)
{

    (*(dispdev->Clear))();

}

void DrawLine(int x1, int y1, int x2, int y2)
{
    (*(dispdev->DrawLine))(x1, y1, x2, y2);

}

void Arc(int x0, int y0, int radius, double theta1, double theta2)
{

    (*(dispdev->Arc))(x0, y0, radius, theta1, theta2);

}

void Text(char *text, int x, int y)
{

    (*(dispdev->Text))(text, x, y);

}

void DefineColor(int colorid, double red, double green, double blue)
{

    (*(dispdev->DefineColor))(colorid, red, green, blue);

}

void DefineLinestyle(int linestyleid, int mask)
{

    (*(dispdev->DefineLinestyle))(linestyleid, mask);

}

void SetLinestyle(int linestyleid)
{

    (*(dispdev->SetLinestyle))(linestyleid);

}

void SetColor(int colorid)
{

    (*(dispdev->SetColor))(colorid);

}

void Update(void)
{

    if (dispdev)
	    (*(dispdev->Update))();

}

/* note: screen coordinates are relative to window
    so need to add viewport offsets */
static void
gen_DatatoScreen(GRAPH *graph, double x, double y, int *screenx, int *screeny)
{

    double low, high;

    /* note: may want to cache datawindowsize/viewportsize */ /* done */

    /* note: think this out---Is 1 part of the viewport? Do we handle
        this correctly? */

    /* have to handle several types of grids */

    /* note: we can't compensate for X's demented y-coordinate system here
        since the grid routines use DrawLine w/o calling this routine */
    if ((graph->grid.gridtype == GRID_LOGLOG) ||
            (graph->grid.gridtype == GRID_YLOG)) {
      low = mylog10(graph->datawindow.ymin);
      high = mylog10(graph->datawindow.ymax);
      *screeny = (mylog10(y) - low) / (high - low) * graph->viewport.height
	  + 0.5 + graph->viewportyoff;
    } else {
      *screeny = ((y - graph->datawindow.ymin) / graph->aspectratioy)
            + 0.5 + graph->viewportyoff;
    }

    if ((graph->grid.gridtype == GRID_LOGLOG) ||
            (graph->grid.gridtype == GRID_XLOG)) {
      low = mylog10(graph->datawindow.xmin);
      high = mylog10(graph->datawindow.xmax);
      *screenx = (mylog10(x) - low) / (high - low) * graph->viewport.width
            + 0.5 + graph ->viewportxoff;
    } else {
      *screenx = (x - graph->datawindow.xmin) / graph->aspectratiox
            + 0.5 + graph ->viewportxoff;
    }

}

void DatatoScreen(GRAPH *graph, double x, double y, int *screenx, int *screeny)
{

    (*(dispdev->DatatoScreen))(graph, x, y, screenx, screeny);

}

void Input(REQUEST *request, RESPONSE *response)
{

    (*(dispdev->Input))(request, response);

}

static int
gen_Input(REQUEST *request, RESPONSE *response)
{

    switch (request->option) {
      case char_option:
        response->reply.ch = inchar(request->fp);
        response->option = request->option;
        break;
      default:
        /* just ignore, since we don't want a million error messages */
	if (response)
	    response->option = error_option;
        break;
    }
    return 0;
}

/* no operation, do nothing */
static int nop(void)
{
    return(1);  /* so NewViewport will fail */
}

static int nodev(void)
{

    sprintf(ErrorMessage,
        "This operation is not defined for display type %s.",
        dispdev->name);
    internalerror(ErrorMessage);
    return(1);

}

void SaveText(GRAPH *graph, char *text, int x, int y)
{

    struct _keyed *keyed;

    keyed = (struct _keyed *) tmalloc(sizeof(struct _keyed));

    if (!graph->keyed) {
      graph->keyed = keyed;
    } else {
      keyed->next = graph->keyed;
      graph->keyed = keyed;
    }

    keyed->text = tmalloc(strlen(text) + 1);
    strcpy(keyed->text, text);

    keyed->x = x;
    keyed->y = y;

    keyed->colorindex = graph->currentcolor;

}

/* if given name of a hardcopy device, finds it and switches devices
   if given NULL, switches back */
int DevSwitch(char *devname)
{

    static DISPDEVICE *lastdev = NULL;

    if (devname != NULL) {
      if (lastdev != NULL) {
        internalerror("DevSwitch w/o changing back");
        return (1);
      }
      lastdev = dispdev;
      dispdev = FindDev(devname);
      if (!strcmp(dispdev->name, "error")) {
        internalerror("no hardcopy device");
        dispdev = lastdev;  /* undo */
        lastdev = NULL;
        return (1);
      }
      (*(dispdev->Init))();
    } else {
      (*(dispdev->Close))();
      dispdev = lastdev;
      lastdev = NULL;
    }
    return(0);

}
