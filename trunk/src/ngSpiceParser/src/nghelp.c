/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 1986 Wayne A. Christopher, U. C. Berkeley CAD Group 
**********/

/*
 * The main routine for the help system in stand-alone mode.
 */

#include <config.h>
#include "ngspice.h"
#include "cpdefs.h"
#include "hlpdefs.h"

#include "frontend/variable.h"

#ifndef X_DISPLAY_MISSING
Widget toplevel;
#endif


FILE *cp_in, *cp_out, *cp_err;


/* dummy declaration so CP.a doesn't pull in lexical.o and other objects */
bool cp_interactive = FALSE;

char *hlp_filelist[] = { "ngspice", 0 };

extern void ivars();

#ifdef HAS_WINDOWS
bool oflag = FALSE;         /* Output over redefined I/O functions */
FILE *flogp;  /* hvogt 15.12.2001 */
#endif /* HAS_WINDOWS */

int
#ifdef HAS_WINDOWS
xmain(int ac, char **av)
#else
main(int ac, char **av)
#endif /* HAS_WINDOWS */
{
    wordlist *wl = NULL;

#ifndef X_DISPLAY_MISSING
    char *displayname;
    /* grrr, Xtk forced contortions */
    char *argv[2];
    int argc = 2;
    char buf[512];
#endif /* X_DISPLAY_MISSING */

    ivars( );

    cp_in = stdin;
    cp_out = stdout;
    cp_err = stderr;

#ifndef X_DISPLAY_MISSING

    if (cp_getvar("display", VT_STRING, buf)) {
      displayname = buf;
    } else if (!(displayname = getenv("DISPLAY"))) {
      fprintf(stderr, "Can't open X display.");
      goto out;
    }

    argv[0] = "nutmeg";
    argv[1] = displayname;
    /* initialize X toolkit */
    toplevel = XtInitialize("nutmeg", "Nutmeg", NULL, 0, &argc, argv);
    
out:
#endif /* X_DISPLAY_MISSING */

    if (ac > 1)
        wl = wl_build(av + 1);
    hlp_main(Help_Path, wl);

#ifndef X_DISPLAY_MISSING
    if (hlp_usex) {
	printf("Hit control-C when done.\n");		/* sigh */
	XtMainLoop();
    }
#endif /* X_DISPLAY_MISSING */

#ifdef HAS_WINDOWS
	/* Keep window open untill a key is pressed */
	printf("Press a key to quit\n");
	while( getchar() == EOF) {}
#endif /* HAS_WINDOWS */

    return EXIT_NORMAL;
}

void
fatal(char *s)
{
    fprintf(stderr, "fatal error: %s\n", s);
    exit(1);
}

/* There is a conflict with another cp_printword in cp/quote.c 
static void
cp_printword(s)
    char *s;
{
    printf("%s", s);
    return;
}

*/

bool
cp_getvar(char *n, int t, void *r)
{
    return (FALSE);
}

char *
cp_tildexpand(char *s)
{
	return tildexpand(s);
}
