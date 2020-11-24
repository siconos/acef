/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 1986 Wayne A. Christopher, U. C. Berkeley CAD Group
**********/

/*
 * Routines to handle "more"d output.  There are some serious system
 * dependencies in here, and it isn't clear that versions of this stuff
 * can be written for every possible machine...
 */
#include <config.h>

#ifdef HAVE_SGTTY_H
#include <sgtty.h>
#endif

#if 0
/* Bad interaction with bool type in bool.h because curses also
   defines this symbol. */
#ifdef HAVE_TERMCAP
#include <curses.h>
#include <term.h>
#endif
#endif

#ifdef HAVE_TERMCAP
#include <termcap.h>
#endif

#include <ngspice.h>
#include <cpdefs.h>

#include "variable.h"
#include "terminal.h"

static char *motion_chars;
static char *clear_chars;
static char *home_chars;
static char *cleol_chars;


#define DEF_SCRHEIGHT   24
#define DEF_SCRWIDTH    80

bool out_moremode = TRUE;
bool out_isatty = TRUE;

static int xsize, ysize;
static int xpos, ypos;
static bool noprint, nopause;


/* out_printf doesn't handle double arguments correctly, so we
    sprintf into this buf and call out_send w/ it */
char out_pbuf[BSIZE_SP];

/* Start output... */

void
out_init(void)
{
#ifdef TIOCGWINSZ
    struct winsize ws;
#endif
    bool moremode;

    noprint = nopause = FALSE;

    if (cp_getvar("nomoremode", VT_BOOL, (char *) &moremode))
        out_moremode = FALSE;
    else
    out_moremode = TRUE;
    if (!out_moremode || !cp_interactive)
        out_isatty = FALSE;

    if (!out_isatty)
        return;

    xsize = ysize = 0;

    /* Figure out the screen size.  We try, in order, TIOCGSIZE,
     * tgetent(), and cp_getvar(height).  Default is 24 x 80.
     */

#ifdef TIOCGWINSZ
    if (!xsize || !ysize) {
        (void) ioctl(fileno(stdout), TIOCGWINSZ, (char *) &ws);
        xsize = ws.ws_col;
        ysize = ws.ws_row;
    }
#endif

    if (!xsize)
        (void) cp_getvar("width", VT_NUM, (char *) &xsize);
    if (!ysize)
        (void) cp_getvar("height", VT_NUM, (char *) &ysize);

    if (!xsize)
        xsize = DEF_SCRWIDTH;
    if (!ysize)
        ysize = DEF_SCRHEIGHT;
    ysize -= 2; /* Fudge room... */
    xpos = ypos = 0;

    return;
}

/* Putc may not be buffered (sp?), so we do it ourselves. */

static char staticbuf[BUFSIZ];
struct {
    int count;
    char *ptr;
} ourbuf = { BUFSIZ, staticbuf };

/* send buffer out */
void
outbufputc(void)
{

    if (ourbuf.count != BUFSIZ) {
      fputs(staticbuf, cp_out);
      memset(staticbuf, 0, BUFSIZ-ourbuf.count);
      ourbuf.count = BUFSIZ;
      ourbuf.ptr = staticbuf;
    }

}

static void
bufputc(char c)
{
    if (--ourbuf.count >= 0) {
	*ourbuf.ptr++ = c;
    } else {
	/* Flush and reset the buffer */
	outbufputc();
	/* and store the character. */
	ourbuf.count--;
	*ourbuf.ptr++ = c;
    }
}


/* prompt for a return */
void
promptreturn(void)
{
    char buf[16];
moe:
    fprintf(cp_out,
        "\n\t-- hit return for more, ? for help -- ");
    if (!fgets(buf, 16, cp_in)) {
        clearerr(cp_in);
        *buf = 'q';
    }
    switch (*buf) {
        case '\n':
            break;
        case 'q':
            noprint = TRUE;
            break;
        case 'c':
            nopause = TRUE;
            break;
        case ' ':
            break;
        case '?':
            fprintf(cp_out,
"\nPossible responses:\n\
\t<cr>   : Print another screenful\n\
\tq <cr> : Discard the rest of the output\n\
\tc <cr> : Continuously print the rest of the output\n\
\t? <cr> : Print this help message\n");
            goto moe;
        default:
            fprintf(cp_out, "Character %d is no good\n", *buf);
            goto moe;
    }

}

/* Print a string to the output.  If this would cause the screen to scroll,
 * print "more".
 */

void
out_send(char *string)
{

    if (noprint)
        return;
    if (!out_isatty || nopause) {
        fputs(string, cp_out);
        return;
    }
    while (*string) {
        switch (*string) {
            case '\n':
                xpos = 0;
                ypos++;
                break;
            case '\f':
                ypos = ysize;
                xpos = 0;
                break;
            case '\t':
                xpos = xpos / 8 + 1;
                xpos *= 8;
                break;
            default:
                xpos++;
                break;
        }
        while (xpos >= xsize) {
            xpos -= xsize;
            ypos++;
        }
        if (ypos >= ysize) {
            outbufputc();       /* out goes buffer */
            promptreturn();
            (void) fflush(cp_out);
            ypos = xpos = 0;
        }
        bufputc(*string);   /* we need to buffer these */
        string++;
    }
    (void) outbufputc();
    return;
}

/* Printf some stuff using more mode. */

#define MAXLEN 4096


void
out_printf(char *fmt, char *s1, char *s2, char *s3, char *s4, char *s5, char *s6, char *s7, char *s8, char *s9, char *s10)
{
    char buf[MAXLEN];

    sprintf(buf, fmt, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10);

    out_send(buf);
    return;
}

static int
outfn(int c)
{
	putc(c, stdout);
	return c;
}


void 
tcap_init(void)
{
    char *s;
#ifdef HAVE_TERMCAP
    char tbuf[1025];
    static char buf2[100];
    char *charbuf;

    charbuf = buf2;

    if ((s = getenv("TERM"))) {
	if (tgetent(tbuf, s) != -1) {
	    xsize = tgetnum("co");
	    ysize = tgetnum("li");
	    if ((xsize <= 0) || (ysize <= 0))
		xsize = ysize = 0;
	    clear_chars = (char *) tgetstr("cl", &charbuf);
	    motion_chars = (char *) tgetstr("cm", &charbuf);
	    home_chars = (char *) tgetstr("ho", &charbuf);
	    cleol_chars = (char *) tgetstr("ce", &charbuf);
	}
    }
#endif

    if (!xsize) {
        if ((s = getenv("COLS")))
            xsize = atoi(s);
        if (xsize <= 0)
            xsize = 0;
    }

    if (!ysize) {
        if ((s = getenv("LINES")))
            ysize = atoi(s);
        if (ysize <= 0)
            ysize = 0;
    }
}


void 
term_clear(void)
{
#ifdef HAVE_TERMCAP
    if (*clear_chars)
	tputs(clear_chars, 1, outfn);
    else
	fputs("\n", stdout);
#endif
}


void 
term_home(void)
{
#ifdef HAVE_TERMCAP
    if (*home_chars)
	tputs(home_chars, 1, outfn);
    else if (*motion_chars)
	tputs(tgoto(motion_chars, 1, 1), 1, outfn);
    else
	fputs("\n", stdout);
#endif
}


void 
term_cleol(void)
{
#ifdef HAVE_TERMCAP
    if (*cleol_chars)
	tputs(cleol_chars, 1, outfn);
#endif
}
