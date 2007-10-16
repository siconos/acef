/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 1985 Wayne A. Christopher, U. C. Berkeley CAD Group
**********/

/*
 *
 * Various things for quoting words. If this is not ascii, quote and
 * strip are no-ops, so '' and \ quoting won't work. To fix this, sell
 * your IBM machine and buy a vax.
 */

#include "ngspice.h"
#include "cpdefs.h"
#include "quote.h"


/* Strip all the 8th bits from a string (destructively). */

void
cp_wstrip(char *str)
{
    char c, d;

    if (str)
      while ((c = *str)) {   /* assign and test */
	    d = strip(c);
	    if (c != d)
		    *str = d;
	    str++;
	}
    return;
}

/* Quote all characters in a word. */

void
cp_quoteword(char *str)
{
    if (str)
	while (*str) {
	    *str = quote(*str);
	    str++;
	}
    return;
}

/* Print a word (strip the word first). */

void
cp_printword(char *string, FILE *fp)
{
    char *s;

    if (string)
        for (s = string; *s; s++)
            (void) putc((strip(*s)), fp);
    return;
}

/* (Destructively) strip all the words in a wlist. */

void
cp_striplist(wordlist *wlist)
{
    wordlist *wl;

    for (wl = wlist; wl; wl = wl->wl_next)
        cp_wstrip(wl->wl_word);
    return;
}

/* Remove the "" from a string. */

char *
cp_unquote(char *string)
{
    char *s;
    int l;
    if (string) {
	l = strlen(string);
	s = MALLOC(l+1);
	
	if (*string == '"' && string[l-1] == '"') {
	    strncpy(s,string+1,l-2);
	    s[l-2] = '\0';
	} else
	    strcpy(s,string);

	return (s);
    } else
	return 0;
}
