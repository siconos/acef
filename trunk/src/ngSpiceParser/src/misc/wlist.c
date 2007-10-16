/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 1985 Wayne A. Christopher, U. C. Berkeley CAD Group
**********/

/* Wordlist manipulation stuff.  */

#include <config.h>
#include <ngspice.h>
#include <bool.h>
#include <wordlist.h>


/* Determine the length of a word list. */
int
wl_length(wordlist *wlist)
{
    int i = 0;
    wordlist *wl;

    for (wl = wlist; wl; wl = wl->wl_next)
        i++;
    return (i);
}


/* Free the storage used by a word list. */
void
wl_free(wordlist *wlist)
{
    wordlist *wl, *nw;

    for (wl = wlist; wl; wl = nw) {
        nw = wl->wl_next;
        tfree(wl->wl_word);
        tfree(wl);
    }
    return;
}


/* Copy a wordlist and the words. */
wordlist *
wl_copy(wordlist *wlist)
{
    wordlist *wl, *nwl = NULL, *w = NULL;

    for (wl = wlist; wl; wl = wl->wl_next) {
        if (nwl == NULL) {
            nwl = w = alloc(struct wordlist);
	    w->wl_prev = NULL;
	    w->wl_next = NULL;
        } else {
            w->wl_next = alloc(struct wordlist);
            w->wl_next->wl_prev = w;
            w = w->wl_next;
	    w->wl_next = NULL;
        }
        w->wl_word = copy(wl->wl_word);
    }
    return (nwl);
}


/* Substitute a wordlist for one element of a wordlist, and return a
 * pointer to the last element of the inserted list.  */
wordlist *
wl_splice(wordlist *elt, wordlist *list)
{

    if (list)
        list->wl_prev = elt->wl_prev;
    if (elt->wl_prev)
        elt->wl_prev->wl_next = list;
    if (list) {
        while (list->wl_next)
            list = list->wl_next;
        list->wl_next = elt->wl_next;
    }
    if (elt->wl_next)
        elt->wl_next->wl_prev = list;
    tfree(elt->wl_word);
    tfree(elt);
    return (list);
}



static void
printword(char *string, FILE *fp)
{
    char *s;

    if (string)
        for (s = string; *s; s++)
            putc((strip(*s)), fp);
    return;
}


/* Print a word list. (No \n at the end...) */
void
wl_print(wordlist *wlist, FILE *fp)
{
    wordlist *wl;

    for (wl = wlist; wl; wl = wl->wl_next) {
        printword(wl->wl_word, fp);
        if (wl->wl_next)
            putc(' ', fp);
    }
    return;
}


/* Turn an array of char *'s into a wordlist. */
wordlist *
wl_build(char **v)
{
    wordlist *wlist = NULL;
    wordlist *wl = NULL;
    wordlist *cwl;

    while (*v) {
        cwl = alloc(struct wordlist);
        cwl->wl_prev = wl;
        if (wl)
            wl->wl_next = cwl;
        else {
            wlist = cwl;
            cwl->wl_next = NULL;
	}
        cwl->wl_word = copy(*v);
        wl = cwl;
        v++;
    }
    return (wlist);
}

char **
wl_mkvec(wordlist *wl)
{
    int len, i;
    char **v;

    len = wl_length(wl);
    v = (char **) tmalloc((len + 1) * sizeof (char **));
    for (i = 0; i < len; i++) {
        v[i] = copy(wl->wl_word);
        wl = wl->wl_next;
    }
    v[i] = NULL;
    return (v);
}


/* Nconc two wordlists together. */
wordlist *
wl_append(wordlist *wlist, wordlist *nwl)
{
    wordlist *wl;
    if (wlist == NULL)
        return (nwl);
    if (nwl == NULL)
        return (wlist);
    for (wl = wlist; wl->wl_next; wl = wl->wl_next);
    wl->wl_next = nwl;
    nwl->wl_prev = wl;
    return (wlist);
}


/* Reverse a word list. */
wordlist *
wl_reverse(wordlist *wl)
{
    wordlist *w, *t;

    for (w = wl; ; w = t) {
         t = w->wl_next;
         w->wl_next = w->wl_prev;
         w->wl_prev = t;
         if (t == NULL)
            break;
    }
    return (w);
}


/* Convert a wordlist into a string. */
char *
wl_flatten(wordlist *wl)
{
    char *buf;
    wordlist *tw;
    int i = 0;

    for (tw = wl; tw; tw = tw->wl_next)
        i += strlen(tw->wl_word) + 1;
    buf = tmalloc(i + 1);
    *buf = 0;

    while (wl != NULL) {
        (void) strcat(buf, wl->wl_word);
        if (wl->wl_next)
            (void) strcat(buf, " ");
        wl = wl->wl_next;
    }
    return (buf);
}


/* Return the nth element of a wordlist, or the last one if n is too
 * big.  Numbering starts at 0...  */
wordlist *
wl_nthelem(int i, wordlist *wl)
{
    wordlist *ww = wl;

    while ((i-- > 0) && ww->wl_next)
        ww = ww->wl_next;
    return (ww);
}



static int
wlcomp(const void *a, const void *b)
{
    char **s = (char **) a;
    char **t = (char **) b;
    return (strcmp(*s, *t));
}


void
wl_sort(wordlist *wl)
{
    int i = 0;
    wordlist *ww = wl;
    char **stuff;

    for (i = 0; ww; i++)
        ww = ww->wl_next;
    if (i < 2)
        return;
    stuff = (char **) tmalloc(i * sizeof (char *));
    for (i = 0, ww = wl; ww; i++, ww = ww->wl_next)
        stuff[i] = ww->wl_word;
    qsort((char *) stuff, i, sizeof (char *), wlcomp);
    for (i = 0, ww = wl; ww; i++, ww = ww->wl_next)
        ww->wl_word = stuff[i];
    tfree(stuff);
    return;
}


/* Return a range of wordlist elements... */
wordlist *
wl_range(wordlist *wl, int low, int up)
{
    int i;
    wordlist *tt;
    bool rev = FALSE;

    if (low > up) {
        i = up;
        up = low;
        low = i;
        rev = TRUE;
    }
    up -= low;
    while (wl && (low > 0)) {
        tt = wl->wl_next;
        tfree(wl->wl_word);
        tfree(wl);
        wl = tt;
        if (wl)
            wl->wl_prev = NULL;
        low--;
    }
    tt = wl; 
    while (tt && (up > 0)) {
        tt = tt->wl_next;
        up--; 
    } 
    if (tt && tt->wl_next) {
        wl_free(tt->wl_next);
        tt->wl_next = NULL;
    }
    if (rev)
        wl = wl_reverse(wl);
    return (wl);
}

