/*************
 * Header file for postcoms.c
 * 1999 E. Rouat
 ************/

#ifndef POSTCOMS_H_INCLUDED
#define POSTCOMS_H_INCLUDED

void com_let(wordlist *wl);
void com_unlet(wordlist *wl);
void com_load(wordlist *wl);
void com_print(wordlist *wl);
void com_write(wordlist *wl);
void com_transpose(wordlist *wl);
void com_setscale(wordlist *wl);
void com_display(wordlist *wl);
void com_cross(wordlist *wl);
void com_destroy(wordlist *wl);
void com_splot(wordlist *wl);


#endif
