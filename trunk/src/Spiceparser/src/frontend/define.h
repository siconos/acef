/*************
 * Header file for define.c
 * 1999 E. Rouat
 ************/

#ifndef DEFINE_H_INCLUDED
#define DEFINE_H_INCLUDED

void com_define(wordlist *wlist);
struct pnode * ft_substdef(char *name, struct pnode *args);
void com_undefine(wordlist *wlist);



#endif
