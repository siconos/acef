/* $Id: pnode.h,v 1.2 2005/05/26 19:30:53 sjborley Exp $ */
#ifndef _PNODE_H
#define _PNODE_H

struct pnode {
    char *pn_name;		/* If non-NULL, the name. */
    struct dvec *pn_value;	/* Non-NULL in a terminal node. */
    struct func *pn_func;	/* Non-NULL is a function. */
    struct op *pn_op;		/* Operation if the above two NULL. */
    struct pnode *pn_left;	/* Left branch or function argument. */
    struct pnode *pn_right;	/* Right branch. */
    struct pnode *pn_next;	/* For expression lists. */
    int pn_use;			/* usage counter */
} ;

#endif
