#include <stddef.h>

#include <dvec.h>
#include <ngspice.h>
#include <fteext.h>
#include <parse.h>
#include <cpextern.h>

#include "com_let.h"
#include "com_display.h"
#include "completion.h"


void
com_let(wordlist *wl)
{
    char *p, *q, *s;
    int	indices[MAXDIMS];
    int	numdims;
    wordlist fake_wl;
    int need_open;
    int offset, length;
    struct pnode *nn;
    struct dvec *n, *t;
    int i, cube;
    int j, depth;
    int newvec;
    char *rhs;

    fake_wl.wl_next = NULL;

    if (!wl) {
        com_display((wordlist *) NULL);
        return;
    }

    p = wl_flatten(wl);

    /* extract indices */
    numdims = 0;
    if ((rhs =strchr(p, '='))) {
	*rhs++ = 0;
    } else {
	fprintf(cp_err, "Error: bad let syntax\n");
	tfree(p);
	return;
    }
    if ((s =strchr(p, '['))) {
	need_open = 0;
	*s++ = 0;
	while (!need_open || *s == '[') {
	    depth = 0;
	    if (need_open)
		s++;
	    for (q = s; *q && (*q != ']' && (*q != ',' || depth > 0)); q++) {
		switch (*q) {
		case '[':
		    depth += 1;
		    break;
		case ']':
		    depth -= 1;
		    break;
		}
	    }

	    if (depth != 0 || !*q) {
		printf("syntax error specifying index\n");
		tfree(p);
		return;
	    }

	    if (*q == ']')
		need_open = 1;
	    else
		need_open = 0;
	    if (*q)
	        *q++ = 0;

	    /* evaluate expression between s and q */
	    /* va, indexing */
	    fake_wl.wl_word = s;
	    nn = ft_getpnames(&fake_wl, TRUE);
	    t = ft_evaluate(nn);

	    if (!isreal(t) || t->v_link2 || t->v_length != 1 || !t->v_realdata)
	    {
		fprintf(cp_err, "Error: index is not a scalar.\n");
		goto quit;
	    }
	    j = (int)floor(t->v_realdata[0]+0.5); /* ignore sanity checks for now, va, which checks? */

	    if (j < 0) {
		printf("negative index (%d) is not allowed\n", j);
		goto quit;
	    }

	    indices[numdims++] = j;

	    /* va: garbage collection for t, if pnode nn is no simple value */
	    if (nn!=NULL && nn->pn_value==NULL && t!=NULL) vec_free(t);
	    free_pnode(nn); /* frees also t, if pnode nn is simple value */

	    for (s = q; *s && isspace(*s); s++)
		;
	}
    }
    /* vector name at p */

    for (q = p + strlen(p) - 1; *q <= ' ' && p <= q; q--)
	;

    *++q = 0;

    /* sanity check */
    if (eq(p, "all") ||strchr(p, '@')) {
        fprintf(cp_err, "Error: bad variable name %s\n", p);
        tfree(p);
        return;
    }

    /* evaluate rhs */
    fake_wl.wl_word = rhs;
    nn = ft_getpnames(&fake_wl, TRUE);
    if (nn == NULL) {
	/* XXX error message */
        tfree(p);
        return;
    }
    t = ft_evaluate(nn);
    if (!t) {
	fprintf(cp_err, "Error: Can't evaluate %s\n", rhs);
        free_pnode(nn);
        tfree(p);
        return;
    }

    if (t->v_link2)
        fprintf(cp_err, "Warning: extra wildcard values ignored\n");

    n = vec_get(p);

    if (n) {
	/* re-allocate? */
	/* vec_free(n); */
	newvec = 0;
    } else {
	if (numdims) {
	    fprintf(cp_err, "Can't assign into a subindex of a new vector\n");
	    goto quit;
	}

	/* create and assign a new vector */
	n = alloc(struct dvec);
	ZERO(n, struct dvec);
	n->v_name = copy(p);
	n->v_type = t->v_type;
	n->v_flags = (t->v_flags | VF_PERMANENT);
	n->v_length = t->v_length;

	if (!t->v_numdims) {
	    n->v_numdims = 1;
	    n->v_dims[0] = n->v_length;
	} else {
	    n->v_numdims = t->v_numdims;
	    for (i = 0; i < t->v_numdims; i++)
		n->v_dims[i] = t->v_dims[i];
	}

	if (isreal(t))
	    n->v_realdata = (double *) tmalloc(n->v_length * sizeof(double));
	else
	    n->v_compdata = (complex *) tmalloc(n->v_length * sizeof(complex));
	newvec = 1;
	vec_new(n);
    }

    /* fix-up dimensions; va, also for v_dims */
    if (n->v_numdims < 1 || n->v_dims[0]==0 ) {
	n->v_numdims = 1;
	n->v_dims[0] = n->v_length;
    }

    /* Compare dimensions */
    offset = 0;
    length = n->v_length;

    cube = 1;
    for (i = n->v_numdims - 1; i >= numdims; i--)
	cube *= n->v_dims[i];

    for (i = numdims - 1; i >= 0; i--) {
	offset += cube * indices[i];
	if (i < n->v_numdims) {
	    cube *= n->v_dims[i];
	    length /= n->v_dims[i];
	}
    }

    /* length is the size of the unit refered to */
    /* cube ends up being the length */

    if (length > t->v_length) {
	fprintf(cp_err, "left-hand expression is too small (need %d)\n",
		length * cube);
	if (newvec)
	    n->v_flags &= ~VF_PERMANENT;
	goto quit;
    }
    if (isreal(t) != isreal(n)) {
	fprintf(cp_err,
		"Types of vectors are not the same (real vs. complex)\n");
	if (newvec)
	    n->v_flags &= ~VF_PERMANENT;
	goto quit;
    } else if (isreal(t)) {
	bcopy((char *) t->v_realdata, (char *) (n->v_realdata + offset),
		length * sizeof (double));
    } else {
	bcopy((char *) t->v_compdata, (char *) (n->v_compdata + offset),
		length * sizeof (complex));
    }

    n->v_minsignal = 0.0; /* How do these get reset ??? */
    n->v_maxsignal = 0.0;

    n->v_scale = t->v_scale;

    if (newvec)
	cp_addkword(CT_VECTOR, n->v_name);

quit:
    /* va: garbage collection for t, if pnode nn is no simple value */
    if (nn!=NULL && nn->pn_value==NULL && t!=NULL) vec_free(t);
    free_pnode(nn); /* frees also t, if pnode nn is simple value */
    tfree(p);
    return;
}
