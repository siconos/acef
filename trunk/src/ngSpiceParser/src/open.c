
 /* Copyright 1990
   Regents of the University of California.
   All rights reserved.

   Author: 1985 Wayne A. Christopher

   The main routine for ngspice
   $Id: main.c,v 1.29 2005/05/31 16:47:48 sjborley Exp $
*/

#include <ngspice.h>

#include <stdio.h>

#ifdef HAVE_STRING_H
#include <string.h>
#endif /* HAVE_STRING_H */

#ifdef __MINGW32__
#define  srandom(a) srand(a) /* srandom */
#endif /* __MINGW32__ */

#include <setjmp.h>
#include <signal.h>
#include <sys/types.h>
#include "perform.h"

#ifdef HAVE_PWD_H
#include <pwd.h>
#endif
#ifdef HAVE_GNUREADLINE
/* Added GNU Readline Support 11/3/97 -- Andrew Veliath <veliaa@rpi.edu> */
/* from spice3f4 patch to ng-spice. jmr */
/*#include <readline/readline.h>
  #include <readline/history.h> OLIVIER */
#endif  /* HAVE_GNUREADLINE */

#ifdef HAVE_BSDEDITLINE
/* SJB added editline support 2005-05-05 */
#include <editline/readline.h>
extern VFunction *rl_event_hook;    /* missing from editline/readline.h */
extern int rl_catch_signals;	    /* missing from editline/readline.h */
#endif /* HAVE_BSDEDITLINE */

#ifndef HAVE_GETRUSAGE
#ifdef HAVE_FTIME
#include <sys/timeb.h>
#endif
#endif

#include "iferrmsg.h"/* declarations de constantes E_...*/
#include "ftedefs.h"/*declaration de la structure circ .h .h .h*/
#include "devdefs.h"/*SPICEDev declaration structure*/
#include "spicelib/devices/dev.h"
/*#include "spicelib/analysis/analysis.h" OLIVIER */
#include "misc/ivars.h"
#include "misc/getopt.h"
/*#include "frontend/resource.h" OLIVIER */
#include "frontend/variable.h"/* variables pour le parser ?*/
/* OLIVIER #include "frontend/display.h"*/  /*  added by SDB to pick up Input() fcn  */
/*#include "frontend/signal_handler.h" OLIVIER */


#if defined(HAVE_GNUREADLINE) || defined(HAVE_BSDEDITLINE)
char history_file[512];
static char *application_name;
#endif  /* HAVE_GNUREADLINE || HAVE_BSDEDITLINE */

/* Undefine this next line for dubug tracing */
/* #define TRACE */

/* Main options */
static bool ft_servermode = FALSE;
static bool ft_batchmode = FALSE;

/* Frontend options */
bool ft_intrpt = FALSE;     /* Set by the (void) signal handlers. TRUE = we've been interrupted. */
bool ft_setflag = FALSE;    /* TRUE = Don't abort simulation after an interrupt. */
char *ft_rawfile = "rawspice.raw";

#ifdef HAS_WINDOWS
bool oflag = FALSE;         /* Output over redefined I/O functions */
FILE *flogp;  /* hvogt 15.12.2001 */
#endif /* HAS_WINDOWS */

/* Frontend and circuit options */
IFsimulator *ft_sim = NULL;

/* (Virtual) Machine architecture parameters */
int ARCHme;
int ARCHsize;

char *errRtn;
char *errMsg;
char *cp_program;


struct variable *(*if_getparam)( );

static int started = FALSE;

/* static functions */
static int SIMinit(IFfrontEnd *frontEnd, IFsimulator **simulator);
static int shutdown(int exitval);
static void app_rl_readlines(void);

#if defined(HAVE_GNUREADLINE) || defined(HAVE_BSDEDITLINE)
static char * prompt(void);
#endif /* HAVE_GNUREADLINE || HAVE_BSDEDITLINE */

#ifdef HAVE_GNUREADLINE
static int rl_event_func(void) ;
#endif /* HAVE_GNUREADLINE || HAVE_BSDEDITLINE */
#ifdef HAVE_BSDEDITLINE
static void rl_event_func(void) ;
#endif /* HAVE_BSDEDITLINE */

static void show_help(void);
static void show_version(void);
static bool read_initialisation_file(char * dir, char * name);
#ifdef SIMULATOR
static void append_to_stream(FILE *dest, FILE *source);
#endif /* SIMULATOR */


#ifndef HAVE_GETRUSAGE
#ifdef HAVE_FTIME
extern struct timeb timebegin;		/* for use w/ ftime */
#endif
#endif

extern IFsimulator SIMinfo;

#ifdef SIMULATOR

bool ft_nutmeg = FALSE;
extern struct comm spcp_coms[ ];
struct comm *cp_coms = spcp_coms;

#else /* SIMULATOR */

bool ft_nutmeg = TRUE;
extern struct comm nutcp_coms[ ];
struct comm *cp_coms = nutcp_coms;
static IFfrontEnd nutmeginfo;

/* -------------------------------------------------------------------------- */
int
if_run(char *t, char *w, wordlist *s, char *b)
{
    return (0);
}

/* -------------------------------------------------------------------------- */
int
if_sens_run(char *t, char *w, wordlist *s, char *b)
{
    return (0);
}

/* -------------------------------------------------------------------------- */
void
if_dump(char *ckt, FILE *fp)
{}

/* -------------------------------------------------------------------------- */
char *
if_inpdeck(struct line *deck, char **tab)
{
    return ((char *) 0);
}

/* -------------------------------------------------------------------------- */
int
if_option(char *ckt, char *name, int type, char *value)
{
    return 0;
}

/* -------------------------------------------------------------------------- */
void if_cktfree(char *ckt, char *tab)
{}

/* -------------------------------------------------------------------------- */
void if_setndnames(char *line)
{}

/* -------------------------------------------------------------------------- */
char *
if_errstring(int code)
{
    return ("spice error");
}

/* -------------------------------------------------------------------------- */
void
if_setparam(char *ckt, char *name, char *param, struct variable *val)
{}

/* -------------------------------------------------------------------------- */
bool
if_tranparams(struct circ *ckt, double *start, double *stop, double *step)
{
    return (FALSE); 
}

/* -------------------------------------------------------------------------- */
struct variable *
if_getstat(char *n, char *c)
{
    return (NULL);
}

#ifdef EXPERIMENTAL_CODE
void com_loadsnap(wordlist *wl) { return; }
void com_savesnap(wordlist *wl) { return; }
#endif

#endif /* SIMULATOR */

#ifndef SIMULATOR


#endif /* SIMULATOR */

char *hlp_filelist[] = { "ngspice", 0 };


/* allocate space for global constants in 'CONST.h' */

double CONSTroot2;
double CONSTvt0;
double CONSTKoverQ;
double CONSTe;
IFfrontEnd *SPfrontEnd = NULL;
int DEVmaxnum = 0;

/* -------------------------------------------------------------------------- */
static int
SIMinit(IFfrontEnd *frontEnd, IFsimulator **simulator)
{
#ifdef SIMULATOR
    spice_init_devices();
    SIMinfo.numDevices = DEVmaxnum = num_devices();
    SIMinfo.devices = devices_ptr();
    SIMinfo.numAnalyses = spice_num_analysis();
     SIMinfo.analyses = (IFanalysis **)spice_analysis_ptr(); /* va: we recast, because we use 
                                                             * only the public part 
							     */
							
#ifdef CIDEROLIVIER
/* Evaluates limits of machine accuracy for CIDER */
     /* OLIVIER evalAccLimits();*/
#endif /* CIDER */  
   						     
#endif /* SIMULATOR */
     /*UTILISE DANS LES Devices*/
    SPfrontEnd = frontEnd;
    *simulator = &SIMinfo;
    CONSTroot2 = sqrt(2.);
    CONSTvt0 = CONSTboltz * (27 /* deg c */ + CONSTCtoK ) / CHARGE;
    CONSTKoverQ = CONSTboltz / CHARGE;
    CONSTe = exp((double)1.0);
    return(OK);
}


/* -------------------------------------------------------------------------- */
/* Shutdown gracefully. */
static int
shutdown(int exitval)
{
    cleanvars();
#ifdef PARALLEL_ARCH
    if (exitval == EXIT_BAD) {
	Error("Fatal error in SPICE", -1);
    } else {
	PEND_();
    }
#endif /* PARALLEL_ARCH */
    exit (exitval);
}

/* -------------------------------------------------------------------------- */

#if defined(HAVE_GNUREADLINE) || defined(HAVE_BSDEDITLINE)
/* Adapted ../lib/cp/lexical.c:prompt() for GNU Readline -- Andrew Veliath <veliaa@rpi.edu> */
static char *
prompt(void)
{
    static char pbuf[128];
    char *p = pbuf, *s;

    if (cp_interactive == FALSE)
        return NULL;	/* NULL means no prompt */
    
    s = get_alt_prompt();
    if(s==NULL)
	s = cp_promptstring;
    if(s==NULL)
	s = "->";
    
    while (*s) {
	switch (strip(*s)) {
	    case '!':
#ifdef HAVE_BSDEDITLINE
		{
		    /* SJB In the present version of editline (v2.9)
		      it seems that where_history() is broken.
		      This is a hack that works round this problem.
		      WARNING: It may fail to work in the future
		      as it relies on undocumented structure */
		    int where = 0;
		    HIST_ENTRY * he = current_history();
		    if(he!=NULL) where = *(int*)(he->data);
		    p += sprintf(p, "%d", where + 1);
		}
#else
		p += sprintf(p, "%d", where_history() + 1);
#endif	/* HAVE_BSDEDITLINE*/
		break;
	    case '\\':
		if (*(s + 1))
		    p += sprintf(p, "%c", strip(*++s));
	    default:
		*p = strip(*s); ++p;
		break;
	}
        s++;
    }
    *p = 0;
    return pbuf;
}
#endif /* HAVE_GNUREADLINE || HAVE_BSDEDITLINE */

#ifdef HAVE_GNUREADLINE
/* -------------------------------------------------------------------------- */
/* Process device events in Readline's hook since there is no where
   else to do it now - AV */
static int
rl_event_func()  
/* called by GNU readline periodically to know what to do about keypresses */
{
    static REQUEST reqst = { checkup_option, 0 };
    Input(&reqst, NULL);
    return 0;
}
#endif /* HAVE_GNUREADLINE */

#ifdef HAVE_BSDEDITLINE
/* -------------------------------------------------------------------------- */
/* Process device events in Editline's hook.
   similar to the readline function above but returns void */
static void
rl_event_func()  
/* called by GNU readline periodically to know what to do about keypresses */
{
    static REQUEST reqst = { checkup_option, 0 };
    Input(&reqst, NULL);
}
#endif /* HAVE_BSDEDITLINE */

/* -------------------------------------------------------------------------- */
/* This is the command processing loop for spice and nutmeg.
   The function is called even when GNU readline is unavailable, in which
   case it falls back to repeatable calling cp_evloop()
   SJB 26th April 2005 */
static void
app_rl_readlines()
{
#if defined(HAVE_GNUREADLINE) || defined(HAVE_BSDEDITLINE)
    /* GNU Readline Support -- Andrew Veliath <veliaa@rpi.edu> */
    char *line, *expanded_line;
    
    /* ---  set up readline params --- */
    strcpy(history_file, getenv("HOME"));
    strcat(history_file, "/.");
    strcat(history_file, application_name);
    strcat(history_file, "_history");
    
    using_history();
    read_history(history_file);
    
    rl_readline_name = application_name;
    rl_instream = cp_in;
    rl_outstream = cp_out;
    rl_event_hook = rl_event_func;
    rl_catch_signals = 0;   /* disable signal handling  */
    
    /* sjb - what to do for editline?
       This variable is not supported by editline. */   
#if defined(HAVE_GNUREADLINE) 
    rl_catch_sigwinch = 1;  /* allow readline to respond to resized windows  */
#endif  
    
    /* note that we want some mechanism to detect ctrl-D and expand it to exit */
    while (1) {
       history_set_pos(history_length);

       SETJMP(jbuf, 1);    /* Set location to jump to after handling SIGINT (ctrl-C)  */

       line = readline(prompt());
       if (line && *line) {
           int s = history_expand(line, &expanded_line);

           if (s == 2) {
               fprintf(stderr, "-> %s\n", expanded_line);
           } else if (s == -1) {
               fprintf(stderr, "readline: %s\n", expanded_line);
           } else {
               cp_evloop(expanded_line);
               add_history(expanded_line);
           }
           free(expanded_line);
       }
       if (line) free(line);
    }
    /* History gets written in ../fte/misccoms.c com_quit */
    
#else
    while (cp_evloop((char *) NULL) == 1) ;
#endif /* defined(HAVE_GNUREADLINE) || defined(HAVE_BSDEDITLINE) */
}


/* -------------------------------------------------------------------------- */
static void
show_help(void)
{
    printf("Usage\n");
}

/* -------------------------------------------------------------------------- */
static void
show_version(void)
{
    printf("");
}

#ifdef SIMULATOR
/* -------------------------------------------------------------------------- */
static void
append_to_stream(FILE *dest, FILE *source)
{
    char *buf[BSIZE_SP];
    int i;

    while ((i = fread(buf, 1, BSIZE_SP, source)) > 0)
	fwrite(buf, i, 1, dest);
}
#endif /* SIMULATOR */

/* -------------------------------------------------------------------------- */
/* Read an initialisation file.
   dir    is the directory (use NULL or "" for current directory)
   name   is the initialisation file's name
   Return true on success
   SJB 25th April 2005 */
static bool
read_initialisation_file(char * dir, char * name)
{
#ifndef HAVE_ASPRINTF
    FILE * fp = NULL;
#endif /* not HAVE_ASPRINTF */
    char * path;
    bool result = FALSE;
    
    /* check name */
    if(name==NULL || name[0]=='\0')
    	return FALSE;	/* Fail; name needed */
    
    /* contruct the full path */
    if(dir == NULL || dir[0]=='\0') {
	path = name;
    } else {
#ifdef HAVE_ASPRINTF
	asprintf(&path, "%s" DIR_PATHSEP "%s", dir,name);
	if(path==NULL) return FALSE;	/* memory allocation error */
#else /* ~ HAVE_ASPRINTF */
	path=(char*)tmalloc(1 + strlen(dir)+strlen(name));
	if(path==NULL) return FALSE;	/* memory allocation error */
	sprintf(path,"%s" DIR_PATHSEP "%s",dir,name);
#endif /* HAVE_ASPRINTF */
    }

    /* now access the file */
#ifdef HAVE_UNISTD_H
    if (access(path, R_OK) == 0) {		
#else
    if ((fp = fopen(path, "r")) != NULL) {
	(void) fclose(fp);
#endif /* HAVE_UNISTD_H */
	inp_source(path);
#ifdef TRACE
	printf("Init file: '%s'\n",path);
#endif /* TRACE */	
	result = TRUE;	/* loaded okay */
    }
    
    /* if dir was not NULL and not empty then we allocated memory above */
    if(dir!=NULL && dir[0] !='\0')
#ifdef HAVE_ASPRINTF
    	free(path);
#else
	tfree(path);
#endif /* HAVE_ASPRINTF */
    
    return result;
}

/* -------------------------------------------------------------------------- */

#ifdef SIMULATOR
extern int OUTpBeginPlot(), OUTpData(), OUTwBeginPlot(), OUTwReference();
extern int OUTwData(), OUTwEnd(), OUTendPlot(), OUTbeginDomain();
extern int OUTendDomain(), OUTstopnow(), OUTerror(), OUTattributes();
#endif /* SIMULATOR */
 int init(){
          int c;
    int		err;
    bool	gotone = FALSE;
    char*       copystring;/*DG*/
#ifdef SIMULATOR
    int error2;
    
    static IFfrontEnd nutmeginfo = {
	IFnewUid,
	IFdelUid,
	OUTstopnow,
	seconds,
	OUTerror,
	OUTpBeginPlot,
	OUTpData,
	OUTwBeginPlot,
	OUTwReference,
	OUTwData,
	OUTwEnd,
	OUTendPlot,
	OUTbeginDomain,
	OUTendDomain,
	OUTattributes
    };
#else  /* ~ SIMULATOR */
    bool gdata = TRUE;
#endif /* ~ SIMULATOR */

    char buf[BSIZE_SP];
    bool readinit = TRUE;
    bool rflag = FALSE;
    bool istty = TRUE;
    bool iflag = FALSE;
    bool qflag = FALSE;
    FILE *fp;
    FILE *circuit_file;
    

    /* MFB tends to jump to 0 on errors.  This tends to catch it. */
    if (started) {
        fprintf(cp_err, "main: Internal Error: jump to zero\n");
        shutdown(EXIT_BAD);
    }
    started = TRUE;

    ARCHme = 0;
    ARCHsize = 1;

    ivars( );

    cp_in = stdin;
    cp_out = stdout;
    cp_err = stderr;

    circuit_file = stdin;
    printf("rep1\n");
#ifdef MALLOCTRACE
    mallocTraceInit("malloc.out");
#endif
#if defined(HAVE_ISATTY) && !defined(HAS_WINDOWS)
    istty = (bool) isatty(fileno(stdin));
#endif

    /*init_time( ); OLIVIER */
     /*initialise les pointer de fcts, et qq constantes utiles lors la conversion*/
    err = SIMinit(&nutmeginfo,&ft_sim);
    if(err != OK) {
        ft_sperror(err,"SIMinit");
        shutdown(EXIT_BAD);
    }
    cp_program = ft_sim->simulator;

    srandom(getpid());
    printf("rep2\n");
    /* Have to initialize cp now. */
    ft_cpinit();
 }
 
int readFile(char *file){

    int		err;
    init();
    cp_interactive = FALSE;
    err = 0;
    printf("readfile\n");
#ifdef SIMULATOR
   
	FILE *tempfile;

	tempfile = tmpfile();
	    char *arg;
	    FILE *tp;

	    /* Copy all the arguments into the temporary file */
	    
	    tp = fopen(file, "r");
	    printf("main rep6 fopen : ");
	    printf(file);
	    printf("\n");
	    if (!tp) {
		perror(file);
		err = 1;
		return shutdown(EXIT_NORMAL);
	    }
	    append_to_stream(tempfile, tp);
	    fclose(tp);
	fseek(tempfile, (long) 0, 0);

        if (tempfile ) {
	  inp_spsource(tempfile, FALSE, NULL);
	
        }
#endif
	printf("readfile end\n");
    
 }
