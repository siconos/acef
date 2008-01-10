/* auxiliary functions that pass call to common CPLEX functions,
   and raise an error in case of failure.
   Some functions are activated by the OUTPUT preprocessor definition.
   -------------------------------------------------------------------
*/

#include <ilcplex/cplex.h>

/* open environment */
CPXENVptr openCPLEX (void);

/* create a problem instance */
CPXLPptr createprob (CPXENVptr env);

/* load data into a problem instance */
void copylp(CPXENVptr env, CPXLPptr lp, int ncols, int nrows, int minOrMax, double *objective,
	    double *rhs, const char *sense, const int *matbeg, const int *matcnt, const int *matind, const double *matval, double *lb, double *ub);
/* output the problem to a file in MPS format */
void writeprob(CPXENVptr env, CPXLPptr lp);

/* solve problem (optimize) */
void lpopt(CPXENVptr env, CPXLPptr lp);
/* get current solution */
void getsolution(CPXENVptr env, CPXLPptr lp, double *objval, double x[]);

/* get current basis */
void getbase(CPXENVptr env, CPXLPptr lp, int cstat[], int rstat[]);

/* reload basis */
void copybase(CPXENVptr env, CPXLPptr lp, const int cstat[], const int rstat[]);

/* change bounds */
void chgbds (CPXENVptr env, CPXLPptr lp, int nbds, int vwIndices[], char upper[], double newUppBnd[]);

/* change objective */
void chgobj (CPXENVptr env, CPXLPptr lp, int nInds, int chgObjInd[], const double chgObjVal[]);
