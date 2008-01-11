/* auxiliary functions that pass call to common CPLEX functions,
   and raise an error in case of failure.
   Some functions are activated by the OUTPUT preprocessor definition.
   -------------------------------------------------------------------
*/

#include <ilcplex/cplex.h>

/* open environment */
CPXENVptr openCPLEX (void) {
	int status;
	CPXENVptr env = CPXopenCPLEX (&status);
	if (env == NULL) {
		char  errmsg[1024];
		printf ("Could not open CPLEX environment.\n");
		CPXgeterrorstring (env, status, errmsg);
		printf("%s", errmsg);
		throw(1);
	}
	return env;
}

/* create a problem instance */
CPXLPptr createprob (CPXENVptr env) { 
	int status;
	CPXLPptr lp;
	lp = CPXcreateprob (env, &status, "test");
	if (lp == NULL) {
		printf ("Failed to create LP.\n");
		throw(1);
	}
	return lp;
}

/* load data into a problem instance */
void copylp(CPXENVptr env, CPXLPptr lp, int ncols, int nrows, int minOrMax, double *objective,
	double *rhs, const char *sense, const int *matbeg, const int *matcnt, const int *matind, const double *matval, double *lb, double *ub) {
	int status;
	status = CPXcopylp(env, lp, ncols, nrows, minOrMax, objective, rhs, sense,
		matbeg, matcnt, matind, matval, lb, ub, NULL);
  if (status) {
		printf ("Failed to load problem.\n");
		throw(1); 
	}
}

/* output the problem to a file in MPS format */
void writeprob(CPXENVptr env, CPXLPptr lp) {
	int status;
	/* write to file */
	status = CPXwriteprob(env, lp, "test.mps", "MPS");
  if (status) {
		printf ("Failed to write problem.\n");
		return; 
	}
}

/* solve problem (optimize) */
void lpopt(CPXENVptr env, CPXLPptr lp) {
	int status;
	status = CPXlpopt(env, lp);
	if (status) {
		printf ("Failed to optimize LP.\n");
		throw(1);
	}
}

/* get current solution */
void getsolution(CPXENVptr env, CPXLPptr lp, double *objval, double x[]) {
	int status, solstat;
	status = CPXsolution (env, lp, &solstat, objval, x, NULL, NULL, NULL);
	/* todo: do something with solstat */
	if (status) {
		printf ("Failed to obtain solution.\n");
		writeprob(env,lp);
		throw(1);
	}
#ifdef OUTPUT
	printf ("\nSolution status = %d\n", solstat);
	printf ("Solution value  = %f\n", *objval);
#endif
}

/* get current basis */
void getbase(CPXENVptr env, CPXLPptr lp, int cstat[], int rstat[]) {
	int status;
	status = CPXgetbase(env, lp, cstat, rstat);
	if (status) {
		printf ("Failed to obtain base.\n");
		throw(1);
	}
}

/* reload basis */
void copybase(CPXENVptr env, CPXLPptr lp, const int cstat[], const int rstat[]) {
	int status;
	status = CPXcopybase (env, lp, cstat, rstat);
	if (status) {
		printf ("Failed to reload basis.\n");
		throw(1);
	}
}

/* change bounds */
void chgbds (CPXENVptr env, CPXLPptr lp, int nbds, int vwIndices[], char upper[], double newUppBnd[]) {
	int status;
	status = CPXchgbds (env, lp, nbds, vwIndices, upper, newUppBnd);
	if (status) {
		printf ("Failed to change bounds.\n");
		throw(1);
	}
}

/* change objective */
void chgobj (CPXENVptr env, CPXLPptr lp, int nInds, int chgObjInd[], const double chgObjVal[]) {
	int status;
	status = CPXchgobj (env, lp, nInds, chgObjInd, chgObjVal);
	if (status) {
		printf ("Failed to update objective.\n");
		throw(1);
	}
}
/* change rhs */
void chgrhs (CPXENVptr env, CPXLPptr lp, int nInds, int chgObjInd[], const double chgObjVal[]) {
  int status;
  status = CPXchgrhs (env, lp, nInds, chgObjInd, chgObjVal);
  if (status) {
    printf ("Failed to update rhs.\n");
    throw(1);
  }
}
