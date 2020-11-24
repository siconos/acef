
#ifndef EXTERN_MLCP_SIMPLEX
#define EXTERN_MLCP_SIMPLEX

/**/
int extern_mlcp_simplex( double *a, double *b, double *u, double *v, double *w , int *info ,	 int *iparamMLCP , double *dparamMLCP  );
unsigned long extern_getConfigLCP();
void extern_mlcp_simplex_stop();
void extern_mlcp_simplex_init(int *nn , int* mm, double *A , double *B , double *C , double *D);
void extern_mlcp_simplex_init_with_M(int *nn , int* mm, double *M );
void extern_setVerbose(int v);

#endif //EXTERN_MLCP_SIMPLEX
