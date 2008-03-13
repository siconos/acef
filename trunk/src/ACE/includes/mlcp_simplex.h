
#ifndef EXTERN_MLCP_SIMPLEX
#define EXTERN_MLCP_SIMPLEX

/**/
extern "C"int extern_mlcp_simplex( double *a, double *b, double *u, double *v, double *w , int *info ,	 int *iparamMLCP , double *dparamMLCP  );
extern "C" unsigned long extern_getConfigLCP();
extern "C" void extern_mlcp_simplex_stop();
extern "C" void extern_mlcp_simplex_init(int *nn , int* mm, double *A , double *B , double *C , double *D);
extern "C" void extern_mlcp_simplex_init_with_M(int *nn , int* mm, double *M );
extern "C" void extern_setVerbose(int v);

#endif //EXTERN_MLCP_SIMPLEX
