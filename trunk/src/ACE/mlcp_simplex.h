
#ifndef MLCP_SIMPLEX
#define MLCP_SIMPLEX

/**/
int mlcp_simplex( double *a, double *b, double *u, double *v, double *w , int *info ,	 int *iparamMLCP , double *dparamMLCP  );
unsigned long getConfigLCP();
void mlcp_simplex_stop();
void mlcp_simplex_init(int *nn , int* mm, double *A , double *B , double *C , double *D);

#endif //MLCP_SIMPLEX
