/* Siconos-Numerics version 2.1.1, Copyright INRIA 2005-2007.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.	
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr 
*/
/*!\file mlcp_simplex.cpp
 
  This subroutine allows the resolution of MLCP (Mixed Linear Complementary Problem).\n
  Try \f$(u,v,w)\f$ such that:\n
  \f$
   \left\lbrace
    \begin{array}{l}
    A u + Cv +a =0\\
    D u + Bv +b = w
    0 \le v \perp  w \ge 0\\
    \end{array}
   \right.
  \f$
  where  A is an (\f$ n \times n\f$ ) matrix, B is an (\f$ m \times m\f$ ) matrix,  C is an (\f$ n \times m\f$ ) matrix, D is an (\f$ m \times n\f$ ) matrix,    a and u is an (\f$ n \f$ ) vectors b,v and w is an (\f$ m \f$ ) vectors.
*/
/*!\fn   void mlcp_simplex( int *nn , int* mm, double *A , double *B , double *C , double *D , double *a  double *b, double *u, double *v, double *w , int *info ,	 int *iparamMLCP , double *dparamMLCP );


  mlcp_pgs (Projected Gauss-Seidel) is a basic Projected Gauss-Seidel solver for MLCP.\n
 
  \param A            On enter, a (\f$n \times n\f$)-vector of doubles which contains the components of the "A" MLCP matrix with a Fortran storage.
  \param B            On enter, a (\f$m \times m\f$)-vector of doubles which contains the components of the "B" MLCP matrix with a Fortran storage.
  \param C            On enter, a (\f$n \times m\f$)-vector of doubles which contains the components of the "C" MLCP matrix with a Fortran storage.
  \param D            On enter, a (\f$m \times n\f$)-vector of doubles which contains the components of the "D" MLCP matrix with a Fortran storage.
  \param a            On enter, a n-vector of doubles which contains the components of the constant right hand side vector.
  \param b            On enter, a m-vector of doubles which contains the components of the constant right hand side vector.
  \param nn           On enter, an integer which represents one the dimension of the MLCP problem.
  \param mm           On enter, an integer which represents one the dimension of the MLCP problem.
  \n \n
  \param u            On return, a n-vector of doubles which contains the solution of the problem.
  \param v            On return, a m-vector of doubles which contains the solution of the problem.
  \param w            On return, a m-vector of doubles which contains the complementary solution of the problem.
  \param info    On return, an integer which returns the termination value:\n
                 0 : convergence\n
                 1 : iter = itermax\n
                 2 : negative diagonal term
 
  \param iparamMLCP  On enter/return a vector of integers:\n
                - iparamMLCP[0] = itermax On enter, the maximum number of iterations allowed.
                - iparamMLCP[1] = verbose  On enter, the output log identifiant:\n
                        0 : no output\n
                        >0: active screen output\n
                - iparamMLCP[2] = it_end  On enter, the number of iterations performed by the algorithm.
 
  \param dparamMLCP  On enter/return a vector of doubles:\n
                - dparamMLCP[0] = tol     On enter, the tolerance required.
                - dparamMLCP[1] = omega   On enter, the relaxation parameter (not yet available).
                - dparamMLCP[2] = res     On return, the final error value.
 
  \author F. Cadoux O. Bonnefon
 
 */
#include "mlcp_simplex.h"
#include <ilcplex/cplex.h>
#include <cstdlib>
#include <set>
#include <cmath>
#include <ctime>
#include "ace.h"
static int n=0;
static int m=0;
static int nrows = 0;
static int ncols = 0;
static int nnz=0;
static  double tCum=0.0;

struct node {
	int *cstat;		/* [ncols] current col basis (for restart) */
	int *rstat;		/* [nrows] current row basis (for restart) */
	char *choice;			/* [m] for choices made so far: n = 'no choice', v = 'v->0', w = 'w->0' */
	void print() const {
	  if (ACE_MUET_LEVEL >= ACE_MUET)
	    return;

		int i;
		printf("Choice= ");
		for (i=0; i<m; i++) printf("%c", choice[i]);
		printf("\n");
	};
};
void buildNode(node* _node){
  _node->cstat = (int*)calloc(ncols,sizeof(int));
  _node->rstat = (int*)calloc(nrows,sizeof(int));
  _node->choice = (char*)calloc(m,sizeof(char));
}
void deleteNode(node* _node){
  free(_node->cstat);
  free(_node->rstat);
  free(_node->choice);
}
/* simple comparator for ordering the nodes.
   Could probably be smarter... */
struct compNodes {
  bool operator()(const node& n1, const node& n2) const {
    return true; /* find a better idea ! */
  }
};

typedef std::set<node,compNodes> mySet;
/* auxiliary functions used to call cplex */
#include "auxFun.h"

void print_double_array(char * name,double *t, int nn){
  if (ACE_MUET_LEVEL >= ACE_MUET)
    return;
  int sizeLine = 10;
  int i = 0;
  printf("%s[%d] = (",name,nn);
  for (i=0;i<nn;i++){
    if (i%sizeLine==sizeLine-1) printf("\n");
    printf("%g ,",t[i]);
  }
  printf(")\n");
}
void print_int_array(char * name,int *t, int nn){
  if (ACE_MUET_LEVEL >= ACE_MUET)
    return;
  int sizeLine = 10;
  int i = 0;
  printf("%s[%d] = (",name,nn);
  for (i=0;i<nn;i++){
    if (i%sizeLine==sizeLine-1) printf("\n");
    printf("%d ,",t[i]);
  }
  printf(")\n");
}
void print_sparse_matrix(int *matbeg,int *matcnt,int *matind,double *matval,double *objective,double *rhs,char *sense,double *lb,double *ub){
  if (ACE_MUET_LEVEL >= ACE_MUET)
    return;
  print_int_array("matbeg",matbeg,ncols);
  print_int_array("matcnt",matcnt,ncols);
  print_int_array("matind",matind,nnz);
  print_double_array("matval",matval,nnz);
  print_double_array("rhs",rhs,nrows);
  print_double_array("lb",lb,ncols);
  print_double_array("ub",ub,ncols);
}


/* function that processes each child node and either prunes it, gets a solution, or adds it to the tree */
int solvesubproblem(CPXENVptr env, CPXLPptr lp, mySet::iterator& currentNode, int branchingIndex, char v_or_w, mySet& tree, double x[], const double tolVar, const double tolComp) {
	int i, solstat;
	double objval;

	/* reload basis */
	copybase (env, lp, currentNode->cstat, currentNode->rstat);
	/* minimize v_i */
	lpopt(env, lp);
	/* retrieve sol */
	getsolution (env, lp, &objval, x);

	if (std::fabs(objval)<tolVar) { /* minimization succeeded */
		double infeasComp = 0;
		for (i=0; i<m; i++) {
			infeasComp+=x[n+i]*x[n+m+i];
		}
		if (infeasComp<tolComp) { /* we have a solution !!! */
		  if (ACE_MUET_LEVEL < ACE_MUET){
			printf("\n----- Getting a solution:\n");
			for (i=0; i<ncols; i++) {
				printf("[%2.i] %10.7f", i, x[i]);
				if (i%5==4) printf("\n");
			}
			printf("\n----- End of solution.\n");
		  }
			return 1;
		} else { /* complementarity doesn't hold, keep going */
			node child;
			buildNode(&child);
			getbase(env, lp, child.cstat, child.rstat);
			memcpy(child.choice, currentNode->choice, m*sizeof(char));
			child.choice[branchingIndex]=v_or_w;
			tree.insert(child);
		}
	} else {
		/*printf("LP pruned node (objval=%f).\n", objval);*/
	}
	return 0;
}



int mlcp_simplex( int *nn , int* mm, double *A , double *B , double *C , double *D , double *a, double *b, double *u, double *v, double *w , int *info ,	 int *iparamMLCP , double *dparamMLCP  )
{ 
  /* Parameters */
  const double tolVar = 1e-3;			/* tolerance to consider that a var is null */
  const double tolComp = 1e-3; 		/* tolerance to consider that complementarity holds */
  const int nIterMax = 1000000;		/* max number of nodes to consider in tree search */
  const int logFrequency = 10000;	/* print log every ... iterations */

  /* Variables declaration */
  int curbeg;
  int *matbeg=0;
  int *matcnt=0;
  int *matind;
  double *matval;
  double *objective;
  double *rhs;
  char *sense;
  double *lb;
  double *ub;

  int linNumber;
  int colNumber;

  int i,j;
  
  /* Variables declaration */
  mySet tree;
  int  solstat;
  double objval;
  double *x;//[ncols];
  CPXENVptr env = NULL;
  CPXLPptr lp  = NULL;

  n = *nn;
  m = *mm;
  nnz=m;
  nrows = n+m;
  ncols = n+m+m;
  x=(double*) calloc(ncols,sizeof(double));


  /*SPARSE matrix, count non nul values*/
  /*i line*/
  /*j col*/
  linNumber = n;
  colNumber = n;
  for (i=0;i<linNumber;i++)
    for (j=0;j<colNumber;j++){
      if (A[j*linNumber+i]!=0.0)
	nnz++;
    }
  linNumber = m;
  colNumber = m;
  for (i=0;i<linNumber;i++)
    for (j=0;j<colNumber;j++){
      if (B[j*linNumber+i]!=0.0)
	nnz++;
    }
  linNumber = n;
  colNumber = m;
  for (i=0;i<linNumber;i++)
    for (j=0;j<colNumber;j++){
      if (C[j*linNumber+i]!=0.0)
	nnz++;
    }
  linNumber = m;
  colNumber = n;
  for (i=0;i<linNumber;i++)
    for (j=0;j<colNumber;j++){
      if (D[j*linNumber+i]!=0.0)
	nnz++;
    }
  /*nnz contains the number of non nul values*/
  if (!nnz || !nrows || !ncols)
    return 0;

  matbeg = (int*)calloc(ncols+1,sizeof(int));/*The last value is not used.*/
  matcnt = (int*)calloc(ncols,sizeof(int));
  objective = (double*)calloc(ncols,sizeof(double));
  sense = (char*)calloc(nrows,sizeof(char));
  lb = (double*)calloc(ncols,sizeof(double));
  ub = (double*)calloc(ncols,sizeof(double));
  matind= (int*)calloc(nnz,sizeof(int));
  matval= (double*)calloc(nnz,sizeof(double));
  rhs=(double*)calloc(nrows,sizeof(double));

  /*INITIAL */
  for (i=0; i<ncols; i++) objective[i]=0;
  
  for (i=0; i<n; i++) { /* first n vars = u, no bounds */
    lb[i] = -CPX_INFBOUND;
    ub[i] =  CPX_INFBOUND;
  }
  for (i=n; i<n+2*m; i++) {
    lb[i] = 0;
    ub[i] = CPX_INFBOUND;
  }
  for (i=n; i<nrows; i++) {
    sense[i] = 'E';
  }
  /*fill sprase matrix*/
  curbeg=0;
  matbeg[0]=0;
  colNumber = n;
  /*Value from A and D*/
  for (j=0;j<colNumber;j++){
    linNumber=n;
    /*Values from A*/
    for (i=0;i<linNumber;i++)
      if(A[j*linNumber+i]!=0.0){
	matval[curbeg] = A[j*linNumber+i];
	matind[curbeg]=i;
	curbeg++;
      }
    linNumber=m;
    /*Values from D*/
    for (i=0;i<linNumber;i++)
      if(D[j*linNumber+i]!=0.0){
	matval[curbeg] = D[j*linNumber+i];
	matind[curbeg]=i+n;
	curbeg++;
      }
    matbeg[j+1]=curbeg;
    matcnt[j]=matbeg[j+1]-matbeg[j];
  }

  colNumber = m;
  /*Value from C and B*/
  for (j=0;j<colNumber;j++){
    linNumber=n;
    /*Values from C*/
    for (i=0;i<linNumber;i++)
      if(C[j*linNumber+i]!=0.0){
	matval[curbeg] = C[j*linNumber+i];
	matind[curbeg]=i;

	curbeg++;
      }
    linNumber=m;
    /*Values from B*/
    for (i=0;i<linNumber;i++)
      if(B[j*linNumber+i]!=0.0){
	matval[curbeg] = B[j*linNumber+i];
	matind[curbeg]=i+n;
	curbeg++;
      }
    matbeg[j+n+1]=curbeg;
    matcnt[j+n]=matbeg[j+n+1]-matbeg[j+n];
  }
  matcnt[n+m]=1;
  for (i=1;i< m;i++){
    matbeg[n+m+i]=curbeg+i;
    matcnt[n+m+i]=1;
  }
  for (i=0;i< m;i++){
    matval[nnz-i-1]=-1;
    matind[nnz-i-1]=n+m-1-i;
  }

  /*fill rhs*/
  for (i=0;i<n;i++)
    rhs[i]=-a[i];
  for (i=0;i<m;i++)
    rhs[n+i]=-b[i];

  //print_sparse_matrix( matbeg,matcnt,matind,matval,objective,rhs,sense,lb,ub);

  /*READY TO SOLVE*/
  	/* to compute time spent */
  time_t tStart = time(NULL);
	
  /* open cplex, create prob, load data, write to file */
  env = openCPLEX();
  lp = createprob(env);
  copylp(env, lp, ncols, nrows, CPX_MIN, objective, rhs, sense, matbeg, matcnt, matind, matval, lb, ub);
  writeprob(env, lp);
	
  /* optimize with null objective to find feasible point */
  lpopt(env, lp);
  getsolution (env, lp, &objval, x);
	
  /* create root node */
  node root;
  buildNode(&root);
  getbase(env, lp, root.cstat, root.rstat);
  for (i=0; i<m; i++) root.choice[i]='n';
  tree.insert(root);
	
  /* variables for tree search */	
  int lastBranchingIndex=0;		/* last index we minimized */
  int branchingIndex=0;				/* new index */
  int nFixed;										/* number of fixed variables (0...m) in current node */
  char *upper;							/* lower or upper ? (for changing bounds) */
  upper = (char*) calloc(2*m,sizeof(char));
  memset(upper,'U',2*m*sizeof(char));
  int *vwIndices;						/* indices of variables v,w */
  vwIndices = (int*)calloc( 2*m,sizeof(int));
  for (i=0; i<2*m; i++) vwIndices[i]=n+i;
  double *newUppBnd;				/* new upper bounds */
  newUppBnd = (double*)calloc( 2*m,sizeof(double));
  int nIter=0;									/* iteration count */
  int chgObjInd[2];							/* indices in obj to change */
  const double chgObjVal[2] = {0., 1.};
  int res=0;
  /* new obj val */

  while(!tree.empty() && nIter < nIterMax) {
    nIter++;
    mySet::iterator currentNode = tree.begin();
    //currentNode->print();//DEBUG
    /* update LP prob and choose index to branch on */
    nFixed=0;
    branchingIndex = -1;
    for (i=0; i<m; i++) {
      switch (currentNode->choice[i]) {
      case 'n':
	// 1st branching strategy: normal order
	//if (branchingIndex==-1)	branchingIndex = i;
					
	// 2nd branching strategy : reverse order
	//branchingIndex = i; /* last index */
					
	// 3eme idee : random (find some better way to randomly decide...)
	if (branchingIndex==-1) {
	  branchingIndex = i;
	} else {
	  branchingIndex = (2*rand() < RAND_MAX ? branchingIndex : i);
	}
					
	newUppBnd[i] = CPX_INFBOUND;
	newUppBnd[m+i] = CPX_INFBOUND;
	break;
      case 'v':
	newUppBnd[i] = 0;
	newUppBnd[m+i] = CPX_INFBOUND;
	nFixed++;
	break;
      case 'w':
	newUppBnd[i] = CPX_INFBOUND;
	newUppBnd[m+i] = 0;
	nFixed++;
	break;
      default:
	if (ACE_MUET_LEVEL < ACE_MUET)
	  printf("Choice should be either n,v or w (here: %c)\n", currentNode->choice[i]);
	return -1;
      }
    }
    //printf("bI=%i, TS=%i\n", branchingIndex, tree.size());//DEBUG
    if (branchingIndex == -1) {
      /* This situation should be rather rare :
	 all ones out of two variables were fixed to approx. 0,
	 but complementarity still doesn't hold up to tolComp
	 (though it cannot be much larger with all the zero fixed vars.)
	 The current node is probably a solution, up to roundoff errors...	*/
      if (ACE_MUET_LEVEL < ACE_MUET){
	printf("Node without any free variable was found. It might be a solution, or close.\n");
	for (i=0; i<m; i++) printf("%c", currentNode->choice[i]);
	printf("\n");
      }
      deleteNode((node*)&(*currentNode));
      tree.erase(currentNode);
      break;
    }
    chgbds (env, lp, 2*m, vwIndices, upper, newUppBnd);
				
    /* reset obj and solve left-child subproblem */
    chgObjInd[0]=lastBranchingIndex;
    chgObjInd[1]=n+branchingIndex;
    lastBranchingIndex=chgObjInd[1];
    chgobj (env, lp, 2, chgObjInd, chgObjVal);
    if (solvesubproblem(env, lp, currentNode, branchingIndex, 'v', tree, x, tolVar, tolComp)){
      for( i = 0 ; i < m ; ++i ) {v[i] = x[i+n]; w[i]=x[n+m+i];}
      for( i = 0 ; i < n ; ++i ) {u[i] = x[i]; }
      break;
    }

    /* reset obj and solve right-child subproblem */
    chgObjInd[0]=lastBranchingIndex;
    chgObjInd[1]=n+m+branchingIndex;
    lastBranchingIndex=chgObjInd[1];
    chgobj (env, lp, 2, chgObjInd, chgObjVal);
    if (solvesubproblem(env, lp, currentNode, branchingIndex, 'w', tree, x, tolVar, tolComp)){
      for( i = 0 ; i < m ; ++i ) {v[i] = x[i+n]; w[i]=x[n+m+i];}
      for( i = 0 ; i < n ; ++i ) {u[i] = x[i]; }
      res =1;
      break;
    }

    /* remove current node */
    deleteNode((node*)&(*currentNode));
    tree.erase(currentNode);
		
    /* disp some log */
    if (nIter%logFrequency==0 && ACE_MUET_LEVEL < ACE_MUET) printf("iter %i : tree size is now %i \n", nIter, tree.size());
  }



  if (upper) free(upper);
  if (vwIndices) free(vwIndices);
  if (newUppBnd) free(newUppBnd);
  if (x) free(x);
  
  /* display results */
  time_t tEnd = time(NULL);
  tCum+=difftime(tEnd, tStart);
  if (ACE_MUET_LEVEL < ACE_MUET)
    printf("%i iterations taken (%f sec.), tree size is %i\n", nIter,tCum , tree.size());
  while (!tree.empty()){
    mySet::iterator currentNode = tree.begin();
    deleteNode((node*)&(*currentNode));
    tree.erase(currentNode);
    res =1;
  }

  free(matbeg);
  free(matcnt);
  free(objective);
  free(sense);
  free(lb);
  free(ub);
  free(matval);
  free(matind);
  free(rhs);
  return res;
}
