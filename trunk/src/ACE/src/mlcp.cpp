/************************************************************************
  			mlcp.cpp 
**************************************************************************/
#include "mlcp.h"
#include "SimpleVector.h"
/*
{M11 M12} Z1  a  W1
{       }*  +  =
{M21 M22} Z2  b  0
*/
mlcp::mlcp(unsigned int Dlcp,unsigned int Dlin,int solverType){

  mDlcp = Dlcp;
  mDlin = Dlin;
  ACE_CHECK_IERROR(mDlcp+mDlin>0,"mlcp::mlcp dim null");
//   mA=0;
//   mB=0;
//   mC=0;
//   mD=0;
//   ma=0;
//   mb=0;
  mu=0;
  mv=0;
  mw=0;

  mW1=0;
  mZ1=0;
  mZ2=0;
  mQ1=0;
  mQ2=0;
  mM11=0;
  mM12=0;
  mM21=0;
  mM22=0;
  mMd =0;
  mSolverType=solverType;
  mTryOnlyGuess = false;
  mTringGuess = false;
  mUseGuess = true;
  //  if (mSolverType !=ACE_SOLVER_ENUM){
    mTryOnlyGuess = true;
    
    mMd=(double *)calloc((mDlcp+mDlin)*(mDlcp+mDlin),sizeof(double));

//     mA=(double *)calloc(mDlin*mDlin,sizeof(double));
//     mC=(double *)calloc(mDlin*mDlcp,sizeof(double));
//     ma=(double *)calloc(mDlin,sizeof(double));
    mu=(double *)calloc(mDlin+mDlcp,sizeof(double));
    mw=(double *)calloc(mDlcp+mDlin,sizeof(double));
//     if (mDlcp){
//       mb=(double *)calloc(mDlcp,sizeof(double));
//       mD=(double *)calloc(mDlcp*mDlin,sizeof(double));
//       mB=(double *)calloc(mDlcp*mDlcp,sizeof(double));
//     }
     mv=mu+mDlin;
    //  }
   
  mW1Z1=0;
  mM=0;
  mCurEnum=8388600;
  mCmp=0;
  mMaxEnum = (unsigned long) powl(2,(long)mDlcp);
  mPourCent=0;

  mM = new aceMatrix(Dlcp+Dlin,Dlcp+Dlin);
  mIwork=0;
  mDwork=0;

  mTryM = false;

  if (mDlcp){
    mW1Z1 = (int*)calloc(mDlcp,sizeof(int));
    //mW1Z1prev = (int*)calloc(mDlcp,sizeof(int));
    mW1 = new aceVector(mDlcp);
    mZ1 = new aceVector(mDlcp);
    mQ1=new aceVector(mDlcp);
    mM11 = new aceMatrix(mDlcp,mDlcp);
  }
   
  if (mDlin){
    mM22 = new aceMatrix(mDlin,mDlin);
    mZ2 = new aceVector(mDlin);
    mQ2 = new aceVector(mDlin);     
  }
   
  if (mDlin && mDlcp){
    mM21 = new aceMatrix(mDlin,mDlcp);
    mM12 = new aceMatrix(mDlcp,mDlin);
  }
  //   mGuess.clear();
  mGuess.push_back(0);

  dim.resize(2);
  start.resize(4);
  end.resize(2);
  dim[0]=mDlin;
  dim[1]=1;
  start[0]=mDlcp;
  start[1]=0;
  start[2]=0;
  start[3]=0;


}
bool mlcp::solveLinearSystem(){
  try{
    mM22->PLUForwardBackwardInPlace(*mQ2);
    for (int lin=0;lin<mDlin;lin++)
      mZ2->setValue(lin,mQ2->getValue(lin));
    return true;
  }
  catch(SiconosException e)
    {
      std::cout << e.report() << endl;
      ACE_ERROR("mlcp::solveLinearSystem linear system no solution\n");
    }
  catch(...)
    {
      std::cout << "Exception caught." << endl;
      ACE_ERROR("mlcp::solveLinearSystem linear system no solution\n");
    }
  return false;
}
/*
 *input mQ1 and mQ2  =========> mProblem.q
 *
 * call mlcp_drivers
 *
 *
 *output mu and mw ==========> mZ1, mZ2, mW1
 */
bool mlcp::solveWithNumerics(){
  int n = mDlin;
  int m = mDlcp;
  bool res=false;
  int info;
  mQ2->VectorToFortran(mProblem.q);
  mQ1->VectorToFortran(mProblem.q+mDlin);
  /*  for (int i=0;i<mDlin+mDlcp;i++)
      mProblem.q[i]=-mProblem.q[i];*/

  ACE_times[ACE_TIMER_SOLVE_NUMERICS].start();
//   printf("\nbegin\n");
//   displayMLCP(&mProblem);
//   printf("\nend\n");
//    printInPut();

  info=mlcp_driver( &mProblem, mu , mw , &mOptions,&mNumericsOptions);
  ACE_times[ACE_TIMER_SOLVE_NUMERICS].stop();
  
  mZ1->FortranToVector(mv);
  mZ2->FortranToVector(mu);
  mW1->FortranToVector(mw+mDlin);
  return ! info;

}
void mlcp::stopSolver(){
  printf("number of case direct solver failed %d \n", mOptions.iparam[7]);
  free(mProblem.M);
  free(mProblem.q);
  free(mOptions.iparam);
  free(mOptions.dparam);
  if (mDlcp){
    mlcp_driver_reset(&mProblem,&mOptions);
    if (mIwork)
      free(mIwork);
    if (mDwork)
      free(mDwork);
    mIwork=0;
    mDwork=0;
  }
}
/*
 *
 *input M11,M12,M21 and M22 =======> mM =====> problem.M.matrix0 (ie mMd)
 *
 *
 */
void mlcp::update(){
  mM->setBlock(0,0,*mM22);
  if (mDlcp){
    mM->setBlock(0,mDlin,*mM21);
    mM->setBlock(mDlin,0,*mM12);
    mM->setBlock(mDlin,mDlin,*mM11);
  }
  //  printInPutABCDab();
  //  cout<<*mM;
  mM->MatrixToFortran(mMd);

}
bool mlcp::initSolver(){
  int n = (int)mDlin;
  int m = (int)mDlcp;


//   mM22->MatrixToFortran(mA);  
//   if (mDlcp){
//     mM21->MatrixToFortran(mC);
//     mM12->MatrixToFortran(mD);
//     mM11->MatrixToFortran(mB);
//   }
  update();
  if (ACE_MUET_LEVEL)
    mNumericsOptions.verboseMode=0;
  else
    mNumericsOptions.verboseMode=1;
  //build the numerics mProblem:
  mProblem.n=n;
  mProblem.m=m;
  mProblem.q=(double*)malloc((n+m)*sizeof(double));
  mProblem.problemType=0;
  mProblem.M= (NumericsMatrix*)malloc(sizeof(NumericsMatrix));
  mProblem.M->storageType=0;
  mProblem.M->size0 = n+m;
  mProblem.M->size1 = n+m;
  mProblem.M->matrix0 = mMd;
  mProblem.M->matrix1 = 0;
  
  mOptions.isSet = 1;
  mOptions.iSize = 9;
  mOptions.iparam = (int*) malloc(10*sizeof(int));
  mOptions.dSize = 9;
  mOptions.dparam = (double*) malloc(10*sizeof(double));
  mOptions.filterOn = 0;
  mOptions.iparam[5]=15;/*Number of registered configurations*/
  mOptions.iparam[7]=0;
  //If adaptive time stepping, then mlcp formulation will change
  if(ACE_WITH_ADAPTATIVE_TIME_STEPPING)
    mOptions.iparam[8]=1;
  else
    mOptions.iparam[8]=0;
  
  if (ACE_SOLVER_TYPE == ACE_SOLVER_SIMPLEX){
    strcpy(mOptions.solverName,"DIRECT_SIMPLEX");
    mOptions.iparam[0]=1000000;
    mOptions.iparam[1]= mNumericsOptions.verboseMode;/*VERBOSE*/
    mOptions.iparam[6]= mNumericsOptions.verboseMode;/*VERBOSE*/
    mOptions.dparam[0]=1e-12;
    mOptions.dparam[1]=1e-12;
    mOptions.dparam[2]=1e-9;

  }else if (ACE_SOLVER_TYPE == ACE_SOLVER_FB){
    strcpy(mOptions.solverName,"DIRECT_FB");

    mOptions.iparam[0]=20000;
    mOptions.iparam[1]= mNumericsOptions.verboseMode;/*VERBOSE*/
    mOptions.iparam[6]= mNumericsOptions.verboseMode;/*VERBOSE*/
    mOptions.iparam[8]=0;/*update prb*/
    mOptions.dparam[0]=1e-12;
    mOptions.dparam[1]=1e-12;
    mOptions.dparam[2]=1e-12;
    mOptions.dparam[6]=1e-12;

  }else if (ACE_SOLVER_TYPE == ACE_SOLVER_PATH){
    strcpy(mOptions.solverName,"DIRECT_PATH");
    mOptions.iparam[0]= mNumericsOptions.verboseMode;/*VERBOSE*/
    mOptions.iparam[6]= mNumericsOptions.verboseMode;/*VERBOSE*/
    mOptions.dparam[0]=1e-12;
    mOptions.dparam[5]=1e-12;
    mOptions.dparam[6]=1e-12;
  }else{
    strcpy(mOptions.solverName,"DIRECT_ENUM");
    mOptions.iparam[0]= mNumericsOptions.verboseMode;/*VERBOSE*/
    mOptions.iparam[6]= mNumericsOptions.verboseMode;/*VERBOSE*/
    mOptions.dSize=6;
    mOptions.dparam[0]=1e-12;
    mOptions.dparam[5]=1e-12;
    mOptions.dparam[6]=1e-12;

  }

  if (mDlcp){
    int nbInts = mlcp_driver_get_iwork(&mProblem,&mOptions);
    int nbDoubles = mlcp_driver_get_dwork(&mProblem,&mOptions);
    mIwork=(int*)malloc( nbInts*sizeof(int));
    mOptions.iWork = mIwork;
    mDwork = (double*)malloc( nbDoubles*sizeof(double));
    mOptions.dWork = mDwork;
    mlcp_driver_init(&mProblem,&mOptions);
  }

  return true;
}
bool mlcp::solve(){
  bool res =false;
  if (!mDlcp){
    res = solveLinearSystem();
  }else{
    res = solveWithNumerics();
    //    printOutPut();
  }
  return res;
}



mlcp::~mlcp(){
  if (mW1Z1)
    free(mW1Z1);
  if (mM)
    delete mM;

    
  if (mW1)
    delete mW1;
  if (mZ1)
    delete mZ1;
  if (mZ2)
    delete mZ2;
  if (mQ1)
    delete mQ1;
  if (mQ2)
    delete mQ2;
  if (mM11)
    delete mM11;
  if (mM12)
    delete mM12 ;
  if (mM21)
    delete mM21;
  if (mM22)
    delete mM22;
  mGuess.clear();
  if (mMd)
    free(mMd);
//   if (mA)
//     free(mA);
//   if(mB)
//     free(mB);
//   if(mC)
//     free(mC);
//   if (mD)
//     free(mD);
//   if(ma)
//     free(ma);
//   if(mb)
//     free(mb);
   if(mu)
     free(mu);
  
////   if(mv)
////     free(mv);
  if(mw)
    free(mw);


}


void mlcp::printInPutABCDab(ostream& os)
{
//   if (ACE_MUET_LEVEL == ACE_MUET)
//     return;
  os<<"mlcp print input"<<endl;
  os<<"dim lcp:\t"<<mDlcp<<"\tlin\t"<<mDlin;
  os<<"A\n";
  if (mM22)
    os<<(*mM22);
  os<<"B\n";
  if (mM11)
    os<<(*mM11);
  os<<"C\n";
  if (mM21)
    os<<(*mM21);
  os<<"D\n";
  if (mM12)
    os<<(*mM12);
  os<<"a\n";
  if (mQ2)
    os<<(*mQ2);
  os<<"b\n";
  if (mQ1)
    os<<(*mQ1);
  
}
void mlcp::printInPut(ostream& os)
{
//   if (ACE_MUET_LEVEL == ACE_MUET)
//     return;
  displayMLCP(&mProblem);
}
void mlcp::printOutPut(ostream& os){
   if (ACE_MUET_LEVEL == ACE_MUET )
     return;
  os<<"mlcp print output"<<endl;
  os<<"v ie Z1:\n";
  if (mZ1)
    os<<(*mZ1);
  os<<"u ie Z2:\n";
  if (mZ2)
    os<<(*mZ2);
  os<<"w ie W1:\n";
  if (mW1)
    os<<(*mW1);
}
