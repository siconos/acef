/************************************************************************
  			mlcp.cpp 
**************************************************************************/
#include "mlcp.h"
#include "SimpleVector.h"
#include "mlcp_simplex.h"
/*
{M11 M12} Z1  a  W1
{       }*  =  +
{M21 M22} Z2  b  0
*/
mlcp::mlcp(unsigned int Dlcp,unsigned int Dlin,int solverType){

  mDlcp = Dlcp;
  mDlin = Dlin;
  ACE_CHECK_IERROR(mDlcp+mDlin>0,"mlcp::mlcp dim null");
  mA=0;
  mB=0;
  mC=0;
  mD=0;
  ma=0;
  mb=0;
  mu=0;
  mv=0;
  mw=0;

  mW1=0;
  mZ1=0;
  mZ2=0;
  mQ1=0;
  mQ2=0;
  mQ=0;
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

    mA=(double *)calloc(mDlin*mDlin,sizeof(double));
    mC=(double *)calloc(mDlin*mDlcp,sizeof(double));
    mD=(double *)calloc(mDlcp*mDlin,sizeof(double));
    mB=(double *)calloc(mDlcp*mDlcp,sizeof(double));
    ma=(double *)calloc(mDlin,sizeof(double));
    mb=(double *)calloc(mDlcp,sizeof(double));
    mu=(double *)calloc(mDlin+mDlcp,sizeof(double));
    mv=mu+mDlin;
    mw=(double *)calloc(mDlcp+mDlin,sizeof(double));
    //  }
   
  mW1Z1=0;
  mM=0;
  mCurEnum=8388600;
  mCmp=0;
  mMaxEnum = (unsigned long) powl(2,(long)mDlcp);
  mPourCent=0;

  mM = new aceMatrix(Dlcp+Dlin,Dlcp+Dlin);

  mTryM = false;
  mQ= new aceVector(Dlcp+Dlin);

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
// void mlcp::addGuess(unsigned long l){
//   Itulongs aux = mGuess.begin();
//   bool find = false;
//   while (aux != mGuess.end()){
//     unsigned long laux = *aux;
//     if (*aux == l){
//       find = true;
//       break;
//     }
//     aux++;
//   }
//   if (  !find){
//     mGuess.push_back(l);
//   }
//   cout <<"add guess "<<l<<endl;
// }
// void mlcp::affectW1Z1(unsigned long ll){
//   unsigned long aux = ll;
//   mCase = ll;
//   for (unsigned int i =0; i <mDlcp; i++){
//     mW1Z1[i]=aux & 1;
//     aux = aux >> 1;
//   }

// }
// void mlcp::setCurrentConfig(unsigned long l){
//   (*(mGuess.begin()))=l;
// }
// void mlcp::initGuess(){
//   mItGuess = mGuess.begin();
//   ACE_MESSAGE("mlcp::initEnum : start guess\n");
// }

// bool mlcp::tryGuess(){
//   if (mItGuess == mGuess.end()){
//     ACE_MESSAGE("mlcp::tryGuess : stop guess\n");
//     mTringGuess = false;
//     ACE_times[ACE_TIMER_SOLVE_GUESS].stop();
//     return false;
//   }
//   if (ACE_MUET_LEVEL != ACE_MUET)
//     printf("try guess %d\n",(int)*mItGuess);
//   affectW1Z1(*mItGuess);
//   mItGuess++;
//   return true;
// }
// void mlcp::initEnum(){

//   mPourCent=0;
//   mCmp=0;
//   if (mUseGuess){
//     ACE_times[ACE_TIMER_SOLVE_GUESS].start();
//     mTringGuess = true;
//     initGuess();
//   }else{
//     mTringGuess = false;
//   }
// }

// bool mlcp::nextEnum(){
//   if (mTringGuess && tryGuess())
//     return true;
//   if (mTryOnlyGuess)
//     return false;
//   ACE_times[ACE_TIMER_SOLVE_ENUM].start();
//   if (mCmp == mMaxEnum)
//     return false;
//   if (mCmp > mPourCent*mMaxEnum){
//     mPourCent+=0.001;
//     printf(" %f %d \n",mPourCent,(int) mCurEnum);
//   }
//   if (mCurEnum >= mMaxEnum){
//     mCurEnum=0;
//   }
//   affectW1Z1(mCurEnum);
//   mCurEnum++;
//   //  if (mCurEnum > 8988600)
//   //  mCurEnum = 4000000;
//   mCmp++;
  
//   return true;
// }
bool mlcp::solveWithNumerics(){
  int n = mDlin;
  int m = mDlcp;
  bool res=false;
  int info;
  mQ2->VectorToFortran(mProblem.q);  
  mQ1->VectorToFortran(mProblem.q+mDlin);
  for (int i=0;i<mDlin+mDlcp;i++)
    mProblem.q[i]=-mProblem.q[i];

  ACE_times[ACE_TIMER_SOLVE_PATH].start();
  info=mlcp_driver( &mProblem, mu , mw , &mOptions,&mNumericsOptions);
  ACE_times[ACE_TIMER_SOLVE_PATH].stop();
  
  mZ1->FortranToVector(mv);
  mZ2->FortranToVector(mu);
  mW1->FortranToVector(mw+mDlin);
  return ! info;

}
void mlcp::stopSolver(){
  free(mProblem.M);
  free(mProblem.q);
  free(mOptions.iWork);
  free(mOptions.dWork);
  mlcp_driver_reset(&mProblem,&mOptions);
}
bool mlcp::initSolver(){
  int n = (int)mDlin;
  int m = (int)mDlcp;


  mM22->MatrixToFortran(mA);  
  mM11->MatrixToFortran(mB);  
  mM21->MatrixToFortran(mC);  
  mM12->MatrixToFortran(mD);
  
  mM->setBlock(0,0,*mM22);
  mM->setBlock(0,mDlin,*mM21);
  mM->setBlock(mDlin,0,*mM12);
  mM->setBlock(mDlin,mDlin,*mM11);
  //  printInPutABCDab();
  //  cout<<*mM;
  mM->MatrixToFortran(mMd);
  mNumericsOptions.verboseMode=0;
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
  mOptions.iparam[5]=5;/*Number of registered configurations*/
  
  if (ACE_SOLVER_TYPE == ACE_SOLVER_SIMPLEX){
    strcpy(mOptions.solverName,"DIRECT_SIMPLEX");
    mOptions.iparam[0]=1000000;
    mOptions.iparam[1]=0;/*VERBOSE*/
    mOptions.dparam[0]=1e-12;
    mOptions.dparam[1]=1e-12;
    mOptions.dparam[2]=1e-9;

  }else if (ACE_SOLVER_TYPE == ACE_SOLVER_PATH){
    strcpy(mOptions.solverName,"DIRECT_PATH");
    mOptions.iparam[0]=0;/*VERBOSE*/
    mOptions.iparam[6]=0;/*VERBOSE*/
    mOptions.dparam[0]=1e-12;
    mOptions.dparam[5]=1e-12;
    mOptions.dparam[6]=1e-12;
  }else{
    strcpy(mOptions.solverName,"DIRECT_ENUM");
    mOptions.iparam[0]=0;/*VERBOSE*/
    mOptions.dSize=6;
    mOptions.dparam[0]=1e-12;
    mOptions.dparam[5]=1e-12;
    mOptions.dparam[6]=1e-12;

  }

  int nbInts = mlcp_driver_get_iwork(&mProblem,&mOptions);
  int nbDoubles = mlcp_driver_get_dwork(&mProblem,&mOptions);
  mOptions.iWork = (int*)malloc( nbInts*sizeof(int));
  mOptions.dWork = (double*)malloc( nbDoubles*sizeof(double));
  mlcp_driver_init(&mProblem,&mOptions);

  return true;
}
// bool mlcp::solveWithSimplex(){
//   int n = mDlin;
//   int m = mDlcp;
//   bool res;

//   mQ2->VectorToFortran(ma);  
//   mQ1->VectorToFortran(mb);
//   //  cout<<"**mQ1\n"<<*mQ1;
//   //  cout<<"**mQ2\n"<<*mQ2;
//   int info;
//   int iparamMLCP[2];
//   double dparamMLCP[3];
//   iparamMLCP[0]= 1000000;
//   iparamMLCP[1]=1;
//   dparamMLCP[0]=1e-12;
//   dparamMLCP[1]=1e-12;
//   dparamMLCP[2]=1e-9;
//   ACE_times[ACE_TIMER_SOLVE_SIMPLEX].start();
//   for (int i=0; i<n;i++)
//     ma[i]=-ma[i];
//   for (int i=0; i <m;i++)
//     mb[i]=-mb[i];
//   res= extern_mlcp_simplex(ma, mb, mu, mv, mw , &info , iparamMLCP , dparamMLCP  );
//   mCase = extern_getConfigLCP();
//   addGuess(mCase);
//   ACE_times[ACE_TIMER_SOLVE_SIMPLEX].stop();

//   mZ1->FortranToVector(mv);
//   mZ2->FortranToVector(mu);
//   mW1->FortranToVector(mw);
  
//   return res;
// }
bool mlcp::solve(){
  bool res =false;
  ACE_times[ACE_TIMER_SOLVER].start();

  res = solveWithNumerics();
//   //1) start with guess point
//   ACE_MESSAGE("mlcp::solve : 1 try with guess point.\n");
//   res = solveGuessAndIt();
//   if (res){
//     ACE_times[ACE_TIMER_SOLVER].stop();
//     return res;
//   }
//   //2) start other algo
//   ACE_MESSAGE("mlcp::solve : 2 start other algo\n");
//   switch (mSolverType) {
//   case ACE_SOLVER_SIMPLEX:
//     res= solveWithSimplex();
//     break;
//   case ACE_SOLVER_PATH:
//     res= solveWithPath();
//     break;
//   case ACE_SOLVER_ENUM:
//     ;//res = solveGuessAndIt();
//     break;
//   default:
//     ACE_ERROR("mlcp::solve, bad mSolverType value");
//   }
//   if (!res){
//     ACE_MESSAGE("mlcp::solve, failed.\n");
//   }else{
//     ACE_MESSAGE("mlcp::solve : 3 Add a new guess point.\n");
//     //    addGuess(mZ1);

//   }
  ACE_STOP_SOLVER_TIME();
  return res;
}

// bool mlcp::solveGuessAndIt(){
//   bool onlyOne = true;

//   bool find = false;
//   unsigned int col =0;
//   unsigned int lin =0;
//   if (mDlin==0){//LCP case
//     ACE_ERROR("mlcp::solveGuessAndIt : LCP not implemented");
//   }else if(mDlcp==0){//Linear system case
//     ACE_MESSAGE("mlcp::solveGuessAndIt : Linear system\n");
//     try{
//       //(*mQ2)=-1*(*mQ2);
//       //cout <<(*mM22);
//       //cout <<(*mQ2);
//       mM22->PLUForwardBackwardInPlace(*mQ2);
//       for (lin=0;lin<mDlin;lin++)
// 	mZ2->setValue(lin,mQ2->getValue(lin));
//       return true;
//     }
//     catch(SiconosException e)
//     {
//       std::cout << e.report() << endl;
//       ACE_ERROR("linearSystem::solveGuessAndIt linear system no solution");
//     }
//     catch(...)
//     {
//       std::cout << "Exception caught." << endl;
//       ACE_ERROR("linearSystem::solveGuessAndIt linear system no solution\n");
//     }
//   }
//   if (mTryM){
//     ACE_times[ACE_TIMER_DIRECT].start();
//     mQ->setBlock(0,*mQ1);
//     mQ->setBlock(mDlcp,*mQ2);
//     //scal(-1,*mQ,*mQ);
//     //    *mQ=-1*(*mQ);
//     ACE_times[ACE_TIMER_LU_DIRECT].start();
//     //    cout<<(*mM)<<endl;
//     mM->PLUForwardBackwardInPlace(*mQ);
//     ACE_times[ACE_TIMER_LU_DIRECT].stop();
//     bool check=true;
//     for (lin = 0 ; lin <mDlcp; lin++){
//       //     double aux = mQ->getValue(lin,0);
//       if (mQ->getValue(lin) < - ACE_NULL){
// 	ACE_MESSAGE("mlcp::solveGuessAndIt not in the cone\n");
// 	check=false;
// 	break;//pas dans le cone!
//       }
//     }
//     if(check){
//       fillSolution();
//       ACE_times[ACE_TIMER_DIRECT].stop();
//       return true;
//     }
//     ACE_times[ACE_TIMER_DIRECT].stop();
//   }
  
//   if(ACE_MUET_LEVEL !=10 && onlyOne){
//     onlyOne=false;
//     mM->setBlock(0,0,*mM11);
//     mM->setBlock(0,mDlcp,*mM12);
//     mM->setBlock(mDlcp,0,*mM21);
//     mM->setBlock(mDlcp,mDlcp,*mM22);
//     mQ->setBlock(0,*mQ1);
//     mQ->setBlock(mDlcp,*mQ2);
//     printInPutABCDab();
//   }
//   initEnum();
//   while(nextEnum()){
//     //if(!mTringGuess){
//     mM->setBlock(0,mDlcp,*mM12);
//     mM->setBlock(mDlcp,mDlcp,*mM22);
//     mQ->setBlock(0,*mQ1);
//     mQ->setBlock(mDlcp,*mQ2);
    
//     if(ACE_MUET_LEVEL !=10)
//       printInPut();
//    //build mM11_ and mM21_
//    for (col =0; col<mDlcp; col++){
//      if (mW1Z1[col]==0){
//        //M11
//        for (lin =0; lin < mDlcp; lin++)
// 	 mM->setValue(lin,col,mM11->getValue(lin,col));
       
//        //M21
//        for (lin =0; lin < mDlin; lin++)
// 	 mM->setValue(mDlcp+lin,col,mM21->getValue(lin,col));
       
//      }else{
//        //M11
//        for (lin =0; lin < mDlcp; lin++)
// 	 mM->setValue(lin,col,0);
//        mM->setValue(col,col,-1);
       
//        //M21
//        for (lin =0; lin < mDlin; lin++)
// 	 mM->setValue(mDlcp+lin,col,0);
//      }
//    }//build mM is built
   
//   try{//solve the system
//     if (ACE_MUET_LEVEL != ACE_MUET){
//       ACE_MESSAGE("mlcp:: try: ");
//       for (lin = 0 ; lin<mDlcp;lin++)
// 	cout<<mW1Z1[lin]<<" ";
//       cout<<endl;
//       printInPut();
//     }
//     //
//     //*mQ=-1*(*mQ);
//    ACE_times[ACE_TIMER_SOLVE_LU].start();
//    mM->PLUForwardBackwardInPlace(*mQ);
//    ACE_times[ACE_TIMER_SOLVE_LU].stop();
//    mTryM=true;
//    bool check=true;
//    for (lin = 0 ; lin <mDlcp; lin++){
//      double aux = mQ->getValue(lin);
//      if (mQ->getValue(lin) < - ACE_NULL){
//        ACE_MESSAGE("mlcp::solveGuessAndIt not in the cone\n");
//        check=false;
//        break;//pas dans le cone!
//      }
//    }
//    if(!check)
//      continue;
//    ACE_MESSAGE("mlcp::solveGuessAndIt succes\n");
//    if (ACE_MUET_LEVEL != ACE_MUET)
//      printf("val de step %d \n",(int)mCase);
//    if (mUseGuess && ! mTringGuess)
//      addGuess(mCase);
//    setCurrentConfig(mCase);
//    find = true;
//    break;// gagne.
//   }
//   catch(SiconosException e)
//     {
//       if (ACE_MUET_LEVEL != ACE_MUET)
// 	std::cout << e.report() << endl;      
//       continue;
//     }
//   catch(...)
//     {
//       std::cout << "Exception caught." << endl;
//       continue;
//     }
   
//  };
//  if (find){
//    fillSolution();
//    return true;
//  }else{
//    ACE_MESSAGE("mlcp::solveGuessAndIt failed\n",9);
//    return false;
//  }
// }
// void   mlcp::fillSolution(){
//   unsigned int lin;
//   for (lin=0;lin<mDlcp;lin++){
//     if (mW1Z1[lin]==0){
//       //double aux =mQ->getValue(lin,0);
//       mW1->setValue(lin,0);
//       mZ1->setValue(lin,mQ->getValue(lin));
//     }else{
//       //double aux =mQ->getValue(lin,0);
//       mZ1->setValue(lin,0);
//       mW1->setValue(lin,mQ->getValue(lin));
//     }    
//   }

//   //  setBlock(mQ,mZ2,dim,start);
//   for (lin=0;lin<mDlin;lin++)
//     mZ2->setValue(lin,mQ->getValue(mDlcp+lin));
//   if((ACE_MUET_LEVEL != ACE_MUET )&& !mTringGuess)
//     printOutPut();
//  }
// void mlcp::addGuess(aceVector *Z){
//   unsigned long cur2pow=1;
//   unsigned long res=0;
  
//   for (unsigned int lin=0;lin<mDlcp;lin++){
//     if (Z->getValue(lin)<ACE_INF)
//       res+=cur2pow;
//     cur2pow=2*cur2pow;
//   }
//   cout<<"addguess :"<<res<<endl;
//   Z->display();
//   setCurrentConfig(res);
//   addGuess(res);
// }

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
  if (mQ)
    delete mQ;
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
  if (mA)
    free(mA);
  if(mB)
    free(mB);
  if(mC)
    free(mC);
  if (mD)
    free(mD);
  if(ma)
    free(ma);
  if(mb)
    free(mb);
  if(mu)
    free(mu);
//   if(mv)
//     free(mv);
  if(mw)
    free(mw);


}

// void mlcp::printGuess(ostream& os){
//   if (!mUseGuess || ACE_MUET_LEVEL == ACE_MUET)
//     ;//return;
//   ACE_MESSAGE("mlcp::printGuess\n");
//   initGuess();
//   int i=1;
//   while (tryGuess()){
//     os<<i<<" guess : "<<mCase<<endl;
//     i++;
//   }
// }

void mlcp::printInPutABCDab(ostream& os)
{
//   if (ACE_MUET_LEVEL == ACE_MUET)
//     return;
  os<<"mlcp print input"<<endl;
  os<<"dim lcp:\t"<<mDlcp<<"\tlin\t"<<mDlin;
  os<<"A\n";
  os<<(*mM22);
  os<<"B\n";
  os<<(*mM11);
  os<<"C\n";
  os<<(*mM21);
  os<<"D\n";
  os<<(*mM12);
  os<<"a\n";
  os<<(*mQ2);
  os<<"b\n";
  os<<(*mQ1);
  
}
void mlcp::printInPut(ostream& os)
{
  if (ACE_MUET_LEVEL == ACE_MUET)
    return;
  os<<"mlcp print input"<<endl;
  os<<"dim lcp:\t"<<mDlcp<<"\tlin\t"<<mDlin;
  os<<"M:\n";
  os<<(*mM);
  os<<"Q:\n"<<endl;
  os<<(*mQ);
}
void mlcp::printOutPut(ostream& os){
//   if (ACE_MUET_LEVEL == ACE_MUET )
//     return;
  os<<"mlcp print output"<<endl;
  os<<"Z1:\n";
  os<<(*mZ1);
  os<<"Z2:\n";
  os<<(*mZ2);
  os<<"W1:\n";
  os<<(*mW1);
}
