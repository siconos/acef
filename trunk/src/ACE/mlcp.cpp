/************************************************************************
  			mlcp.cpp 
**************************************************************************/
#include "mlcp.h"
#include "SimpleVector.h"
#include "mlcp_simplex.h"


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
   mSolverType=solverType;
   mTryOnlyGuess = false;
   mTringGuess = false;
   mUseGuess = true;
   if (mSolverType !=ACE_SOLVER_ENUM){
     mTryOnlyGuess = true;
     mA=(double *)calloc(mDlin*mDlin,sizeof(double));
     mC=(double *)calloc(mDlin*mDlcp,sizeof(double));
     mD=(double *)calloc(mDlcp*mDlin,sizeof(double));
     mB=(double *)calloc(mDlcp*mDlcp,sizeof(double));
     ma=(double *)calloc(mDlin,sizeof(double));
     mb=(double *)calloc(mDlcp,sizeof(double));
     mu=(double *)calloc(mDlin,sizeof(double));
     mv=(double *)calloc(mDlcp,sizeof(double));
     mw=(double *)calloc(mDlcp,sizeof(double));
   }
   
   mW1Z1=0;
   mM=0;
   mCurEnum=8388600;
   mCmp=0;
   mMaxEnum = (unsigned long) powl(2,(long)mDlcp);
   mPourCent=0;

   mM = new aceMatrix(Dlcp+Dlin,Dlcp+Dlin);
   mQ= new aceMatrix(Dlcp+Dlin,1);

   if (mDlcp){
     mW1Z1 = (int*)calloc(mDlcp,sizeof(int));
     //mW1Z1prev = (int*)calloc(mDlcp,sizeof(int));
     mW1 = new aceMatrix(mDlcp,1);
     mZ1 = new aceMatrix(mDlcp,1);
     mQ1=new aceMatrix(mDlcp,1);
     mM11 = new aceMatrix(mDlcp,mDlcp);
   }
   
   if (mDlin){
     mM22 = new aceMatrix(mDlin,mDlin);
     mZ2 = new aceMatrix(mDlin,1);
     mQ2 = new aceMatrix(mDlin,1);     
   }
   
   if (mDlin && mDlcp){
     mM21 = new aceMatrix(mDlin,mDlcp);
     mM12 = new aceMatrix(mDlcp,mDlin);
   }
   mGuess.clear();
   mGuess.push_back(0);
}
void mlcp::addGuess(unsigned long l){
  mGuess.push_back(l);
  cout <<"add guess "<<l<<endl;
}
void mlcp::affectW1Z1(unsigned long ll){
  unsigned long aux = ll;
  mCase = ll;
  for (unsigned int i =0; i <mDlcp; i++){
    mW1Z1[i]=aux & 1;
    aux = aux >> 1;
  }

}
void mlcp::setCurrentConfig(unsigned long l){
  (*(mGuess.begin()))=l;
}
void mlcp::initGuess(){
  mItGuess = mGuess.begin();
  ACE_MESSAGE("mlcp::initEnum : start guess\n");
}

bool mlcp::tryGuess(){
  if (mItGuess == mGuess.end()){
    ACE_MESSAGE("mlcp::tryGuess : stop guess\n");
    mTringGuess = false;
    ACE_times[ACE_TIMER_SOLVE_GUESS].stop();
    return false;
  }
  if (ACE_MUET_LEVEL != ACE_MUET)
    printf("try guess %d\n",(int)*mItGuess);
  affectW1Z1(*mItGuess);
  mItGuess++;
  return true;
}
void mlcp::initEnum(){
   mCurEnum=8000000;

  mPourCent=0;
  mCmp=0;
  if (mUseGuess){
    ACE_times[ACE_TIMER_SOLVE_GUESS].start();
    mTringGuess = true;
    initGuess();
  }else{
    mTringGuess = false;
  }
}

bool mlcp::nextEnum(){
  if (mTringGuess && tryGuess())
    return true;
  if (mTryOnlyGuess)
    return false;
  ACE_times[ACE_TIMER_SOLVE_ENUM].start();
  if (mCmp == mMaxEnum)
    return false;
  if (mCmp > mPourCent*mMaxEnum){
    mPourCent+=0.001;
    printf(" %f %d \n",mPourCent,(int) mCurEnum);
  }
  if (mCurEnum >= mMaxEnum){
    mCurEnum=0;
  }
  affectW1Z1(mCurEnum);
  mCurEnum++;
  mCmp++;
  
  return true;
}
bool mlcp::solveWithPath(){
  int n = mDlin;
  int m = mDlcp;
  bool res=false;
  
  mM22->MatrixToFortran(mA);  
  mM11->MatrixToFortran(mB);  
  mM21->MatrixToFortran(mC);  
  mM12->MatrixToFortran(mD);  
  mQ2->MatrixToFortran(ma);  
  mQ1->MatrixToFortran(mb);  
  ACE_times[ACE_TIMER_SOLVE_PATH].start();
  //  res= mlcp_path(&n , &m, mA , mB , mC , mD , ma, mb, mu, mv, mw , 0 , 0 , 0  );
  ACE_times[ACE_TIMER_SOLVE_PATH].stop();
  
  mZ1->FortranToMatrix(mv);
  mZ2->FortranToMatrix(mu);
  mW1->FortranToMatrix(mw);
  return res;

}
bool mlcp::solveWithSimplex(){
  int n = mDlin;
  int m = mDlcp;
  bool res;

  mM22->MatrixToFortran(mA);  
  mM11->MatrixToFortran(mB);  
  mM21->MatrixToFortran(mC);  
  mM12->MatrixToFortran(mD);  
  mQ2->MatrixToFortran(ma);  
  mQ1->MatrixToFortran(mb);
  
  ACE_times[ACE_TIMER_SOLVE_SIMPLEX].start();
  res= mlcp_simplex(&n , &m, mA , mB , mC , mD , ma, mb, mu, mv, mw , 0 , 0 , 0  );
  mCase = getConfigLCP();
  ACE_times[ACE_TIMER_SOLVE_SIMPLEX].stop();

  mZ1->FortranToMatrix(mv);
  mZ2->FortranToMatrix(mu);
  mW1->FortranToMatrix(mw);
  
  return res;
}
bool mlcp::solve(){
  bool res =false;
  //1) start with guess point
  ACE_MESSAGE("mlcp::solve : 1 try with guess point.\n");
  res = solveGuessAndIt();
  if (res){
    ACE_STOP_SOLVER_TIME();
    return true;
  }
  //2) start other algo
  ACE_MESSAGE("mlcp::solve : 2 start other algo\n");
  switch (mSolverType) {
  case ACE_SOLVER_SIMPLEX:
    res= solveWithSimplex();
    addGuess(mCase);
    break;
  case ACE_SOLVER_PATH:
    res= solveWithPath();
    break;
  case ACE_SOLVER_ENUM:
    break;
  default:
    ACE_ERROR("mlcp::solve, bad mSolverType value");
  }
  if (!res){
    ACE_MESSAGE("mlcp::solve, failed.\n");
  }else{
    ACE_MESSAGE("mlcp::solve : 3 Add a new guess point.\n");
    addGuess(mZ1);

  }
  ACE_STOP_SOLVER_TIME();
  return res;
}
bool mlcp::solveGuessAndIt(){
  bool onlyOne = true;

  bool find = false;
  unsigned int col =0;
  unsigned int lin =0;
  if (mDlin==0){//LCP case
    ACE_ERROR("mlcp::solveGuessAndIt : LCP not implemented");
  }else if(mDlcp==0){//Linear system case
    ACE_MESSAGE("mlcp::solveGuessAndIt : Linear system\n");
    try{
      (*mQ2)=-1*(*mQ2);
      //cout <<(*mM22);
      //cout <<(*mQ2);
      mM22->PLUForwardBackwardInPlace(*mQ2);
      for (lin=0;lin<mDlin;lin++)
	mZ2->setValue(lin,0,mQ2->getValue(lin,0));
      return true;
    }
    catch(SiconosException e)
    {
      std::cout << e.report() << endl;
      ACE_ERROR("linearSystem::solveGuessAndIt linear system no solution");
    }
    catch(...)
    {
      std::cout << "Exception caught." << endl;
      ACE_ERROR("linearSystem::solveGuessAndIt linear system no solution\n");
    }
  }
  
  if(ACE_MUET_LEVEL !=10 && onlyOne){
    onlyOne=false;
    mM->setBlock(0,0,*mM11);
    mM->setBlock(0,mDlcp,*mM12);
    mM->setBlock(mDlcp,0,*mM21);
    mM->setBlock(mDlcp,mDlcp,*mM22);
    mQ->setBlock(0,0,*mQ1);
    mQ->setBlock(mDlcp,0,*mQ2);
    printInPutABCDab();
  }
  initEnum();
  while(nextEnum()){
    //if(!mTringGuess){
    mM->setBlock(0,mDlcp,*mM12);
    mM->setBlock(mDlcp,mDlcp,*mM22);
    mQ->setBlock(0,0,*mQ1);
    mQ->setBlock(mDlcp,0,*mQ2);
    
    if(ACE_MUET_LEVEL !=10)
      printInPut();
   //build mM11_ and mM21_
   for (col =0; col<mDlcp; col++){
     if (mW1Z1[col]==0){
       //M11
       for (lin =0; lin < mDlcp; lin++)
	 mM->setValue(lin,col,mM11->getValue(lin,col));
       
       //M21
       for (lin =0; lin < mDlin; lin++)
	 mM->setValue(mDlcp+lin,col,mM21->getValue(lin,col));
       
     }else{
       //M11
       for (lin =0; lin < mDlcp; lin++)
	 mM->setValue(lin,col,0);
       mM->setValue(col,col,-1);
       
       //M21
       for (lin =0; lin < mDlin; lin++)
	 mM->setValue(mDlcp+lin,col,0);
     }
   }//build mM is built
   
  try{//solve the system
    if (ACE_MUET_LEVEL != ACE_MUET){
      ACE_MESSAGE("mlcp:: try: ");
      for (lin = 0 ; lin<mDlcp;lin++)
	cout<<mW1Z1[lin]<<" ";
      cout<<endl;
      printInPut();
    }
    //
   *mQ=-1*(*mQ);
   ACE_times[ACE_TIMER_DIRECT].start();
   mM->PLUForwardBackwardInPlace(*mQ);
   ACE_times[ACE_TIMER_DIRECT].stop();
   bool check=true;
   for (lin = 0 ; lin <mDlcp; lin++){
     double aux = mQ->getValue(lin,0);
     if (mQ->getValue(lin,0) < - ACE_NULL){
       ACE_MESSAGE("mlcp::solveGuessAndIt not in the cone\n");
       check=false;
       break;//pas dans le cone!
     }
   }
   if(!check)
     continue;
   ACE_MESSAGE("mlcp::solveGuessAndIt succes\n");
   if (ACE_MUET_LEVEL != ACE_MUET)
     printf("val de step %d \n",(int)mCase);
   if (mUseGuess && ! mTringGuess)
     addGuess(mCase);
   setCurrentConfig(mCase);
   find = true;
   break;// gagne.
  }
  catch(SiconosException e)
    {
      if (ACE_MUET_LEVEL != ACE_MUET)
	std::cout << e.report() << endl;      
      continue;
    }
  catch(...)
    {
      std::cout << "Exception caught." << endl;
      continue;
    }
   
 };
 if (find){
   for (lin=0;lin<mDlcp;lin++){
     if (mW1Z1[lin]==0){
       double aux =mQ->getValue(lin,0);
       mW1->setValue(lin,0,0);
       mZ1->setValue(lin,0,aux);
     }else{
       double aux =mQ->getValue(lin,0);
       mZ1->setValue(lin,0,0);
       mW1->setValue(lin,0,aux);
     }

   }
   for (lin=0;lin<mDlin;lin++)
     mZ2->setValue(lin,0,mQ->getValue(mDlcp+lin,0));
   if(ACE_MUET_LEVEL != ACE_MUET || true)
     printOutPut();
   return true;
 }else{
   ACE_MESSAGE("mlcp::solveGuessAndIt failed\n",9);
   return false;
 }
}
void mlcp::addGuess(aceMatrix *Z){
  unsigned long cur2pow=1;
  unsigned long res=0;
  
  for (unsigned int lin=0;lin<mDlcp;lin++){
    if (Z->getValue(lin,0)<ACE_INF)
      res+=cur2pow;
    cur2pow=2*cur2pow;
  }
  cout<<"addguess :"<<res<<endl;
  Z->display();
  setCurrentConfig(res);
  addGuess(res);
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
  if(mv)
    free(mv);
  if(mw)
    free(mw);


}

void mlcp::printGuess(ostream& os){
  if (!mUseGuess || ACE_MUET_LEVEL == ACE_MUET)
    return;
  ACE_MESSAGE("mlcp::printGuess\n");
  initGuess();
  while (tryGuess())
    os<<"guess : "<<mCase<<endl;
}

void mlcp::printInPutABCDab(ostream& os)
{
  if (ACE_MUET_LEVEL == ACE_MUET)
    return;
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
  if (ACE_MUET_LEVEL == ACE_MUET )
    return;
  os<<"mlcp print output"<<endl;
  os<<"Z1:\n";
  os<<(*mZ1);
  os<<"Z2:\n";
  os<<(*mZ2);
  os<<"W1:\n";
  os<<(*mW1);
}
