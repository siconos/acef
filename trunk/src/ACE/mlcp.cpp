/************************************************************************
  			mlcp.cpp 
**************************************************************************/
#include "mlcp.h"
#include "SimpleVector.h"


mlcp::mlcp(unsigned int Dlcp,unsigned int Dlin){

   mDlcp = Dlcp;
   mDlin = Dlin;
   ACE_CHECK_IERROR(mDlcp+mDlin>0,"mlcp::mlcp dim null");
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
   //mFirst=true;
   
   mW1Z1=0;
   mM=0;
   mCurEnum=0;
   mCmp=0;
   mMaxEnum = (unsigned long) powl(2,(long)mDlcp);

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
}

void mlcp::initEnum(){
  unsigned long aux;
  mCmp=0;
  aux = mCurEnum;
  for (unsigned int i =0; i <mDlcp; i++){
    mW1Z1[i]=aux & 1;
    aux = aux >> 1;
  }
}

bool mlcp::nextEnum(){
  if (mCmp == mMaxEnum)
    return false;
  unsigned long aux;
  mCmp++;
  mCurEnum++;
  aux = mCurEnum;
  for (unsigned int i =0; i <mDlcp; i++){
    mW1Z1[i]=aux & 1;
    aux = aux >> 1;
  }

  
  return true;
  /*for (unsigned int i =0; i <mDlcp; i++){
    if (mW1Z1[i]==0){
      mW1Z1[i]=1;
      return true;
    }
    mW1Z1[i]=0;
  }    
  return false;*/
}
bool mlcp::solve(){
  bool muet = true;
  bool find = false;
  unsigned int col =0;
  unsigned int lin =0;
  if (mDlin==0){//LCP case
    ACE_ERROR("mlcp::solve : LCP not implemented");
  }else if(mDlcp==0){//Linear system case
    ACE_MESSAGE("mlcp::solve : Linear system\n");
    try{
      (*mQ2)=-1*(*mQ2);
      cout <<(*mM22);
      cout <<(*mQ2);
      mM22->PLUForwardBackwardInPlace(*mQ2);
      for (lin=0;lin<mDlin;lin++)
	mZ2->setValue(lin,0,mQ2->getValue(lin,0));
      return true;
    }
    catch(SiconosException e)
    {
      std::cout << e.report() << endl;
      ACE_ERROR("linearSystem::solve linear system no solution");
    }
    catch(...)
    {
      std::cout << "Exception caught." << endl;
      ACE_ERROR("linearSystem::solve linear system no solution\n");
    }
  }
  
  initEnum();
  if(!muet){
    mM->setBlock(0,0,mM11);
    mM->setBlock(0,mDlcp,mM12);
    mM->setBlock(mDlcp,0,mM21);
    mM->setBlock(mDlcp,mDlcp,mM22);
    mQ->setBlock(0,0,mQ1);
    mQ->setBlock(mDlcp,0,mQ2);
    printInPut();
  }
  do{
    mM->setBlock(0,mDlcp,mM12);
    mM->setBlock(mDlcp,mDlcp,mM22);
    mQ->setBlock(0,0,mQ1);
    mQ->setBlock(mDlcp,0,mQ2);
    
    if(!muet)
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
    if (!muet){
      ACE_MESSAGE("mlcp:: try: ");
      for (lin = 0 ; lin<mDlcp;lin++)
	cout<<mW1Z1[lin]<<" ";
      cout<<endl;
      printInPut();
    }
    //mM->PLUFactorizationInPlace();
    //if(!mM->isFactorized())
    //  continue;

   //cout <<"mM:\n"<<(*mM)<<"\n";
   *mQ=-1*(*mQ);
   mM->PLUForwardBackwardInPlace(*mQ);
   bool check=true;
   for (lin = 0 ; lin <mDlcp; lin++){
     double aux = mQ->getValue(lin,0);
     if (mQ->getValue(lin,0) < - ACE_NULL){
       if(!muet)
	 ACE_MESSAGE("mlcp::solve not in the cone\n");
       check=false;
       break;//pas dans le cone!
     }
   }
   if(!check)
     continue;
   if(!muet)
     ACE_MESSAGE("mlcp::solve succes\n");
   find = true;
   break;// gagne.
  }
  catch(SiconosException e)
    {
      if (!muet)
	std::cout << e.report() << endl;      
      continue;
    }
  catch(...)
    {
      std::cout << "Exception caught." << endl;
      continue;
    }
   
 }while(nextEnum());
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
   if(!muet)
     printOutPut();
   return true;
 }else{
   ACE_MESSAGE("mlcp::solve failed\n");
   return false;
 }
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
}

void mlcp::printInPut(ostream& os)
{
  os<<"mlcp print input"<<endl;
  os<<"dim lcp:\t"<<mDlcp<<"\tlin\t"<<mDlin;
  os<<"M:\n";
  os<<(*mM);
  os<<"Q:\n"<<endl;
  os<<(*mQ);
}
void mlcp::printOutPut(ostream& os){
  os<<"mlcp print input"<<endl;
  os<<"Z1:\n";
  os<<(*mZ1);
  os<<"Z2:\n";
  os<<(*mZ2);
  os<<"W1:\n";
  os<<(*mW1);
}
