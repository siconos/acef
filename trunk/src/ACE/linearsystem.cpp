/************************************************************************
  			linearsystem.cpp 
**************************************************************************/


#include "linearsystem.h"
#include <fstream>
#include "SimpleMatrix.h"

using namespace std;
//|--------x'---------|----------x--------|-------Zs------|------Zns------|

linearSystem::linearSystem(){
  
  mReAlloc = false;
  mDynIndex=0;
  mNbNodes=0;
  mNbUnknowns=0;
  mNbEquations=0;
  mNbDynEquations=0;
  mNbNonDynEquations=0;
  mDimLambda=0;
  mA=0;
  mB=0;
  mC=0;
  mD=0;
  ms=0;

  mA1x=0;
  mA1zs=0;
  mA1zns=0;
  mA1s=0;

  mB1x=0;
  mB1zs=0;
  mB1zns=0;
  mB1s=0;

  mC1x=0;
  mC1zs=0;
  mC1l=0;
  mC1s=0;

  mD1x=0;
  mD1zs=0;
  mD1zns=0;
  mD1l=0;
  mD1s=0;

  mR=0;
  mA2x=0;
  mA2zs=0;
  mA2s=0;
  mA2sti=0;
  
  mB2x=0;
  mB2zs=0;
  mB2l=0;
  mB2s=0;

  mD2x=0;
  mD2zs=0;
  mD2l=0;
  mD2s=0;

  mDimx=0;
  mDimzs=0;
  mDimzns=0;

  mxti=0;
  mzsti=0;
  mznsti=0;
  mxfree=0;

  mW=0;
  mD3l=0;
  mD3zs=0;
  mB3l=0;
  mB3zs=0;
  mPfree=0;
  mPAux=0;
  mPxAux=0;
  mQfree=0;
  mMLCP=0;

  mTheta = 0.5;
  mThetap = 0.5;
  mH = 1;
  mLogFrequency=0;
  mLogPrint=0;
  mPourMille=0;

  //
  mD2xW=0;
  mB2xW=0;
  mHThetaWA2zs=0;
  mHWR=0;
}
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////MEMORY


void linearSystem::allocMemory(){
  ACE_CHECK_IERROR(!mReAlloc,"linearSystem::allocMemory again");
  mReAlloc = true;
  mNbUnknowns = 2*mx.size()+mZs.size()+mZns.size();
  mRS = mNbUnknowns;
  int i;
  int n=mKCL.size();
  for (i=0;i<n;i++)
    mKCL[i]->allocMemory(mNbUnknowns+1);
  n=mVD.size();
  for (i=0;i<n;i++)
    mVD[i]->allocMemory(mNbUnknowns+1);
  n=mTEN.size();
  for (i=0;i<n;i++)
    mTEN[i]->allocMemory(mNbUnknowns+1);
  n=mIND.size();
  for (i=0;i<n;i++)
    mIND[i]->allocMemory(mNbUnknowns+1);
  n=mCAP.size();
  for (i=0;i<n;i++)
    mCAP[i]->allocMemory(mNbUnknowns+1);


  //Matrix allocation
  allocA1Matrix();
  mNbNonDynEquations = mNbEquations - mNbDynEquations;
  ACE_CHECK_IERROR(mNbNonDynEquations == mDimzs-1,"linearSystem::preparForStamp, mNbNonDynEquations and mDimzs-1 not coherent");
  ACE_CHECK_IERROR(mNbNonDynEquations >=0,"linearSystem::preparForStamp, mNbNonDynEquations <0.");
  allocB1Matrix();
  allocC1Matrix();
  allocD1Matrix();
  if (mDimLambda && mDimx){
    mD2xW=new aceMatrix(mDimLambda,mDimx,ACE_MAT_TYPE);
    mHWR=new aceMatrix(mDimx,mDimLambda,ACE_MAT_TYPE);
  }
  if (mNbNonDynEquations && mDimx)
    mB2xW=new aceMatrix(mNbNonDynEquations,mDimx,ACE_MAT_TYPE);
  if (mDimx)
    mHThetaWA2zs=new aceMatrix(mDimx,mDimzs-1,ACE_MAT_TYPE);

  
}
//x'=A1x*x + A1zs*Zs + A1zns*Zns + A1s
void linearSystem::allocA1Matrix(){
  ACE_CHECK_IERROR(!mA,"linearSystem::buildABCDs : mA not NULL");
  ACE_CHECK_IERROR(!mB,"linearSystem::buildABCDs : mB not NULL");

  
  if (mDimx != 0 ){
    mA = new aceMatrix(mDimx,mDimx,ACE_MAT_TYPE);
    mB= new aceMatrix(mDimx,mDimx,ACE_MAT_TYPE);
    mPxAux=new aceMatrix(mDimx,mDimx);

    if (mDimzs-1>0){
      mC=new aceMatrix(mDimx,mDimzs-1,ACE_MAT_TYPE);
      mMatBuf1=new aceMatrix(mDimx,mDimzs-1,ACE_MAT_TYPE);
    }
    if (mDimzns)
      mD=new aceMatrix(mDimx,mDimzns,ACE_MAT_TYPE);
    ms=new aceVector(mDimx,ACE_MAT_TYPE);
    
    mA1x = mB;
    mA1zs= mC;
    mA1zns= mD;
    mA1s=ms;
  }
}
void linearSystem::freeA1Matrix(){

  if (mA)
    delete mA;
  if (mB)
    delete mB;
  if (mC)
    delete mC;
  if (mMatBuf1)
    delete mMatBuf1;
  if (mD)
    delete mD;
  if (ms)
    delete ms;
}

//0 = B1x*x + B1zs*Zs + B1zns*Zns + B1s
void linearSystem::allocB1Matrix(){
  if (!mNbNonDynEquations)
    return;
  ACE_CHECK_IERROR(mB1zs ==0,"linearSystem::allocB1Matrix again");
  if (mDimx)
    mB1x=new aceMatrix(mNbNonDynEquations,mDimx,ACE_MAT_TYPE);
  if (mDimzs-1>0)
    mB1zs=new aceMatrix(mNbNonDynEquations,mDimzs-1,ACE_MAT_TYPE);
  if (mDimzns)
    mB1zns=new aceMatrix(mNbNonDynEquations,mDimzns,ACE_MAT_TYPE);
  mB1s=new aceVector(mNbNonDynEquations,ACE_MAT_TYPE);
}
void linearSystem::freeB1Matrix(){
  if (mB1x)
    delete mB1x;
  if (mB1zs)
    delete mB1zs;
  if (mB1zns)
    delete mB1zns;
  if (mB1s)
    delete mB1s;
  
}
void linearSystem::allocC1Matrix(){
  ACE_CHECK_IERROR(mC1zs == 0,"linearSystem::allocC1Matrix");
  if (!mDimzns)
    return;
  if (mDimx)
    mC1x=new aceMatrix(mDimzns,mDimx,ACE_MAT_TYPE);
  if (mDimzs-1>0)
    mC1zs=new aceMatrix(mDimzns,mDimzs-1,ACE_MAT_TYPE);
  if(mDimLambda)
    mC1l = new aceMatrix(mDimzns,mDimLambda,ACE_MAT_TYPE);
  mC1s = new aceVector(mDimzns,ACE_MAT_TYPE);

}
void linearSystem::freeC1Matrix(){
  if (mC1x)
    delete mC1x;
  if (mC1zs)
    delete mC1zs;
  if (mC1l)
    delete mC1l;
  if (mC1s)
    delete mC1s;
}
void linearSystem::allocD1Matrix(){
  ACE_CHECK_IERROR(mD1zs == 0,"linearSystem::allocD1Matrix");
  if (!mDimLambda)
    return;
  if (mDimx)
    mD1x=new aceMatrix(mDimLambda,mDimx,ACE_MAT_TYPE);
  if (mDimzs-1>0)
    mD1zs=new aceMatrix(mDimLambda,mDimzs-1,ACE_MAT_TYPE);
  if (mDimzns)
    mD1zns=new aceMatrix(mDimLambda,mDimzns,ACE_MAT_TYPE);
  mD1l = new aceMatrix(mDimLambda,mDimLambda,ACE_MAT_TYPE);
  mD1s = new aceVector(mDimLambda,ACE_MAT_TYPE);
  if (mDimx){
    mR = new aceMatrix(mDimx,mDimLambda,ACE_MAT_TYPE);
    mMatBuf2 = new aceMatrix(mDimx,mDimLambda,ACE_MAT_TYPE);
  }

}
void linearSystem::set2matrix(){

  mA2x=mA1x;
  mA2zs=mA1zs;
  
  mB2x=mB1x;
  mB2zs=mB1zs;
  mB2l=new aceMatrix(mNbNonDynEquations,mDimLambda,ACE_MAT_TYPE);

  mD2x=mD1x;
  mD2zs=mD1zs;
  mD2l=mD1l;
  mA2s=mA1s;
  if (mDimx)
    mA2sti = new aceVector(mDimx,ACE_MAT_TYPE);
  mB2s=mB1s;
  mD2s=mD1s;

  if (mDimzns){
    if (mDimx){
      ACEprod(*mA1zns,*mC1l,*mR);
      //      *mA2x = *mA1x + prod(*mA1zns,*mC1x);
      ACEprod(*mA1zns,*mC1x,*mA2x,false);
      if (mDimzs)
	//*mA2zs = *mA1zs + prod(*mA1zns,*mC1zs);
	ACEprod(*mA1zns,*mC1zs,*mA1zs,false);
    }
    if (mNbNonDynEquations){
      if (mDimx)
	//	*mB2x = *mB1x + prod(*mB1zns,*mC1x);
	ACEprod(*mB1zns,*mC1x,*mB1x,false);
      if (mDimzs)
	//	*mB2zs = *mB1zs + prod(*mB1zns,*mC1zs);
	ACEprod(*mB1zns,*mC1zs,*mB1zs,false);
      if (mDimx)
	//	*mD2x=*mD1x + prod(*mD1zns,*mC1x);
	ACEprod(*mD1zns,*mC1x,*mD1x,false);
      //      *mD2zs=*mD1zs+prod(*mD1zns,*mC1zs);
      ACEprod(*mD1zns,*mC1zs,*mD1zs,false);
      //      *mD2l=*mD1l + prod(*mD1zns,*mC1l);
      ACEprod(*mD1zns,*mC1l,*mD1l,false);

      //      *mB2l=prod(*mB1zns,*mC1l);
      ACEprod(*mB1zns,*mC1l,*mB2l,true);

    }
  }

}
void linearSystem::freeD1Matrix(){
  if (mD1x)
    delete mD1x;
  if (mD1zs)
    delete mD1zs;
  if (mD1zns)
    delete mD1zns;
  if (mD1l)
    delete mD1l;
  if (mD1s)
    delete mD1s;
  if (mR)
    delete mR;
  if (mMatBuf2)
    delete mMatBuf2;
  if (mB2l)
    delete mB2l;
  if (mA2sti)
    delete mA2sti;
}
void linearSystem::allocForInitialValue(){
  if(mDimx){
    mxti = new aceVector(mDimx,ACE_MAT_TYPE);
  }
  mzsti = new aceVector(mDimzs-1,ACE_MAT_TYPE);
}

void linearSystem::allocDiscretisation(){
  allocForInitialValue();
  if(mDimx){
    mxfree = new aceVector(mDimx,ACE_MAT_TYPE);
    mW = new aceMatrix(mDimx,mDimx,ACE_MAT_TYPE);
  }
  
  if(mDimzns)
    mznsti = new aceVector(mDimzns,ACE_MAT_TYPE);
  if (mB2zs)
    mB3zs = new aceMatrix(mNbNonDynEquations,mDimzs-1,ACE_MAT_TYPE);
  if (mNbNonDynEquations && mDimLambda)
    mB3l=new aceMatrix(mNbNonDynEquations,mDimLambda,ACE_MAT_TYPE);
  if (mNbNonDynEquations)
    mQfree = new aceVector(mNbNonDynEquations,ACE_MAT_TYPE);

  if (mDimLambda)
    mD3zs = new aceMatrix(mDimLambda,mDimzs-1,ACE_MAT_TYPE);
  if (mDimLambda)
    mD3l = new aceMatrix(mDimLambda,mDimLambda,ACE_MAT_TYPE);
  if (mDimLambda){
    mPfree = new aceVector(mDimLambda,ACE_MAT_TYPE);
    mPAux = new aceVector(mDimLambda,ACE_MAT_TYPE);

  }
}
void linearSystem::freeForInitialValue(){
  if (mxti)
    delete mxti;
  if (mzsti)
    delete mzsti;
}
void linearSystem::freeDiscretisation(){
  freeForInitialValue();
  if (mxfree)
    delete mxfree;
  if (mznsti)
    delete mznsti;
  if(mW)
    delete mW;
  if (mD3l)
    delete mD3l;
  if (mD3zs)
    delete mD3zs;
  if (mB3l)
    delete mB3l;
  if (mB3zs)
    delete mB3zs;
  if (mPfree)
    delete mPfree;
  if (mPAux)
    delete mPAux;
  if (mPxAux)
    delete mPxAux;
  if (mQfree)
    delete mQfree;
}

linearSystem::~linearSystem(){
  int i =0;
  //free unknows
  int n = mx.size();
  for (i=0; i<n; i++)
    delete mx[i];
  n = mZs.size();
  for (i=0; i<n; i++)
    delete mZs[i];
  n = mZns.size();
  for (i=0; i<n; i++)
    delete mZns[i];

  //free equations
  n=mKCL.size();
  for (i=0;i<n;i++)
    delete mKCL[i];
  n=mVD.size();
  for (i=0;i<n;i++)
    delete mVD[i];
  n=mTEN.size();
  for (i=0;i<n;i++)
    delete mTEN[i];
  n=mIND.size();
  for (i=0;i<n;i++)
    delete mIND[i];
  n=mCAP.size();
  for (i=0;i<n;i++)
    delete mCAP[i];
  freeA1Matrix();
  freeB1Matrix();
  freeC1Matrix();
  freeD1Matrix();
  freeDiscretisation();
  if (mD2xW)
    delete mD2xW;
  if (mB2xW)
    delete mB2xW;
  if (mHThetaWA2zs)
    delete mHThetaWA2zs;
  if (mHWR)
    delete mHWR;

}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////DISCRETISATION
void linearSystem::readInitialValue(){
  int i;
  double aux;
  int nbGuess;
  int useIc;
  unsigned long guess;

  cout << "get init value from Netlist\n";
  try{

    nbGuess=0;
    double stop,start;
    getTransValues(&mH,&stop,&start);
    mStepNumber = (long)(stop/mH);
    cout << "mStepNumber\t"<<mStepNumber<<"\tmH\t"<<mH<<endl;
    for (i=0;i<mDimzs-1;i++){
      mzsti->setValueIfNotNull(i,0);
    }
    initICvalue();
    while(getICvalue(&i,&useIc,&aux)){
      if (i>0){
	mzsti->setValueIfNotNull(i-1,aux);
	cout<<"set value from netlist :v_"<<i<<"="<<aux<<endl;
      }
    }
      
    for (i=0;i<mDimx;i++){
      mxti->setValueIfNotNull(i,0);
    }
    
  }
  catch(...)
    {
      std::cout << "Exception caught." << endl;
      ACE_INTERNAL_ERROR("linearSystem::readInitialValue");
    }

}
void linearSystem::initSimu(){
  allocDiscretisation();
  mMLCP = new mlcp(mDimLambda,mNbNonDynEquations,ACE_SOLVER_TYPE);
  readInitialValue();
  mStepCmp=0;
  mLogFrequency =  mStepNumber/1000;
  mLogPrint = mStepNumber/100000;
  if (mLogPrint==0) mLogPrint=1;
  if (mLogFrequency==0) mLogFrequency = 1;
  mPourMille=0;

  try{
    //ACE_CHECK_IERROR(mDimx >0 && mDimLambda >0,"linearSystem::initSimu case no x or no lambda not yet implemented");
    if (mDimx){
      mW->zero();
      for(int i = 0; i < mDimx;i++)	mW->setValue(i,i,1.0);
      *mW -=  mH*mTheta*(*mA2x);
      if (ACE_MAT_TYPE == SPARSE){
	mPxAux->set(*mW);
	mPxAux->PLUInverseInPlace();
	mW->set(*mPxAux);
      }else{
	mW->PLUInverseInPlace();
      }
      
      //*mB3zs = *mB2zs + mH*mTheta*prod(*mB2x,prod(*mW,*mA2zs)) ;
      ACEprod(*mW,*mA2zs,*mMatBuf1,true);
      ACEprod(*mB2x,*mMatBuf1,*mB3zs,true);
      scal(mH*mTheta,*mB3zs,*mB3zs);
      *mB3zs+=*mB2zs;
      
      if (mDimLambda){
	//*mD3zs = *mD2zs+mH*mTheta*prod(*mD2x,prod(*mW,*mA2zs));
	ACEprod(*mD2x,*mMatBuf1,*mD3zs,true);
	scal(mH*mTheta,*mD3zs,*mD3zs);
	*mD3zs+=*mD2zs;

      
	//*mB3l = mH*prod(*mB2x,prod(*mW,*mR))+*mB2l;
	ACEprod(*mW,*mR,*mMatBuf2,true);
	ACEprod(*mB2x,*mMatBuf2,*mB3l,true);
	scal(mH,*mB3l,*mB3l);
	*mB3l+=*mB2l;
	     
	//*mD3l = mH*prod(*mD2x,prod(*mW,*mR))+*mD2l;
	ACEprod(*mD2x,*mMatBuf2,*mD3l,true);
	scal(mH,*mD3l,*mD3l);
	*mD3l+=*mD2l;
      }
    }else{
      *mD3zs = *mD2zs;
      *mB3zs = *mB2zs;
      if (mDimLambda){
	*mB3l = *mB2l;
	*mD3l = *mD2l;
      }
    }
    *(mMLCP->mM11) = *mD3l;
    *(mMLCP->mM12) = *mD3zs;
    *(mMLCP->mM21) = *mB3l;
    *(mMLCP->mM22) = *mB3zs;
  }
  catch(SiconosException e)
    {
      std::cout << e.report() << endl;
      ACE_INTERNAL_ERROR("linearSystem::initSimu");
    }
  catch(...)
    {
      std::cout << "Exception caught." << endl;
      ACE_INTERNAL_ERROR("linearSystem::initSimu");
    }
  if (mDimLambda && mDimx){
    //  *mD2xW=prod(*mD2x,*mW);
    ACEprod(*mD2x,*mW,*mD2xW,true);
  }
  if (mDimx){
    //*mB2xW=prod(*mB2x,*mW);
    ACEprod(*mB2x,*mW,*mB2xW,true);
    //*mHThetaWA2zs=mH*mTheta*prod(*mW,*mA2zs);
    ACEprod(*mW,*mA2zs,*mHThetaWA2zs,true);
    scal(mH*mTheta,*mHThetaWA2zs,*mHThetaWA2zs);
  }
  //*mHWR=mH*prod(*mW,*mR);
  if (mDimLambda && mDimx){
    ACEprod(*mW,*mR,*mHWR,true);
    scal(mH,*mHWR,*mHWR);
  }

  mMLCP->initSolver();

  
}
void linearSystem::preparStep(){
  if (mDimx && mDimLambda){//both
    //(*mxfree) = (*mxti) + mH*((1-mTheta)*(prod(*mA2x,*mxti)+prod(*mA2zs,*mzsti) + (*mA2sti))+mThetap*(*mA2s));
    //mxfree->display();
    ACEprod(*mA2zs,*mzsti,*mxfree,true);
    ACEprod(*mA2x,*mxti,*mxfree,false);
    *mxfree+=*mA2sti;
    scal((1-mTheta),*mxfree,*mxfree);
    *mxfree+=mThetap*(*mA2s);
    scal(mH,*mxfree,*mxfree);
    *mxfree+=*mxti;
    //mxfree->display();


    //(*mPfree) = prod(*mD2xW,*mxfree) + (*mD2s);
    ACEprod(*mD2xW,*mxfree,*mPfree,true);
    *mPfree+=*mD2s;
    
    //(*mQfree) = prod(*mB2xW,*mxfree) + (*mB2s);
    ACEprod(*mB2xW,*mxfree,*mQfree,true);
    *mQfree+=*mB2s;



    
  } else if(mDimx){//only x
    (*mxfree) = (*mxti) + mH*(1-mTheta)*prod(*mA2x,*mxti)+mH*(1-mTheta)*prod(*mA2zs,*mzsti)+mH*(*mA2s);
    (*mQfree) = prod(*mB2xW,*mxfree) + (*mB2s);
  }else if (mDimLambda){//only lambda
    (*mPfree) = (*mD2s);
    (*mQfree) = (*mB2s);
  }else{
    (*mQfree) = (*mB2s);
  }
//   cout<<"mQfree\n";
//   (mQfree)->display();
  if(ACE_MAT_TYPE == SPARSE){
    (mMLCP->mQ1)->set(*(mPfree));
    scal(-1,*(mMLCP->mQ1),*(mMLCP->mQ1));
    //  *(mMLCP->mQ1)= *mPfree;  
    //*(mMLCP->mQ2)= *mQfree;
    (mMLCP->mQ2)->set(*(mQfree));
    scal(-1,*(mMLCP->mQ2),*(mMLCP->mQ2));
  }else{
    if (mDimLambda)
      scal(-1,*mPfree,*(mMLCP->mQ1));
    scal(-1,*mQfree,*(mMLCP->mQ2));
    //  *(mMLCP->mQ1)= *mPfree;  
    //*(mMLCP->mQ2)= *mQfree;
  }
//   cout<<"mMLCP->mQ1\n";
//   (mMLCP->mQ1)->display();
//   cout<<"mMLCP->mQ2\n";
//   (mMLCP->mQ2)->display();

}
void linearSystem::computeZnstiFromX_Zs(){
  if (mDimLambda==0)
    return;
  if (mDimx){
    //(*mznsti)=prod(*mC1x,*mxti)+prod(*mC1zs,*mzsti)+prod(*mC1l,*(mMLCP->mZ1))+(*mC1s);
    ACE_times[ACE_TIMER_PROD_MAT].start();
    if (ACE_MAT_TYPE==SPARSE)
      ACEprod(*mC1l,*mPAux,*mznsti,true);
    else
      ACEprod(*mC1l,*(mMLCP->mZ1),*mznsti,true);

    ACEprod(*mC1zs,*mzsti,*mznsti,false);
    ACEprod(*mC1x,*mxti,*mznsti,false);
    ACE_times[ACE_TIMER_PROD_MAT].stop();    
    *mznsti+=*mC1s;
  }  else
    (*mznsti)=prod(*mC1zs,*mzsti)+prod(*mC1l,*(mMLCP->mZ1))+(*mC1s);
}
bool linearSystem::step(){
  ACE_times[ACE_TIMER_LS_STEP].start();
  mStepCmp++;
  if (mStepCmp >= mStepNumber){
    ACE_times[ACE_TIMER_LS_STEP].stop();
    return false;
  }
  if (mStepCmp%mLogFrequency==0){
    printf("-->%d\n",mPourMille);
    mPourMille ++;
  }
  bool res = mMLCP->solve();
  ACE_times[ACE_TIMER_COMPUTE_VAR].start();
  if (res){
    //*mzsti=*(mMLCP->mZ2);
    if (ACE_MAT_TYPE == SPARSE)
      mzsti->set(*(mMLCP->mZ2));
    else
      *mzsti=*(mMLCP->mZ2);
    if (mDimx){
      //*mxti=prod(*mW,*mxfree)+prod(*mHThetaWA2zs,*mzsti)+prod(*mHWR,*(mMLCP->mZ1));
      if (ACE_MAT_TYPE == SPARSE){
	mPAux->set(*(mMLCP->mZ1));
	//cout<<"mHWR "<<(*mHWR);
	//cout<<"mPAux "<<(*mPAux);
	ACEprod(*mHWR,*mPAux,*mxti,true);
      }else{
	if (mDimLambda)
	  ACEprod(*mHWR,*(mMLCP->mZ1),*mxti,true);
	else
	  mxti->zero();
      }
      //cout<<"mHThetaWA2zs "<<(*mHThetaWA2zs);
      //cout<<"mW "<<(*mW);
      ACE_times[ACE_TIMER_PROD_MAT].start();
      ACEprod(*mHThetaWA2zs,*mzsti,*mxti,false);
      ACEprod(*mW,*mxfree,*mxti,false);
      ACE_times[ACE_TIMER_PROD_MAT].stop();
  }
    computeZnstiFromX_Zs();
  }
  ACE_times[ACE_TIMER_COMPUTE_VAR].stop();
  ACE_times[ACE_TIMER_LS_STEP].stop();
  return res;
}
void linearSystem::stopSimu(){
  mMLCP->printGuess();
  mMLCP->stopSolver();
  if (mMLCP)
    delete mMLCP;
  mMLCP=0;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////UNKNOWN
unknown* linearSystem::addinx(int type, component* c)
{
  unknown* res = new unknown(type,c);
  res->mIndexInVector = mx.size();
  mx.push_back(res);
  return res;
}
unknown* linearSystem::addinZs(int type, component* c)
{
  unknown* res = new unknown(type,c);
  res->mIndexInVector = mZs.size();
  mZs.push_back(res);
  return res;
}
unknown* linearSystem::addinZns(int type, component* c)
{
  unknown* res = new unknown(type,c);
  res->mIndexInVector = mZns.size();
  mZns.push_back(res);
  return res;
}

int linearSystem::getIndexUnknown (int type,int node) {
  if (type == ACE_TYPE_V){
    return node + 2*mx.size();
  }else{
    ACE_INTERNAL_WARNING("linearSystem::getIndexUnknown not implemented");
    return -1;
  }
  
}
void  linearSystem::addVUnknowns(){
  ACE_CHECK_IERROR(mZs.size()==0,"linearSystem::initAddVUnknowns mZs not empty");
  ACE_CHECK_IERROR(mNbNodes,"linearSystem::initAddVUnknowns no nodes");
 for(int i=0;i<mNbNodes;i++){
   mZs.push_back(new unknown(ACE_TYPE_V,i));
 }
}
/**
 * 
 */
bool linearSystem::isUnknown (int type, component* c,unknown **uout=0) {
  if(type == ACE_TYPE_U){
      int ii;
      for(ii=0; ii <(int) mx.size(); ii++)
	{
	  unknown* u = mx[ii];
	  
	  if (u->mType == ACE_TYPE_U && (
	      (u->mComponent->mNodePos == c->mNodePos && u->mComponent->mNodeNeg == c->mNodeNeg ) ||
	      (u->mComponent->mNodePos == c->mNodeNeg && u->mComponent->mNodeNeg == c->mNodePos ))){
	    if (uout)
	      (*uout)=u;
	    return true;
	  }
	}
      return false;
  }else{
    ACE_INTERNAL_ERROR("linearSystem::isUnknown not implemented");
    return false;
  }
  
}
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////EQUATIONS
equationKCL* linearSystem::KCL(int i){
  if (i >= (int)mKCL.size())
    ACE_ERROR("linearSystem::KCL");
  else
    return (equationKCL*) mKCL[i];
  return 0;
}
void linearSystem::initKCL(){
  ACE_CHECK_IERROR(mKCL.size()==0,"linearSystem::initKCL with mKCL not empty!");
  for (int i=0;i<mNbNodes;i++){
    mKCL.push_back(new equationKCL(i));
  }
  mKCL[0]->mAvailable=false;
  mNbEquations+=mNbNodes-1;
}
equationCAP* linearSystem::addCapEquation(){
  equationCAP* res = new equationCAP();
  mNbEquations++;
  mNbDynEquations++;
  res->mLine = mDynIndex;
  mDynIndex++;
  mCAP.push_back(res);
  return res;
}
equationIND* linearSystem::addIndEquation(){
 equationIND* res = new equationIND();
 mNbEquations++;
 mNbDynEquations++;
 res->mLine = mDynIndex;
 mDynIndex++;
 mIND.push_back(res);
 return res;
}
equationVD* linearSystem::addVdEquation(char* name){
  equationVD* res = new equationVD(name);
  mNbEquations++;
  mVD.push_back(res);
  return res;

}
equationTEN* linearSystem::addTenEquation(){
   equationTEN* eq = new equationTEN();
   mNbEquations++;
   mTEN.push_back(eq);
   return eq;
}

void linearSystem::addKCLinDyn(int j){
  ACE_CHECK_IERROR(j <= mNbNodes,"linearSystem::addKCLinDyn");
  if (!mKCL[j]->mIsDyn)
    mNbDynEquations++;
  else
    ACE_INTERNAL_WARNING("linearSystem::addKCLinDyn kcl in dyn.");
  
  mKCL[j]->mLine = mDynIndex;
  mKCL[j]->mIsDyn = true;
  mKCL[j]->mAvailable = false;
  mKCL[j]->mLine=mDynIndex;
  mDynIndex++;

}

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////STAMP

void linearSystem::preparForStamp(){
  int i =0;
  mDimx = mx.size();
  for (i=0; i<mDimx; i++){
    mx[i]->mDynIndex = i;
    mx[i]->mIndex = mDimx+i;    
  }
  mDimzs=mZs.size();
  for (i=0; i<mDimzs; i++){
    mZs[i]->mIndex = 2*mDimx+i;
  }
  mDimzns=mZns.size();
  for (i=0; i<mDimzns; i++){
    mZns[i]->mIndex = 2*mDimx+mDimzs+i;
  }

  allocMemory();

}





///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////DYNAMICS EQUATIONS  : x'=A1x*x + A1zs*Zs + A1zns*Zns + (not A1s because time's dependent)
//fill matrix: x'=A1x * mx + A1zs * mZs + A1zns * mZns +A1s; 
void linearSystem::computedxdt(){
  if (mx.size() == 0){
    ACE_MESSAGE("No Dynamic\n");
    return;
  }
  
  buildABCDs();
  if(!mA)
    return;
  try{
    cout<<"mA\n"<<(*mA);
    if (ACE_MAT_TYPE == SPARSE){
      //*aaux=*mA;
      mPxAux->set(*mA);
      cout<<"mPxAux\n"<<(*mPxAux);
      mPxAux->PLUInverseInPlace();
      cout<<"mPxAux-1\n"<<(*mPxAux);
      mA->set(*mPxAux);
    }else{
      mA->PLUInverseInPlace();
    }
    cout<<"mA-1\n"<<(*mA);
  }
  catch(SiconosException e)
    {
      std::cout << e.report() << endl;
      ACE_INTERNAL_ERROR("linearSystem::computedxdt");
    }
  catch(...)
    {
      std::cout << "Exception caught." << endl;
      ACE_INTERNAL_ERROR("linearSystem::computedxdt");
    }
  cout<<"inv A:\n-----\n";
  cout<<(*mA);

  try{
    if (mA)
      //(*mA1x)=prod(*mA,*mB);
      ACEprod(*mA,*mB,*mA1x,true);
    if (mC)
      //(*mA1zs)=prod(*mA,*mC);
      ACEprod(*mA,*mC,*mA1zs,true);
    if (mD)
      //(*mA1zns)=prod(*mA,*mD);
      ACEprod(*mA,*mD,*mA1zns,true);
    
    // (*mA1s)=prod(*mA,*ms);
  }
  catch(SiconosException e)
    {
      std::cout << e.report() << endl;
      ACE_INTERNAL_ERROR("linearSystem::computedxdt");
    }
  catch(...)
    {
      std::cout << "Exception caught." << endl;
      ACE_INTERNAL_ERROR("linearSystem::computedxdt");
    }
  

}
void linearSystem::getlinefromdxdt(int line, ACE_DOUBLE * coefs){
  if (!mA)
    return;
  ACE_CHECK_IERROR(coefs && line < mNbDynEquations && line >=0,"linearSystem::getlinefromdxdt");
  int idStart =0;
  int i=0;
  for (i =idStart; i < mDimx;i++){
    coefs[i]= mA1x->getValue(line,i-idStart);
  }
  idStart+=mDimx;
  coefs[idStart]=0;//Because V0
  for (i =idStart+1; i < idStart+mDimzs;i++){
    coefs[i]= mA1zs->getValue(line,i-idStart-1);
  }
  idStart+=mDimzs;
  for (i =idStart; i < idStart+mDimzns;i++){
    coefs[i]= mA1zns->getValue(line,i-idStart);
  }
  idStart+=mDimzns;
  ACE_CHECK_IERROR(idStart+mDimx == mNbUnknowns,"linearSystem::getlinefromdxdt idStart +mDimx == mNbUnknowns");
  coefs[idStart]=mA1s->getValue(line);

}
//Ax'=Bx+CZs+DZns+s
void linearSystem::buildABCDs(){
  int istart=0;
  if (mDimx==0)
    return;
  ACE_CHECK_IERROR( mDimx == mNbDynEquations,"linearSystem::buildABCDs : mDimx != nbDynEquations");
  //BUILD A
  if (mA)
    extractDynBockInMat(mA,istart,istart+mDimx);

  istart+=mDimx;
  //BUILD B
  if (mB)
    extractDynBockInMat(mB,istart,istart+mDimx);
  istart+=mDimx;

  //BUILD C
  if (mC)
    extractDynBockInMat(mC,istart+1,istart+mDimzs);//because V0
  istart+=mDimzs;

  //BUILD D
  if (mD)
    extractDynBockInMat(mD,istart,istart+mDimzns);    
  istart+=mDimzns;

  //BUILD s
  ACE_CHECK_IERROR(istart == mNbUnknowns,"linearSystem::buildABCDs istart == mNbUnknowns");

}
//copy in m the double contains in dynamique equation from IndexBegin to IndexEnd
void linearSystem::extractDynBockInMat(aceMatrix * m, int IndexBegin, int IndexEnd){
  ACE_CHECK_IERROR(m,"linearSystem::extractDynBockInMat m null");
  ACE_CHECK_IERROR(IndexBegin>=0,"linearSystem::extractDynBockInMat IndexBegin");
  ACE_CHECK_IERROR(IndexEnd>=IndexBegin && IndexEnd < mNbUnknowns+2 ,"linearSystem::extractDynBockInMat IndexEnd");
  
  
    int nbInd = mIND.size();
    int i,j=0;
    int line=0;
    for (i=0;i<nbInd;i++){
      ACE_DOUBLE * coefs=mIND[i]->mCoefs;
      line = mIND[i]->mLine;
      ACE_CHECK_IERROR(line>=0 && coefs,"linearSystem::extractDynBockInMat Ind");
      for(j=IndexBegin;j<IndexEnd;j++){
	m->setValueIfNotNull(line,j-IndexBegin,coefs[j]);
      }
    }
    int nbCap = mCAP.size();
    for (i=0;i<nbCap;i++){
      ACE_DOUBLE * coefs=mCAP[i]->mCoefs;    
      line = mCAP[i]->mLine;
      ACE_CHECK_IERROR(line>=0 && coefs,"linearSystem::extractDynBockInMat  Cap");
      for(j=IndexBegin;j<IndexEnd;j++){
	m->setValueIfNotNull(line,j-IndexBegin,coefs[j]);
      }
    }
    int nbKcl = mKCL.size();
    for (i=0;i<nbKcl;i++){
      if (mKCL[i]->mIsDyn){
	ACE_DOUBLE * coefs=mKCL[i]->mCoefs;    
	line = mKCL[i]->mLine;
	ACE_CHECK_IERROR(line>=0 && coefs,"linearSystem::extractDynBockInMat  Kcl");
	for(j=IndexBegin;j<IndexEnd;j++){
	  m->setValueIfNotNull(line,j-IndexBegin,coefs[j]);
	}
      }
    }
}


//copy in m the double contains in dynamique equation from IndexBegin to IndexEnd
void linearSystem::extractDynBockInVect(aceVector * V){
  ACE_CHECK_IERROR(V,"linearSystem::extractDynBockInMat v null");
  
  
    int nbInd = mIND.size();
    int i,j=0;
    int line=0;
    j=mNbUnknowns;
    for (i=0;i<nbInd;i++){
      ACE_DOUBLE * coefs=mIND[i]->mCoefs;
      line = mIND[i]->mLine;
      ACE_CHECK_IERROR(line>=0 && coefs,"linearSystem::extractDynBockInMat Ind");
      V->setValueIfNotNull(line,coefs[j]);
    }
    int nbCap = mCAP.size();
    for (i=0;i<nbCap;i++){
      ACE_DOUBLE * coefs=mCAP[i]->mCoefs;    
      line = mCAP[i]->mLine;
      ACE_CHECK_IERROR(line>=0 && coefs,"linearSystem::extractDynBockInMat  Cap");
      V->setValueIfNotNull(line,coefs[j]);
    }
    int nbKcl = mKCL.size();
    for (i=0;i<nbKcl;i++){
      if (mKCL[i]->mIsDyn){
	ACE_DOUBLE * coefs=mKCL[i]->mCoefs;    
	line = mKCL[i]->mLine;
	ACE_CHECK_IERROR(line>=0 && coefs,"linearSystem::extractDynBockInMat  Kcl");
	V->setValueIfNotNull(line,coefs[j]);
      }
    }
}


///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////SOURCES must be call each step
void linearSystem::extractDynamicSystemSource(){
  if (mDimx){
    extractDynBockInVect(ms);
    //(*mA1s)=prod(*mA,*ms);
    ACEprod(*mA,*ms,*mA1s,true);
  }

}
void linearSystem::extractInteractionSource(){
    extractNonDynBockInVect(mB1s);
}

void linearSystem::extractSources(){
  extractInteractionSource();
  extractDynamicSystemSource();
}
void linearSystem::computeDynamicSystemSource(){
  if (mDimx && mDimLambda){
    //(*mA2s)=(*mA1s)+prod(*mA1zns,*mC1s);
    ACEprod(*mA1zns,*mC1s,*mA2s,false);
    //prod(*mA1zns,*mC1s,*mA2s);
    //(*mA2s)+=(*mA1s);
    
  }
}
void linearSystem::computeInteractionSource(){
  if (mDimx && mDimLambda){
    
    ACEprod(*mB1zns,*mC1s,*mB2s,false);
    //prod(*mB1zns,*mC1s,*mB2s);
    //(*mB2s)+=(*mB1s);

    ACEprod(*mD1zns,*mC1s,*mD2s,false);
    //prod(*mD1zns,*mC1s,*mD2s);
    //(*mD2s)+=(*mD1s);
  }else if (mDimLambda){
    //(*mB2s)=(*mB1s)+prod(*mB1zns,*mC1s);
    ACEprod(*mB1zns,*mC1s,*mB2s,false);
    //(*mD2s)=(*mD1s)+prod(*mD1zns,*mC1s);
    ACEprod(*mD1zns,*mC1s,*mD2s,false);
  }
}

void linearSystem::ExtractAndCompute2Sources(){
  if (mDimx)
    *mA2sti=*mA2s;
  extractSources();
//   cout<<"ExtractAndCompute2Sources input:\n";
//   cout<<"mA1zns\n";
//   cout<<*mA1zns;
//   cout<<"mA2s\n";
//   cout<<*mA2s;
//   cout<<"mC1s\n";
//   cout<<*mC1s;
  computeDynamicSystemSource();
  computeInteractionSource();

//   if (mDimx && mDimLambda){
//     //(*mA2s)=(*mA1s)+prod(*mA1zns,*mC1s);
//     ACEprod(*mA1zns,*mC1s,*mA2s,false);
//     //prod(*mA1zns,*mC1s,*mA2s);
//     //(*mA2s)+=(*mA1s);
    
//     ACEprod(*mB1zns,*mC1s,*mB2s,false);
//     //prod(*mB1zns,*mC1s,*mB2s);
//     //(*mB2s)+=(*mB1s);

//     ACEprod(*mD1zns,*mC1s,*mD2s,false);
//     //prod(*mD1zns,*mC1s,*mD2s);
//     //(*mD2s)+=(*mD1s);
//   }else if (mDimx){
//     (*mA2s)=(*mA1s);
//     (*mB2s)=(*mB1s);    
//   }else if (mDimLambda){
//     //(*mB2s)=(*mB1s)+prod(*mB1zns,*mC1s);
//     ACEprod(*mB1zns,*mC1s,*mB2s,false);
//     //(*mD2s)=(*mD1s)+prod(*mD1zns,*mC1s);
//     ACEprod(*mD1zns,*mC1s,*mD2s,false);
//   }
//   cout<<"ExtractAndCompute2Sources output:\n";
//   cout<<"mA2s\n";
//   cout<<*mA2s;
//   cout<<"mB2s\n";
//   cout<<*mB2s;
//   cout<<"mD2s\n";
//   cout<<*mD2s;
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////LINEAR SYSTEM : 0 = B1x*x + B1zs*Zs + B1zns*Zns + (not B1s because time dependent)
void linearSystem::buildLinearSystem(){
  
  int istart = mDimx;
  //BUILD B1x
  if (mB1x)
    extractNonDynBockInMat(mB1x,istart,istart+mDimx);
  istart+=mDimx;
  
  //BUILD B1zs
  if (mB1zs)
    extractNonDynBockInMat(mB1zs,istart+1,istart+mDimzs);
  istart+=mDimzs;

  //BUILD B1zns
  if (mB1zns)
    extractNonDynBockInMat(mB1zns,istart,istart+mDimzns);    
  istart+=mDimzns;
  ACE_CHECK_IERROR(istart == mNbUnknowns,"linearSystem::buildLinearSystem istart == mNbUnknowns");
}

//copy in m the double contains in none dynamique equation from IndexBegin to IndexEnd
void linearSystem::extractNonDynBockInVect(aceVector * V){
  ACE_CHECK_IERROR(V,"linearSystem::extractNonDynBockInMat v null");
  
  int line=0;
  int i,j;
  j=mNbUnknowns;
  int nbKcl = mKCL.size();
  for (i=1;i<nbKcl;i++){
    if (!mKCL[i]->mIsDyn){
      ACE_DOUBLE * coefs=mKCL[i]->mCoefs;    
      ACE_CHECK_IERROR(coefs,"linearSystem::extractNonDynBockInMat  Kcl");
      V->setValueIfNotNull(line,coefs[j]);
      line++;
    }
  }
  int nVD = mVD.size();
  for (i=0;i<nVD;i++){
      ACE_DOUBLE * coefs=mVD[i]->mCoefs;    
      ACE_CHECK_IERROR(coefs,"linearSystem::extractNonDynBockInMat  VD");
      V->setValueIfNotNull(line,coefs[j]);
      line++;
  }

  int nTEN = mTEN.size();
  for (i=0;i<nTEN;i++){
      ACE_DOUBLE * coefs=mTEN[i]->mCoefs;    
      ACE_CHECK_IERROR(coefs,"linearSystem::extractNonDynBockInMat  TEN");
      V->setValueIfNotNull(line,coefs[j]);
      line++;
  }

  ACE_CHECK_IERROR(mNbNonDynEquations == line,"linearSystem::extractNonDynBockInMat number equation!");
}
//copy in m the double contains in none dynamique equation from IndexBegin to IndexEnd
void linearSystem::extractNonDynBockInMat(aceMatrix * m, int IndexBegin, int IndexEnd){
  ACE_CHECK_IERROR(m,"linearSystem::extractNonDynBockInMat m null");
  ACE_CHECK_IERROR(IndexBegin>=mDimx,"linearSystem::extractNonDynBockInMat IndexBegin");
  ACE_CHECK_IERROR(IndexEnd>=IndexBegin && IndexEnd < mNbUnknowns+2 ,"linearSystem::extractNonDynBockInMat IndexEnd");
  
  int line=0;
  int i,j;

  int nbKcl = mKCL.size();
  for (i=1;i<nbKcl;i++){
    if (!mKCL[i]->mIsDyn){
      ACE_DOUBLE * coefs=mKCL[i]->mCoefs;    
      ACE_CHECK_IERROR(coefs,"linearSystem::extractNonDynBockInMat  Kcl");
      for(j=IndexBegin;j<IndexEnd;j++){
	m->setValueIfNotNull(line,j-IndexBegin,coefs[j]);
      }
      line++;
    }
  }
  int nVD = mVD.size();
  for (i=0;i<nVD;i++){
      ACE_DOUBLE * coefs=mVD[i]->mCoefs;    
      ACE_CHECK_IERROR(coefs,"linearSystem::extractNonDynBockInMat  VD");
      for(j=IndexBegin;j<IndexEnd;j++){
	m->setValueIfNotNull(line,j-IndexBegin,coefs[j]);
      }
      line++;
  }

  int nTEN = mTEN.size();
  for (i=0;i<nTEN;i++){
      ACE_DOUBLE * coefs=mTEN[i]->mCoefs;    
      ACE_CHECK_IERROR(coefs,"linearSystem::extractNonDynBockInMat  TEN");
      for(j=IndexBegin;j<IndexEnd;j++){
	m->setValueIfNotNull(line,j-IndexBegin,coefs[j]);
      }
      line++;
  }

  ACE_CHECK_IERROR(mNbNonDynEquations == line,"linearSystem::extractNonDynBockInMat number equation!");
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////PRINT FONCTIONS
 //
  //R=A1zns*C1l
  //x'=A2x*x + A2zs*Zs + R*lambda+A2s
  //0=B2x*x + B2zs*Zs + B2l*lambda + B2s
  //Y=D2x*x + D2zs*Zs + D2l*lambda + D2s
void linearSystem::printSystem2(ostream& os){
  os<<"R=A1zns*C1l\nx'=A2x*x + A2zs*Zs + R*lambda+A2s\n0=B2x*x + B2zs*Zs + B2l*lambda + B2s\nY=D2x*x + D2zs*Zs + D2l*lambda + D2s\n";
  os<<"R:\n";
  if (mR)
    os<<(*mR);
  os<<"A2x:\n";
  if (mA2x)
    os<<(*mA2x);
  os<<"A2zs:\n";
  if (mA2zs)
    os<<(*mA2zs);
  os<<"A2s:\n";
  if (mA2s)
    os<<(*mA2s);
  os<<"B2x:\n";
  if (mB2x)
    os<<(*mB2x);
  os<<"B2zs:\n";
  if (mB2zs)
    os<<(*mB2zs);
  os<<"B2l:\n";
  if (mB2l)
    os<<(*mB2l);
  os<<"B2s:\n";
  if (mB2s)
    os<<(*mB2s);
  os<<"D2x:\n";
  if (mD2x)
    os<<(*mD2x);
  os<<"D2zs:\n";
  if (mD2zs)
    os<<(*mD2zs);
  os<<"D2l:\n";
  if (mD2l)
    os<<(*mD2l);
  os<<"D2s:\n";
  if (mD2s)
    os<<(*mD2s);
}
void linearSystem::printStep(ostream& os){
  int i;
  bool printALL = true;
  if (printALL){
    os << "xt("<<mStepCmp*mH<<")\t";
    if (mxti)
      for (i=0;i<mDimx;i++)
	os << mxti->getValue(i)<<"\t";
    os << "zs("<<mStepCmp*mH<<")\t";
    if (mzsti)
      for (i=0;i<mDimzs-1;i++)
	os << mzsti->getValue(i)<<"\t";
  
    os << "zns("<<mStepCmp*mH<<")\t";
    if (mznsti)
      for (i=0;i<mDimzns;i++)
	os << mznsti->getValue(i)<<"\t";
    os<<"\n";
  }else{
    if (mStepCmp%mLogPrint==0){
      os <<mStepCmp*mH<<"\t"<<mxti->getValue(4)<<"\n";
    }
  }
}

void linearSystem::printEquations(ostream& os ){
  int i,n =0;
  
  printf("--->linearSystem with %d equations whose %d dynamic equations.\n",mNbEquations,mNbDynEquations);
  printf("x\n");
  for (i=0; i<mDimx; i++)
    mx[i]->print();
  printf("\nZs\n");
  for (i=0; i<mDimzs; i++)
    mZs[i]->print();
  printf("\nZns\n");
  for (i=0; i<mDimzns; i++)
    mZns[i]->print();
  printf("\n---------------------------------------\nequation");
  for (i=0; i<mDimx; i++)
    mx[i]->printdev();
  for (i=0; i<mDimx; i++)
    mx[i]->print();
  for (i=0; i<mDimzs; i++)
    mZs[i]->print();
  for (i=0; i<mDimzns; i++)
    mZns[i]->print();
  printf("\n");
  
  //print dyn equation
  n=mKCL.size();
  for (i=0;i<n;i++)
    if (mKCL[i]->mIsDyn)
      mKCL[i]->print();
  n=mCAP.size();
  for (i=0;i<n;i++)
    mCAP[i]->print();  
  n=mIND.size();
  for (i=0;i<n;i++)
    mIND[i]->print();

  //non dyn equation
  n=mKCL.size();
  for (i=0;i<n;i++)
    if (!mKCL[i]->mIsDyn)
      mKCL[i]->print();
  n=mVD.size();
  for (i=0;i<n;i++)
    mVD[i]->print();
  n=mTEN.size();
  for (i=0;i<n;i++)
    mTEN[i]->print();
}

void linearSystem::printABCDs(ostream& os ){
  os <<"Ax'=Bx+CZs+DZns+s\n";
  os <<"-----------------\n";
  os <<"A:\n";
  if (mA)
    os << (*mA);
  os <<"B:\n";
  if (mB)
    os<< (*mB);
  os<<"C:\n";
  if (mC)
    os <<(*mC);
  os<<"D:\n";
  if (mD)
    os <<(*mD);
  os <<"s:\n";
  if (ms)
    os <<(*ms);
 }
void linearSystem::printA1(ostream& os ){
  os <<"system x'=A1x+A1zs+A1zns+s\n";
  os <<"A1x:\n";
  if (mA1x)
    os <<(*mA1x);
  os <<"A1zs:\n";
  if (mA1zs)
    os <<(*mA1zs);
  os <<"A1zns:\n";
  if (mA1zns)
    os <<(*mA1zns);
  os <<"A1s:\n";
  if (mA1s)
    os << (*mA1s);
 }
void linearSystem::printB1(ostream& os ){
  os <<"0 = B1x*x + B1zs*Zs + B1zns*Zns + B1s\n";
  os <<"B1x:\n";
  if (mB1x)
    os <<(*mB1x);
  os <<"B1zs:\n";
  if (mB1zs)
    os <<(*mB1zs);
  os <<"B1zns:\n";
  if (mB1zns)
    os <<(*mB1zns);
  os <<"B1s:\n";
  if (mB1s)
    os << (*mB1s);
 }
void linearSystem::printC1(ostream& os ){
  os <<"Zns = C1x*x + C1s*Zs + C1l*lamdba + C1s\n";
  os <<"C1x:\n";
  if (mC1x)
    os <<(*mC1x);
  os <<"C1zs:\n";
  if (mC1zs)
    os <<(*mC1zs);
  os <<"C1l:\n";
  if (mC1l)
    os <<(*mC1l);
  os <<"C1s:\n";
  if (mC1s)
    os << (*mC1s);
 }
void linearSystem::printD1(ostream& os ){
  os <<"Y = D1x*x + D1s*Zs + D1ns*Zns + D1l*lambda +D1s\n";
  os <<"D1x:\n";
  if (mD1x)
    os <<(*mD1x);
  os <<"D1zs:\n";
  if (mD1zs)
    os <<(*mD1zs);
  os <<"D1zns:\n";
  if (mD1zns)
    os <<(*mD1zns);
  os <<"D1l:\n";
  if (mD1l)
    os <<(*mD1l);
  os <<"D1s:\n";
  if (mD1s)
    os << (*mD1s);
 }
void linearSystem::printDiscretisation(ostream& os){
  os <<"W(I - h*ThetaA2x)"<<endl;
  if (mW)
    os << mW;
  os <<"xfree"<<endl;
    
  os <<"0=Qfree + B3zs * Zs + B3l * Lambda"<<endl;
  os <<"Qfree"<<endl;
  if (mQfree)
    os << (*mQfree);
  os <<"B3zs"<<endl;
  if (mB3zs)
    os << (*mB3zs);
  os <<"B3l"<<endl;
  if (mB3l)
    os<<(*mB3l);
  os <<"Y=pfree + D3zs*Zs + D3l*Lambda"<<endl;
  os << "pfree\n";
  if (mPfree)
    os<<(*mPfree);
  os <<"D3zs"<<endl;
  if (mD3zs)
    os <<(*mD3zs);
  os <<"D3l"<<endl;
  if (mD3l)
    os <<(*mD3l);
}

void linearSystem::printSystemInTabFile(char * file){
  ofstream pout(file);
  printA1(pout);
  printB1(pout);
  printC1(pout);
  printD1(pout);
  printSystem2(pout);
  
  pout.close();
}
