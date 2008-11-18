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

  //mW=0;
  for (int i=0; i<ACE_NB_ADAPT_STEP+1;i++){
    mhR[i]=0;
    mD3l[i]=0;
    mD3zs[i]=0;
    mB3l[i]=0;
    mW[i]=0;
    mHThetaWA2zs[i]=0;
    mHWR[i]=0;
    mHThetaA2zs[i]=0;
    mD2xW[i]=0;
    mB2xW[i]=0;

  }
  //  mB3zs=0;
  mPfree=0;
  mPAux=0;
  mPxAux=0;
  mQfree=0;
  mMLCP=0;

  mTheta = 0.5;
  mThetap = 0.5;
  mH = 1;
  mHori=1;
  mTstart=0;
  mTstop=0;
  mTcurrent=0;
  //adaptive time stepping
  mUseAdaptiveTimeStepping=false;
  mAlphaMax=2;
  mAlphaMin=0.0001;
  mAlpha=0.9;
  mAdaptCmp=-1;
  mxtiprev=0;
  mzstiprev=0;
  mxticurrent=0;
  mzsticurrent=0;
  mzst_inter=0;
  mxt1=0;
  mzst1=0;
  mxbuf=0;
  mzsbuf=0;

  mNbToSmall=0;
  mNbToBig=0;
  mSerrorX=0;
  mSerrorZs=0;
  mNbBacktrack=0;



  mLogFrequency=0;
  mLogPrint=0;
  mPourMille=0;
  mSommeError=0;
  mCoef=1;

  //
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
    for (int i=0; i<ACE_NB_ADAPT_STEP+1;i++){
      mD2xW[i]=new aceMatrix(mDimLambda,mDimx,ACE_MAT_TYPE);
      mHWR[i]=new aceMatrix(mDimx,mDimLambda,ACE_MAT_TYPE);
    }
  }
  if (mNbNonDynEquations && mDimx)
    for (int i=0; i<ACE_NB_ADAPT_STEP+1;i++)
      mB2xW[i]=new aceMatrix(mNbNonDynEquations,mDimx,ACE_MAT_TYPE);
  if (mDimx)
    for (int i=0; i<ACE_NB_ADAPT_STEP+1;i++)
      mHThetaWA2zs[i]=new aceMatrix(mDimx,mDimzs-1,ACE_MAT_TYPE);

  
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
  for (int i=0; i <ACE_NB_ADAPT_STEP+1; i++)
    if(mHThetaA2zs[i])
      delete mHThetaA2zs[i];
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
    for(int i=0; i <ACE_NB_ADAPT_STEP+1; i++)
      mhR[i] = new aceMatrix(mDimx,mDimLambda,ACE_MAT_TYPE);
    mMatBuf2 = new aceMatrix(mDimx,mDimLambda,ACE_MAT_TYPE);
  }

}
void linearSystem::set2matrix(){

  mA2x=mA1x;
  mA2zs=mA1zs;
  for(int i=0; i <ACE_NB_ADAPT_STEP+1; i++)
    mHThetaA2zs[i] = new aceMatrix(mDimx,mDimzs-1,ACE_MAT_TYPE);
  
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

  if (mDimLambda){
    //(*mD2s)=(*mD1s)+prod(*mD1zns,*mC1s);
    ACEprod(*mD1zns,*mC1s,*mD2s,false);
  }

  
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
  for (int i=0; i <ACE_NB_ADAPT_STEP+1; i++)
    if (mhR[i])
      delete mhR[i];
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
    mxtiprev=new aceVector(mDimx,ACE_MAT_TYPE);
    mxticurrent=new aceVector(mDimx,ACE_MAT_TYPE);
    mxt1=new aceVector(mDimx,ACE_MAT_TYPE);
    mxbuf=new aceVector(mDimx,ACE_MAT_TYPE);
    

  }
  mzsti = new aceVector(mDimzs-1,ACE_MAT_TYPE);
  mzsticurrent=new aceVector(mDimzs-1,ACE_MAT_TYPE);
  mzstiprev=new aceVector(mDimzs-1,ACE_MAT_TYPE);
  mzst_inter=new aceVector(mDimzs-1,ACE_MAT_TYPE);
  mzst1=new aceVector(mDimzs-1,ACE_MAT_TYPE);
  mzsbuf=new aceVector(mDimzs-1,ACE_MAT_TYPE);
}

void linearSystem::allocDiscretisation(){
  if (mxfree)
    return;
  allocForInitialValue();
  if(mDimx){
    mxfree = new aceVector(mDimx,ACE_MAT_TYPE);
    for (int i=0;i<ACE_NB_ADAPT_STEP+1;i++){
      mW[i] = new aceMatrix(mDimx,mDimx,ACE_MAT_TYPE);
      mW[i]->zero();
    }
  }
  
  if(mDimzns)
    mznsti = new aceVector(mDimzns,ACE_MAT_TYPE);
  if (mB2zs){
    for (int i=0;i<ACE_NB_ADAPT_STEP+1;i++)
      mB3zs[i] = new aceMatrix(mNbNonDynEquations,mDimzs-1,ACE_MAT_TYPE);
  }
  if (mNbNonDynEquations && mDimLambda)
    for (int i=0;i<ACE_NB_ADAPT_STEP+1;i++)
      mB3l[i]=new aceMatrix(mNbNonDynEquations,mDimLambda,ACE_MAT_TYPE);
  if (mNbNonDynEquations)
    mQfree = new aceVector(mNbNonDynEquations,ACE_MAT_TYPE);

  if (mDimLambda)
    for (int i=0;i<ACE_NB_ADAPT_STEP+1;i++)
      mD3zs[i] = new aceMatrix(mDimLambda,mDimzs-1,ACE_MAT_TYPE);
  if (mDimLambda)
    for (int i=0;i<ACE_NB_ADAPT_STEP+1;i++)
      mD3l[i] = new aceMatrix(mDimLambda,mDimLambda,ACE_MAT_TYPE);
  if (mDimLambda){
    mPfree = new aceVector(mDimLambda,ACE_MAT_TYPE);
    mPAux = new aceVector(mDimLambda,ACE_MAT_TYPE);

  }
}
void linearSystem::freeForInitialValue(){
  if (mxti)
    delete mxti;
  if (mxtiprev)
    delete mxtiprev;
  if (mxticurrent)
    delete mxticurrent;
  if (mzsti)
    delete mzsti;
  if (mzst_inter)
    delete mzst_inter;
  if (mzst1)
    delete mzst1;
  if (mzsbuf)
    delete mzsbuf;
  if (mxt1)
    delete mxt1;
  if (mxbuf)
    delete mxbuf;
}
void linearSystem::freeDiscretisation(){
  freeForInitialValue();
  if (mxfree)
    delete mxfree;
  if (mznsti)
    delete mznsti;
  for (int i=0; i <ACE_NB_ADAPT_STEP+1; i++)
    if(mW[i])
      delete mW[i];
  for (int i=0; i <ACE_NB_ADAPT_STEP+1; i++)
    if (mD3l[i])
      delete mD3l[i];
  for (int i=0; i <ACE_NB_ADAPT_STEP+1; i++)
    if (mD3zs[i])
      delete mD3zs[i];
  for (int i=0; i <ACE_NB_ADAPT_STEP+1; i++)
    if (mB3l[i])
      delete mB3l[i];
  for (int i=0; i <ACE_NB_ADAPT_STEP+1; i++)
    if (mB3zs[i])
      delete mB3zs[i];
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
  for (int i=0; i <ACE_NB_ADAPT_STEP+1; i++){
    if (mD2xW[i])
      delete mD2xW[i];
    if (mB2xW[i])
      delete mB2xW[i];
    if (mHThetaWA2zs[i])
      delete mHThetaWA2zs[i];
    if (mHWR[i])
      delete mHWR[i];
  }

}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////DISCRETISATION
void linearSystem::readInitialValue(){
  int i;
  double aux;
  int nbGuess;
  unsigned long guess;
  int useIc;

  cout << "read init value\n";
  try{

    ParserGetTransValues(&mH,&mTstop,&mTstart);
    mH = mH/mCoef;
    mHori=mH;
#ifdef PRE_COMPUTE_ADAPTIVE
    mMaxHori = mHori*(1<<(ACE_NB_ADAPT_STEP-1));
#else
    mMaxHori = mHori*1000;
#endif
    long long b1 = llrint(mTstop/mH);;
    long b2 = lrint(mTstop/mH);;
    mStepNumber = b2;

    
      
    for (i=0;i<mDimzs-1;i++){
      mzsti->setValueIfNotNull(i,0);
    }
    
    for (i=0;i<mDimx;i++){
      mxti->setValueIfNotNull(i,0);
    }
    //DIODEBRIDGE :     mxti->setValueIfNotNull(1,10.0/2000);

    ParserInitICvalue();
    while(ParserGetICvalue(&i,&useIc,&aux)){
      if (i>0){
	mzsti->setValueIfNotNull(i-1,aux);
	cout<<"set value from netlist :v_"<<i<<"="<<aux<<endl;
      }
    }
    //    printStep(*mSimuStream,mzsti);
    mzsti->setValueIfNotNull(2,-5e-2);
    /*for each x of type TENSION, compute the initial value with V init*/
    for(i=0; i <(int) mx.size(); i++){
      unknown* u = mx[i];
      if (u->mType == ACE_TYPE_U ){
	aux =  (u->mComponent->mNodeNeg>0)? mzsti->getValue(u->mComponent->mNodeNeg-1):0;
	aux -= (u->mComponent->mNodePos>0)? mzsti->getValue(u->mComponent->mNodePos-1):0;
	mxti->setValueIfNotNull(i,aux);
      }
    }

  }
  catch(...)
    {
      std::cout << "Exception caught." << endl;
      ACE_INTERNAL_ERROR("linearSystem::readInitialValue");
    }

}

void linearSystem::buildMLCP(){
  if (!mMLCP)
    mMLCP = new mlcp(mDimLambda,mNbNonDynEquations,ACE_SOLVER_TYPE);
  for (int i=0;i<ACE_NB_ADAPT_STEP+1;i++){

#ifdef PRE_COMPUTE_ADAPTIVE
    if (ACE_WITH_ADAPTATIVE_TIME_STEPPING)
	mH = mHori*(1<<i);
#endif

    
    try{
      //ACE_CHECK_IERROR(mDimx >0 && mDimLambda >0,"linearSystem::initSimu case no x or no lambda not yet implemented");
      if (mDimx){
	(mW[i])->zero();
	for(int j = 0; j < mDimx;j++)	mW[i]->setValue(j,j,1.0);
	*mW[i] -=  mH*mTheta*(*mA2x);
	if (ACE_MAT_TYPE == SPARSE){
	  mPxAux->set(*mW[i]);
	  mPxAux->PLUInverseInPlace();
	  mW[i]->set(*mPxAux);
	}else{
	  mW[i]->PLUInverseInPlace();
	}
      
	//*mB3zs = *mB2zs + mH*mTheta*prod(*mB2x,prod(*mW,*mA2zs)) ;
	ACEprod(*mW[i],*mA2zs,*mMatBuf1,true);
	ACEprod(*mB2x,*mMatBuf1,*mB3zs[i],true);
	scal(mH*mTheta,*mB3zs[i],*mB3zs[i]);
	scal(mH*mTheta,*mA2zs,*mHThetaA2zs[i]);
	*mB3zs[i]+=*mB2zs;
      
	if (mDimLambda){
	  //*mD3zs = *mD2zs+mH*mTheta*prod(*mD2x,prod(*mW,*mA2zs));
	  ACEprod(*mD2x,*mMatBuf1,*mD3zs[i],true);
	  scal(mH*mTheta,*mD3zs[i],*mD3zs[i]);
	  *mD3zs[i]+=*mD2zs;

      
	  //*mB3l = mH*prod(*mB2x,prod(*mW,*mR))+*mB2l;
	  ACEprod(*mW[i],*mR,*mMatBuf2,true);
	  ACEprod(*mB2x,*mMatBuf2,*mB3l[i],true);
	  scal(mH,*mB3l[i],*mB3l[i]);
	  *mB3l[i]+=*mB2l;
	     
	  //*mD3l = mH*prod(*mD2x,prod(*mW,*mR))+*mD2l;
	  ACEprod(*mD2x,*mMatBuf2,*mD3l[i],true);
	  scal(mH,*mD3l[i],*mD3l[i]);
	  *mD3l[i]+=*mD2l;
	}
      }else{
	*mD3zs[i] = *mD2zs;
	*mB3zs[i] = *mB2zs;
	if (mDimLambda){
	  *mB3l[i] = *mB2l;
	  *mD3l[i] = *mD2l;
	}
      }

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
      ACEprod(*mD2x,*mW[i],*mD2xW[i],true);
    }
    if (mDimx){
      //*mB2xW=prod(*mB2x,*mW);
      ACEprod(*mB2x,*mW[i],*mB2xW[i],true);
      //*mHThetaWA2zs=mH*mTheta*prod(*mW,*mA2zs);
      ACEprod(*mW[i],*mA2zs,*mHThetaWA2zs[i],true);
      scal(mH*mTheta,*mHThetaWA2zs[i],*mHThetaWA2zs[i]);
    }
    //*mHWR=mH*prod(*mW,*mR);
    if (mDimLambda && mDimx){
      ACEprod(*mW[i],*mR,*mHWR[i],true);
      scal(mH,*mHWR[i],*mHWR[i]);
    }
  }
}
void linearSystem::fillMLCP(){
  *(mMLCP->mM11) = *mD3l[ACE_CUR_STEP];
  *(mMLCP->mM12) = *mD3zs[ACE_CUR_STEP];
  *(mMLCP->mM21) = *mB3l[ACE_CUR_STEP];
  *(mMLCP->mM22) = *mB3zs[ACE_CUR_STEP];
  mMLCP->update();

}
void linearSystem::initSimu(ofstream * pstream){
  mSimuStream = pstream;
  mSommeError=0;
  mTcurrent=0;
  allocDiscretisation();
  readInitialValue();
  //compteur for the adaptive time stepping.
  mNbToSmall=0;
  mNbToBig=0;
  mSerrorX=0;
  mSerrorZs=0;
  mNbBacktrack=0;
  mAllStepCmp=0;


  mStepCmp=0;
  mTcurrent=0;
  mLogFrequency =  mStepNumber/1000;
  ACE_DOUBLE aux = log(mStepNumber)/log(10);
  mLogPrint = 1;
  mStepNumber/10000;
  if (mLogPrint==0) mLogPrint=1;
  if (mLogFrequency==0) mLogFrequency = 1;
  mPourMille=0;

  buildMLCP();
  ACE_CUR_STEP=0;
  mH=mHori;
  fillMLCP();
  mMLCP->initSolver();
  for (int i=0; i<ACE_NB_ADAPT_STEP+1;i++)
    ACE_CMP_ADAT[ACE_CUR_STEP]=0;

  
}
double linearSystem::getCurrentTime(){
  return mTcurrent;
}
void linearSystem::preparStep(){
  computeAndAcceptStep();
  preparMLCP();
}

void linearSystem::preparMLCP(){
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
    ACEprod(*(mD2xW[ACE_CUR_STEP]),*mxfree,*mPfree,true);
    *mPfree+=*mD2s;
    //(*mQfree) = prod(*mB2xW,*mxfree) + (*mB2s);
    ACEprod(*mB2xW[ACE_CUR_STEP],*mxfree,*mQfree,true);
    *mQfree+=*mB2s;
    
  } else if(mDimx){//only x
    (*mxfree) = (*mxti) + mH*(1-mTheta)*prod(*mA2x,*mxti)+mH*(1-mTheta)*prod(*mA2zs,*mzsti)+mH*(*mA2s);
    (*mQfree) = prod(*mB2xW[ACE_CUR_STEP],*mxfree) + (*mB2s);
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
    //    scal(-1,*(mMLCP->mQ1),*(mMLCP->mQ1));
    //  *(mMLCP->mQ1)= *mPfree;  
    //*(mMLCP->mQ2)= *mQfree;
    (mMLCP->mQ2)->set(*(mQfree));
    //    scal(-1,*(mMLCP->mQ2),*(mMLCP->mQ2));
  }else{
      //  scal(-1,*mPfree,*(mMLCP->mQ1));
      //scal(-1,*mQfree,*(mMLCP->mQ2));
    if (mDimLambda)
      *(mMLCP->mQ1)= *mPfree;  
    *(mMLCP->mQ2)= *mQfree;
  }
//    cout<<"mMLCP->mQ1\n";
//    (mMLCP->mQ1)->display();
//    cout<<"mMLCP->mQ2\n";
//    (mMLCP->mQ2)->display();

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
  mAllStepCmp++;

  if (mTcurrent+mH >= mTstop){
    ACE_times[ACE_TIMER_LS_STEP].stop();
    return false;
  }
  int aux = mPourMille;
  ACE_DOUBLE aux1 = 1000*mTcurrent/mTstop;
  
  while (aux < 1000 && aux < aux1)
    aux++;
  if (aux != mPourMille){
    mPourMille= aux;
    printf("-->%d : %e to %e\n",mPourMille, mTcurrent,mTcurrent + mH);
  }
  mTcurrent += mH;
  ACE_CMP_ADAT[ACE_CUR_STEP]++;
  //  mMLCP->printInPutABCDab();
  bool res = mMLCP->solve();
  //  mMLCP->printOutPut();
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
	ACEprod(*mHWR[ACE_CUR_STEP],*mPAux,*mxti,true);
      }else{
	if (mDimLambda)
	  ACEprod(*mHWR[ACE_CUR_STEP],*(mMLCP->mZ1),*mxti,true);
	else
	  mxti->zero();
      }
      //cout<<"mHThetaWA2zs "<<(*mHThetaWA2zs);
      //cout<<"mW "<<(*mW);
      ACE_times[ACE_TIMER_PROD_MAT].start();
      ACEprod(*mHThetaWA2zs[ACE_CUR_STEP],*mzsti,*mxti,false);
      ACEprod(*mW[ACE_CUR_STEP],*mxfree,*mxti,false);
      ACE_times[ACE_TIMER_PROD_MAT].stop();
    }
    computeZnstiFromX_Zs();
  }else{
    ACE_GET_LOG_STREAM()<<"linearSystem::step number,"<< mStepCmp<<" solver failled!!!"<<endl;
  }
  ACE_times[ACE_TIMER_COMPUTE_VAR].stop();
  ACE_times[ACE_TIMER_LS_STEP].stop();
  return res;
}
void linearSystem::printLog(){
  ACE_GET_LOG_STREAM()<<"linearSystem::stopSimu, number of ALL steps: "<<mAllStepCmp<<endl;
  ACE_GET_LOG_STREAM()<<"linearSystem::stopSimu, adaptive nb backtrack: "<<mNbBacktrack<<endl;
  ACE_GET_LOG_STREAM()<<"linearSystem::stopSimu, adaptive to small step: "<<mNbToSmall<<endl;
  ACE_GET_LOG_STREAM()<<"linearSystem::stopSimu, adaptive to big step: "<<mNbToBig<<endl;
  ACE_GET_LOG_STREAM()<<"linearSystem::stopSimu, total error voltage step: "<<mSerrorZs<<endl;
  ACE_GET_LOG_STREAM()<<"linearSystem::stopSimu, total error X step: "<<mSerrorX<<endl;
  ACE_GET_LOG_STREAM()<<"linearSystem::stopSimu,  h, error somme, number of step, average: "<<mH<<" "<<mSommeError<<" "<<mStepCmp<<" "<<mSommeError/mStepCmp<<endl;
  ACE_GET_LOG_STREAM()<<"linearSystem::stopSimu,  log10(average), log10(h)"<<log10(mSommeError/mStepCmp)<<" "<<log10(mH)<<endl;
  if (ACE_WITH_ADAPTATIVE_TIME_STEPPING){
    ACE_GET_LOG1_STREAM()<<log10(mSommeError/mStepCmp)<<"\t"<<log10(ACE_RTOL_LOCAL)<<endl;
  }else{
    ACE_GET_LOG1_STREAM()<<log10(mSommeError/mStepCmp)<<"\t"<<log10(mH)<<endl;
  }
  for (int i=0;i < ACE_NB_ADAPT_STEP+1;i++)
    ACE_GET_LOG_STREAM()<<"linearSystem::stopSimu, step number : "<< i << " : "<< ACE_CMP_ADAT[i]<<endl;
 }
void linearSystem::stopSimu(){
  //  mMLCP->printGuess();
  printLog();
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
    if (ACE_MUET_LEVEL != ACE_MUET)
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
    if (ACE_MUET_LEVEL != ACE_MUET)
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
  if (ACE_MUET_LEVEL != ACE_MUET){
    cout<<"inv A:\n-----\n";
    cout<<(*mA);
  }

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
void linearSystem::updateDynamicSystemSource(){
  if (mDimx && mDimLambda){
    //(*mA2s)=(*mA1s)+prod(*mA1zns,*mC1s);
    ACEprod(*mA1zns,*mC1s,*mA2s,false);
    //prod(*mA1zns,*mC1s,*mA2s);
    //(*mA2s)+=(*mA1s);
    
  }
}
void linearSystem::updateInteractionSource(){
  
  if (mDimLambda){
    //(*mB2s)=(*mB1s)+prod(*mB1zns,*mC1s);
    ACEprod(*mB1zns,*mC1s,*mB2s,false);
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
  updateDynamicSystemSource();
  updateInteractionSource();

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
  // cout<<"mB2s\n";
  // cout<<*mB2s;
  // cout<<"mD2s\n";
  // cout<<*mD2s;
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
    if (ACE_MUET_LEVEL == ACE_MUET)
    return;

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
ACE_DOUBLE linearSystem::computeAnalyticSolution(ACE_DOUBLE t){
    ACE_DOUBLE R=1000;
  ACE_DOUBLE L=1e-2;
  ACE_DOUBLE C=1e-6;

  ACE_DOUBLE a = -t/(2.0*R*C);
  ACE_DOUBLE w = t*sqrt((1.0/(L*C)-(1.0/(4.0*R*R*C*C))));
  return fabs(10.0*exp(a)*cos(w));
  /*  ACE_DOUBLE R=100;
  ACE_DOUBLE C=1e-6;

  ACE_DOUBLE a = -5*exp(-t/(R*C));
  return a;*/

  
}
void printHeaderCSV(ostream& os){
  os<<".HEADER"<<endl;
  os<<"..NAMES"<<endl;
  os<<"Time,I(L1)"<<endl;
  os<<"..UNITS"<<endl;
  os<<"s,Current(A)"<<endl;
  os<<"..DATATYPES"<<endl;
  os<<"double,double"<<endl;
  os<<"..WAVEFORM_TYPES"<<endl;
  os<<",analog"<<endl;
  os<<"..AXIS_SPACING"<<endl;
  os<<"linear"<<endl;
  os<<"..FOLDER_NAME"<<endl;
  os<<"BuckConverter, TRAN"<<endl;
  os<<".DATA"<<endl;

}
void linearSystem::printStep(ostream& os,aceVector *pVzs){
//   int i;
   bool printALL = true;
   ACE_DOUBLE curDate = getCurrentTime();
   if (!printALL){
     if (curDate==0.0)
       printHeaderCSV(os);
     os <<curDate<<","<<mxti->getValue(4)<<"\n";
     return;
   }

   //     ACE_DOUBLE v = computeAnalyticSolution(getCurrentTime());

  dataPrint * pPrint;
  double aux;
  ParserInitPrintElem();
  os<<getCurrentTime();
  //  os <<"\t"<<mxti->getValue(4); // for buck with out inverter
  //os <<"\t"<<mxti->getValue(14); // for buck with inverter
//   return;
  while(ParserGetPrintElem((void**)&pPrint)){
    os<<"\t\t";
    aux = pVzs->getValue(pPrint->node1-1);
    if (pPrint->node2 >0)
      aux -=  pVzs->getValue(pPrint->node2-1);
    // mSommeError += fabs(v-aux);
    
      os << aux ;
      //      os << aux << "\t"<<v;
  }
  os<<endl;
}

void linearSystem::printEquations(ostream& os ){
  if (ACE_MUET_LEVEL == ACE_MUET)
    return;
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
  if (ACE_MUET_LEVEL == ACE_MUET)
    return;

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
  if (ACE_MUET_LEVEL == ACE_MUET)
    return;

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
  if (ACE_MUET_LEVEL == ACE_MUET)
    return;

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
  if (ACE_MUET_LEVEL == ACE_MUET)
    return;

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
  if (ACE_MUET_LEVEL == ACE_MUET)
    return;

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
  if (ACE_MUET_LEVEL == ACE_MUET)
    return;

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

//**********************************************
//**********************************************
//**********************************************
//**********************************************
//ADAPTIVE TIME STEPPING



/*
 *mStepCmp%2==0 : nothing.
 *mStepCmp=1 : initialisation.
 *mStepCmp%2==1 : backtrack, and compute new step, remenber current state.
 *
 */
void  linearSystem::computeAndAcceptStep(){
  if (!ACE_WITH_ADAPTATIVE_TIME_STEPPING){
    printStep(*mSimuStream,mzsti);
    return;
  }
  mAdaptiveStepEvaluation=true;
  ACE_DOUBLE timeSav = mTcurrent;
  if (mStepCmp%2==0 ){
    *mzst_inter=*mzsti;
  }else{
    //initialisation.
    if (mStepCmp == 1){
      mNormX0=fmax(1,mxti->norm2());
      mNormZ0=fmax(1,mzsti->norm2());
      //sav current state to be able to come back in tow steps
      *mxtiprev=  *mxti ;
      *mzstiprev=*mzsti ;
      *mzst1 = *mzsti;
      *mxt1 = *mxti;
      printStep(*mSimuStream,mzsti);

    }else{
      mNbBacktrack++;
      //backtrack excepted for the first step
      int cmpSav = mStepCmp;
      ACE_DOUBLE timeSav = mTcurrent;
      ACE_DOUBLE H = mH;
      ACE_DOUBLE newH = mH;
      ACE_DOUBLE errorX, errorZs=0;
      ACE_DOUBLE coef;
      ACE_DOUBLE buf1,buf2;
      //save current state for error evaluation.
      *mxticurrent = *mxti;
      *mzsticurrent = *mzsti;
      //come back in the past.
      *mxti = *mxtiprev;
      *mzsti = *mzstiprev;
      //      cout<<"mzstiprev \n"<<*mzstiprev;
      mTcurrent-=2*H;
      //set step for evaluation
      setStep(2*H);
      //run
      preparMLCP();
      bool l_accept = false;
      if (step()){
	//compute error
	//      cout<<"mzsti \n"<<*mzsti;
	//      cout<<"mzsticurrent \n"<<*mzsticurrent;
	ACE_DOUBLE NormCurrX = mxti->norm2();
	*mxbuf=*mxti-*mxticurrent;
	*mzsbuf=*mzsti-*mzsticurrent;
	buf2=mzsbuf->normInf();
	errorZs=fmax(ACE_INF,buf2);
	buf2=mxbuf->normInf();
	errorX=fmax(ACE_INF,buf2);
	mSerrorX+=errorX;
	mSerrorZs+=errorZs;
	//compute new step
      
	coef = //fmin(
	  mNormZ0*ACE_RTOL_LOCAL/errorZs;//,
	//		    fmax(mNormX0,NormCurrX)*ACE_RTOL_LOCAL/errorX
	//	    );
	if (coef >=1){
	  coef = //fmin(
	    fmin(mAlphaMax,fmax(mAlphaMin,mAlpha*mNormZ0*ACE_RTOL_LOCAL/errorZs));//,
	    //     fmin(mAlphaMax,fmax(mAlphaMin,mAlpha*mNormX0*ACE_RTOL_LOCAL/errorX))
	    //      );
	  l_accept =true;
	}
#ifdef PRE_COMPUTE_ADAPTIVE
	ACE_CHECK_IERROR(ACE_CUR_STEP>=1,"linearSystem::computeAndAcceptStep : ACE_CUR_STEP>=1");
	int curCoef = (1<<(ACE_CUR_STEP-1));
	ACE_DOUBLE newcoef = coef * curCoef;
	newH=newcoef*H;
      }else{
	newH=mHori;
      }
	
#else
      //      cout<<"tol "<<ACE_MAX_LOCAL_ERROR<<endl;
      //      cout<<"error "<<error<<" ie: coef "<<coef<<endl;
      newH = coef * H;
    }else{
	newH=mHori;
    }
#endif
    
    if (H<=mHori)
      l_accept=true;

      if (newH < mHori){
	newH = mHori;
	mNbToSmall++;
      }else if (newH > mMaxHori){
	newH = mMaxHori;
	mNbToBig++;
      }
      //      if (coef<0.9 || coef >1.5){
      if (timeSav + newH > mTstop)
	setStep(mTstop - timeSav);
      else
	setStep(newH);
      //}

    

      if (l_accept){
	//accept the two previous step and go head
	//restore current values.
	*mxti=*mxticurrent;
	*mzsti=*mzsticurrent;
	mStepCmp=cmpSav;
	mTcurrent=timeSav-H;
	printStep(*mSimuStream,mzst_inter);
	mTcurrent=timeSav;
	printStep(*mSimuStream,mzsti);
	
	//sav current state to be able to come back in tow steps
	*mxtiprev=  *mxti ;
	*mzstiprev=*mzsti ;
      }else{
	//Do not eccept the tow previous step, and try with the new adptive step
	//come back in the past.
	*mxti = *mxtiprev;
	*mzsti = *mzstiprev;
	mStepCmp=cmpSav;
	mTcurrent=timeSav-2*H;
      }
      mAdaptiveStepEvaluation=false;
  }
}
}

void  linearSystem::computeAndAcceptStep_2(){
  if (!ACE_WITH_ADAPTATIVE_TIME_STEPPING){
    printStep(*mSimuStream,mzsti);
    return;
  }
//    if (mStepCmp == 5){
//      printStep(*mSimuStream,mzsti);
//      setStep(5*mHori);
//    }else{
//      printStep(*mSimuStream,mzsti);
//    }
//    return;
  mAdaptiveStepEvaluation=true;
  ACE_DOUBLE timeSav = mTcurrent;
  if (mStepCmp%2==0 ){
    *mzst_inter=*mzsti;
  }else{
    //initialisation.
    if (mStepCmp == 1){
      //sav current state to be able to come back in tow steps
      *mxtiprev=  *mxti ;
      *mzstiprev=*mzsti ;
      *mzst1 = *mzsti;
      *mxt1 = *mxti;
      printStep(*mSimuStream,mzsti);
    }else{
      mNbBacktrack++;
      //backtrack excepted for the first step
      int cmpSav = mStepCmp;
      ACE_DOUBLE H = mH;
      ACE_DOUBLE newH = mH;
      ACE_DOUBLE errorX, errorZs=0;
      ACE_DOUBLE coef;
      ACE_DOUBLE buf1,buf2;
      //save current state for error evaluation.
      *mxticurrent = *mxti;
      *mzsticurrent = *mzsti;
      //come back in the past.
      *mxti = *mxtiprev;
      *mzsti = *mzstiprev;
      //      cout<<"mzstiprev \n"<<*mzstiprev;
      mTcurrent-=2*H;
      //set step for evaluation
      setStep(2*H);
      //run
      preparMLCP();
      ACE_DOUBLE bidon = mzsti->getValue(0);
      bool l_accept = false;
      if (step()){
	//compute error
	//      cout<<"mzsti \n"<<*mzsti;
	//      cout<<"mzsticurrent \n"<<*mzsticurrent;
	*mxbuf=*mxti-*mxticurrent;
	*mzsbuf=*mzsti-*mzsticurrent;
	ACE_DOUBLE err=0;
	ACE_DOUBLE aux=0;
	ACE_DOUBLE sci=0;
	int i;
	for (i=0; i < mDimzs-1; i++){
	  sci = ACE_ATOL_LOCAL + ACE_RTOL_LOCAL*fmax(mzst1->getValue(i),mzsti->getValue(i));
	  aux=(mzsbuf->getValue(i)/sci);
	  err+=aux*aux;
	}
	for (i=0; i < mDimx; i++){
	  sci = ACE_ATOL_LOCAL + ACE_RTOL_LOCAL*fmax(mxt1->getValue(i),mxti->getValue(i));
	  aux=(mxbuf->getValue(i)/sci);
	  err+=aux*aux;
	}
	err=err/(mDimzs-1 + mDimx);
// 	err=err/(mDimzs-1);
	err=sqrt(err);

	coef = 1/err;
	if (coef >= 1){
	  l_accept = true;
	  coef = fmin(mAlphaMax,fmax(mAlphaMin,mAlpha/err));
	}
	//	coef=2;
#ifdef PRE_COMPUTE_ADAPTIVE
	ACE_CHECK_IERROR(ACE_CUR_STEP>=1,"linearSystem::computeAndAcceptStep : ACE_CUR_STEP>=1");
	int curCoef = (1<<(ACE_CUR_STEP-1));
	ACE_DOUBLE newcoef = coef * curCoef;
	newH=newcoef*H;
      }else{
	newH=mHori;
      }
	
#else
      //      cout<<"tol "<<ACE_MAX_LOCAL_ERROR<<endl;
      //      cout<<"error "<<error<<" ie: coef "<<coef<<endl;
      newH = coef * H;
    }else{
	newH=mHori;
    }
#endif
//     newH=mHori;
//     if (mNbBacktrack ==1)
//       newH=2*mHori;
//     if (mNbBacktrack ==2)
//       newH=4*mHori;
//     if (mNbBacktrack ==3)
//       newH=10*mHori;
//     if (mNbBacktrack ==4)
//       newH=3*mHori;
    

      if (H<=mHori)
	l_accept=true;
    
      if (newH < mHori){
	newH = mHori;
	mNbToSmall++;
      }else if (newH > mMaxHori){
	newH = mMaxHori;
	mNbToBig++;
      }

      if (timeSav + newH > mTstop)
	setStep(mTstop - timeSav);
      else
	setStep(newH);

    

      if (l_accept){
	//accept the two previous step and go head
	//restore current values.
	*mxti=*mxticurrent;
	*mzsti=*mzsticurrent;
	mStepCmp=cmpSav;
	mTcurrent=timeSav-H;
	printStep(*mSimuStream,mzst_inter);
	mTcurrent=timeSav;
	printStep(*mSimuStream,mzsti);
	
	//sav current state to be able to come back in tow steps
	*mxtiprev=  *mxti ;
	*mzstiprev=*mzsti ;
      }else{
	//Do not eccept the tow previous step, and try with the new adptive step
	//come back in the past.
	*mxti = *mxtiprev;
	*mzsti = *mzstiprev;
	mStepCmp=cmpSav;
	mTcurrent=timeSav-2*H;
      }
      mAdaptiveStepEvaluation=false;
  }
}
}
void linearSystem::setStep(ACE_DOUBLE newH){
#ifdef PRE_COMPUTE_ADAPTIVE
  int aux = lround(floor(newH/mHori));
  aux=aux>>1;
  ACE_CUR_STEP=0;
  while (aux && ACE_CUR_STEP < ACE_NB_ADAPT_STEP){
    aux=aux>>1;
    ACE_CUR_STEP++;
  }
    
  ACE_CHECK_IERROR(ACE_CUR_STEP<=ACE_NB_ADAPT_STEP,"linearSystem::setStep ACE_CUR_STEP<=ACE_NB_ADAPT_STEP ");
  mH = mHori*(1<<ACE_CUR_STEP);
//   cout<<"******************************** ";
//   cout<<"new ACE_CUR_STEP : "<<ACE_CUR_STEP<<" mH: "<<mH<<endl;
  fillMLCP();
#else
  mH=newH;
  buildMLCP();
  fillMLCP();
#endif
  //  ACE_ERROR("linearSystem::setStep not implemeted");
}
