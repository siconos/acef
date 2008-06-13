/************************************************************************
  			linearsystemDAE.cpp 
**************************************************************************/


#include "linearsystemDAE.h"
#include "SimpleMatrix.h"

using namespace std;
//|--------x'---------|----------x--------|-------Zs------|------Zns------|

linearSystemDAE::linearSystemDAE(){
}

linearSystemDAE::~linearSystemDAE(){
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////DISCRETISATION
void linearSystemDAE::buildMLCP(){
    mMLCP = new mlcp(mDimLambda,mNbNonDynEquations+mDimx,ACE_SOLVER_TYPE);
    scal(mH,*mR,*mhR);
    //W=I-h*theta A2x
    mW->zero();
    for(int i = 0; i < mDimx;i++)	mW->setValue(i,i,-1.0);
    *mW += mH*mTheta*(*mA2x);
    
    *(mMLCP->mM11) = *mD2l;
    mMLCP->mM12->setBlock(0,0,*mD2x);
    mMLCP->mM12->setBlock(0,mDimx,*mD2zs);
    
    mMLCP->mM21->setBlock(0,0,*mhR);
    mMLCP->mM21->setBlock(mDimx,0,*mB2l);

    mMLCP->mM22->setBlock(0,0,*mW);
    mMLCP->mM22->setBlock(0,mDimx,*mHThetaA2zs);
    mMLCP->mM22->setBlock(mDimx,0,*mB2x);
    mMLCP->mM22->setBlock(mDimx,mDimx,*mB2zs );

}
void linearSystemDAE::preparMLCP(){
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
 
  } else if(mDimx){//only x
    (*mxfree) = (*mxti) + mH*(1-mTheta)*prod(*mA2x,*mxti)+mH*(1-mTheta)*prod(*mA2zs,*mzsti)+mH*(*mA2s);
  }else if (mDimLambda){//only lambda
    ;
  }else{
    ;
  }
//   cout<<"mQfree\n";
//   (mQfree)->display();
  if(ACE_MAT_TYPE == SPARSE){
    ACE_ERROR("linearSystemDAE::preparStep not implemented\n");
//     (mMLCP->mQ1)->set(*(mPfree));
//     scal(-1,*(mMLCP->mQ1),*(mMLCP->mQ1));
//     (mMLCP->mQ2)->set(*(mQfree));
//     scal(-1,*(mMLCP->mQ2),*(mMLCP->mQ2));
  }else{
    if (mDimLambda)
      scal(-1,*mD2s,*(mMLCP->mQ1));
    mMLCP->mQ2->setBlock(0,*mxfree);
    mMLCP->mQ2->setBlock(mDimx,*mB2s);
    scal(-1,*(mMLCP->mQ2),*(mMLCP->mQ2));
  }

}
void linearSystemDAE::preparStep(){
  computeBestStep();
  preparMLCP();
}

bool linearSystemDAE::step(){
  cout<<"*******************begin step "<<mTcurrent<<" to ";
  ACE_times[ACE_TIMER_LS_STEP].start();
  mStepCmp++;
  if (mTcurrent > mTstop){
    ACE_times[ACE_TIMER_LS_STEP].stop();
    return false;
  }
  mTcurrent += mH;
  cout<<mTcurrent<<endl;
  cout<<"before mzsti"<<*mzsti;
  if (mStepCmp%mLogFrequency==0){
    printf("-->%d : %lf\n",mPourMille, mTcurrent);
    mPourMille ++;
  }
  bool res = mMLCP->solve();
  ACE_times[ACE_TIMER_COMPUTE_VAR].start();
  if (res){
    //*mzsti=*(mMLCP->mZ2);
    if (ACE_MAT_TYPE == SPARSE){
      ACE_ERROR("linearSystemDAE::preparStep not implemented\n");
      mzsti->set(*(mMLCP->mZ2));
    }else{
      for (int i=0; i< mDimzs-1; i++){
	mzsti->setValue(i,mMLCP->mZ2->getValue(mDimx+i));
      }
    }
    if (mDimx){
      for (int i=0; i< mDimx; i++){
	mxti->setValue(i,mMLCP->mZ2->getValue(i));
      }
    }
    computeZnstiFromX_Zs();
  }
  ACE_times[ACE_TIMER_COMPUTE_VAR].stop();
  ACE_times[ACE_TIMER_LS_STEP].stop();
  cout<<"after mzsti "<<*mzsti;
  cout<<"*******************end step : "<<res<<"\n";
  return res;
}

void linearSystemDAE::setStep(ACE_DOUBLE newH){
  
  //update mhR mHThetaA2zs and mW
  scal(newH,*mR,*mhR);
  scal(newH*mTheta,*mA2zs,*mHThetaA2zs);

  //W=I-h*theta A2x
  mW->zero();
  for(int i = 0; i < mDimx;i++)	mW->setValue(i,i,-1.0);
  *mW += mH*mTheta*(*mA2x);

  
  //update the matrix M of the MLCP
  
  mMLCP->mM21->setBlock(0,0,*mhR);
  mMLCP->mM22->setBlock(0,0,*mW);
  mMLCP->mM22->setBlock(0,mDimx,*mHThetaA2zs);

  mH=newH;
}
/*
 *mStepCmp%2==0 : nothing.
 *mStepCmp=1 : initialisation.
 *mStepCmp%2==1 : backtrack, and compute new step, remenber current state.
 *
 */

void  linearSystemDAE::computeBestStep(){
  if (!ACE_WITH_ADAPTATIVE_TIME_STEPPING)
    return;
  
  if (mStepCmp%2==1 ){
    //initialisation.
    if (mStepCmp == 1){
      mNormX0=fmax(1,mxti->norm2());
      mNormZ0=fmax(1,mzsti->norm2());
    }else{
      //backtrack excepted for the first step
      int cmpSav = mStepCmp;
      ACE_DOUBLE timeSav = mTcurrent;
      ACE_DOUBLE H = mH;
      ACE_DOUBLE newH = mH;
      ACE_DOUBLE error=0;
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
      step();
      //compute error
      //      cout<<"mzsti \n"<<*mzsti;
      //      cout<<"mzsticurrent \n"<<*mzsticurrent;
      *mxti=*mxti-*mxticurrent;
      *mzsti=*mzsti-*mzsticurrent;
      buf2=mzsti->normInf();
      error=fmax(ACE_INF,mzsti->normInf());
      //compute new step
      coef = fmin(mAlphaMax,fmax(mAlphaMin,mAlpha*mNormZ0*mLocalErrorTol/error));
      newH = coef * H;
      if (newH < mHori){
	newH = mHori;
	ACE_WARNING("linearSystemDAE::computeBestStep : try to set a to small time stepping");
      }
      if (newH > mMaxHori){
	newH = mMaxHori;
	ACE_WARNING("linearSystemDAE::computeBestStep : try to set a to large time stepping");
      }
      if (timeSav + newH > mTstop)
	setStep(mTstop - timeSav);
      else
	setStep(newH);
      printf("new H : %lf\n",newH);

      //restore current values.
      *mxti=*mxticurrent;
      *mzsti=*mzsticurrent;
      mStepCmp=cmpSav;
      mTcurrent=timeSav;
    }

    //sav current state to be able to come back in tow steps
    *mxtiprev=  *mxti ;
    *mzstiprev=*mzsti ;
 
  }
}
