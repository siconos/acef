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
void linearSystemDAE::preparStep(){
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
bool linearSystemDAE::step(){
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
  return res;
}

