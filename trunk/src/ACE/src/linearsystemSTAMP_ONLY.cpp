/************************************************************************
  			linearsystemDAE.cpp 
**************************************************************************/


#include "linearsystemSTAMP_ONLY.h"
#include "SimpleMatrix.h"

using namespace std;

linearSystemSTAMP_ONLY::linearSystemSTAMP_ONLY():linearSystemMNA_V(){
  mA_A1=0;
  mA1sti=0;

}

linearSystemSTAMP_ONLY::~linearSystemSTAMP_ONLY(){
  if(mA_A1)
    delete mA_A1;
  if (mA1sti)
    delete mA1sti;
}




void linearSystemSTAMP_ONLY::buildMLCP(){
  if (!mMLCP){
    mMLCP = new mlcp(mDimLambda, mNbDynEquations + mDimzns,ACE_SOLVER_TYPE);
    mMLCP->mM22->eye();
    (*mMLCP->mM22)=-1*(*mMLCP->mM22);
    mMLCP->mM21->zero();
    mA_A1 = new aceMatrix(mNbDynEquations,mDimx-mV0x,ACE_MAT_TYPE);
    mA1sti = new aceVector(mNbDynEquations,ACE_MAT_TYPE);
  }
  *mA_A1=*mA+(mH*(1-mTheta))*(*mA1x);

  for (int i=0;i<ACE_NB_ADAPT_STEP+1;i++){

#ifdef PRE_COMPUTE_ADAPTIVE
    if (ACE_WITH_ADAPTATIVE_TIME_STEPPING)
      mH = mHori*(1<<i);
#endif
    try{
      //ACE_CHECK_IERROR(mDimx >0 && mDimLambda >0,"linearSystem::initSimu case no x or no lambda not yet implemented");
      if (mDimx){
	*(mW[i]) = -1*(*mA) + mH*mTheta*(*mA1x);
      }
    }
    catch(SiconosException e)
      {
	std::cout << e.report() << endl;
	ACE_INTERNAL_ERROR("linearSystemSTAMP_ONLY::buildMLCP");
      }
    catch(...)
      {
	std::cout << "Exception caught." << endl;
	ACE_INTERNAL_ERROR("linearSystem::buildMLCP");
      }
  }
}
void linearSystemSTAMP_ONLY::fillMLCP(){
  *(mMLCP->mM11) = *mD1l;
  
  mMLCP->mM12->setBlock(0,0,*mD1x);
  if(mD1zs)
    mMLCP->mM12->setBlock(0,mDimx-mV0x,*mD1zs);
  mMLCP->mM12->setBlock(0,mDimx-mV0x+mDimzs-mV0zs,*mD1zns);

  mMLCP->mM21->setBlock(mNbDynEquations,0,*mC1l);

  mMLCP->mM22->setBlock(0,0,*(mW[ACE_CUR_STEP]));
  if (mA1zs)
    mMLCP->mM22->setBlock(0,mDimx-mV0x,(mH*mTheta)*(*mA1zs));
  mMLCP->mM22->setBlock(0,mDimx-mV0x+mDimzs-mV0zs,(mH*mTheta)*(*mA1zns));
  
  mMLCP->mM22->setBlock(mNbDynEquations,0,*mC1x);
  if (mC1zs)
    mMLCP->mM22->setBlock(mNbDynEquations,mDimx-mV0x,*mC1zs);
}
void linearSystemSTAMP_ONLY::preparMLCP(){
  if (mDimx && mDimLambda){//both
    //(*mxfree) = (mA + h(1-Theta)mA1x)(*mxti) + mH*((1-mTheta)*(prod(*mA1x,*mxti)+prod(*mA1zs,*mzsti)+prod(*mA1zns,*mznsti) + (*mA1sti))+mThetap*(*mA2s));
    //mxfree->display();
    if (mA1zs)
      ACEprod(*mA1zs,*mzsti,*mxfree,true);
    else
      mxfree->zero();
    ACEprod(*mA1zns,*mznsti,*mxfree,false);
    scal((1-mTheta),*mxfree,*mxfree);
    *mxfree+=mThetap*(*mA1s);
    *mxfree+=(1-mThetap)*(*mA1sti);
    scal(mH,*mxfree,*mxfree);
    ACEprod(*mA_A1,*mxti,*mxfree,false);
    
  } else if(mDimx){//only x
      ACE_INTERNAL_ERROR("linearSystemSTAMP_ONLY::preparMLCP: not yet implemented");
  }else if (mDimLambda){//only lambda
      ACE_INTERNAL_ERROR("linearSystemSTAMP_ONLY::preparMLCP: not yet implemented");
  }else{
      ACE_INTERNAL_ERROR("linearSystemSTAMP_ONLY::preparMLCP: not yet implemented");
  }
  
  mMLCP->mQ2->setBlock(0,*mxfree);
  mMLCP->mQ2->setBlock(mNbDynEquations,*mC1s);
  *( mMLCP->mQ1)=*mD1s;
}

bool linearSystemSTAMP_ONLY::step(){
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
  bool res = mMLCP->solve();
  //  mMLCP->printOutPut();
  ACE_times[ACE_TIMER_COMPUTE_VAR].start();
  if (res){
    if (mDimx){
      for (int i=0; i< mxti->size() ; i++){
	mxti->setValue(i,mMLCP->mZ2->getValue(i));
      }
//       cout<<"*************************mxti\n";
//       mxti->display();
      if (mzsti)
	for (int i=0; i< mzsti->size(); i++){
	  mzsti->setValue(i,mMLCP->mZ2->getValue(mDimx-mV0x+i));
	}
      for (int i=0; i< mznsti->size(); i++){
	mznsti->setValue(i,mMLCP->mZ2->getValue(mDimx-mV0x+mDimzs-mV0zs+i));
      }
//       cout<<"*************************mznsti\n";
//       mznsti->display();
//       cout<<"*************************mZ1\n";
//       mMLCP->mZ1->display();
     }

  }else{
    ACE_GET_LOG_STREAM()<<"linearSystemSTAMP_ONLY::step number,"<< mStepCmp<<" solver failled!!!"<<endl;
  }
  ACE_times[ACE_TIMER_COMPUTE_VAR].stop();
  ACE_times[ACE_TIMER_LS_STEP].stop();
  return res;
}

void linearSystemSTAMP_ONLY::ExtractAndCompute2Sources(){
  if (mDimx)
    *mA1sti=*mA1s;
  extractSources();
//   cout<<"linearSystemSTAMP_ONLY::ExtractAndCompute2Sources"<<endl;
//   mA1s->display();
  
}

void linearSystemSTAMP_ONLY::extractDynamicSystemSource(){
  if (mDimx){
    extractDynBockInVect(ms);
  }

}
