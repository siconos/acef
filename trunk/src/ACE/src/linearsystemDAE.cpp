/************************************************************************
  			linearsystemDAE.cpp 
**************************************************************************/


#include "linearsystemDAE.h"
#include "SimpleMatrix.h"

using namespace std;
//|--------x'---------|----------x--------|-------Zs------|------Zns------|

linearSystemDAE::linearSystemDAE(){
  mAdaptiveStepEvaluation = false;
}

linearSystemDAE::~linearSystemDAE(){
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////DISCRETISATION
void linearSystemDAE::buildMLCP(){

  if (!mMLCP)
    mMLCP = new mlcp(mDimLambda,mNbNonDynEquations+mDimx,ACE_SOLVER_TYPE);
  
  for (int i=0;i<ACE_NB_ADAPT_STEP+1;i++){

#ifdef PRE_COMPUTE_ADAPTIVE
    if (ACE_WITH_ADAPTATIVE_TIME_STEPPING)
	mH=mHori*(1<<i);
#endif
    if(mDimLambda){
      scal(mH,*mR,*mhR[i]);
    }
    //W=A-h*theta A2x
    mW[i]->zero();
    *mW[i]-=*mA;
    *mW[i] += mH*mTheta*(*mA2x);
    scal(mH*mTheta,*mA2zs,*mHThetaA2zs[i]);
  }
  *(mMLCP->mM11) = *mD2l;
  mMLCP->mM12->setBlock(0,0,*mD2x);
  mMLCP->mM12->setBlock(0,mDimx,*mD2zs);
  
  mMLCP->mM21->setBlock(mDimx,0,*mB2l);
  
  mMLCP->mM22->setBlock(mDimx,0,*mB2x);
  mMLCP->mM22->setBlock(mDimx,mDimx,*mB2zs );
}
void linearSystemDAE::fillMLCP(){
      mMLCP->mM21->setBlock(0,0,*mhR[ACE_CUR_STEP]);
      mMLCP->mM22->setBlock(0,0,*mW[ACE_CUR_STEP]);
      mMLCP->mM22->setBlock(0,mDimx,*mHThetaA2zs[ACE_CUR_STEP]);
      mMLCP->update();
}
void linearSystemDAE::preparMLCP(){
    if (mDimx && mDimLambda){//both
    //(*mxfree) = mA*(*mxti) + mH*((1-mTheta)*(prod(*mA2x,*mxti)+prod(*mA2zs,*mzsti) + (*mA2sti))+mThetap*(*mA2s));
    //mxfree->display();
    ACEprod(*mA2zs,*mzsti,*mxfree,true);
    ACEprod(*mA2x,*mxti,*mxfree,false);
    *mxfree+=*mA2sti;
    scal((1-mTheta),*mxfree,*mxfree);
    *mxfree+=mThetap*(*mA2s);
    scal(mH,*mxfree,*mxfree);
    ACEprod(*mA,*mxti,*mxfree,false);
 
  } else if(mDimx){//only x
    (*mxfree) =  mH*(1-mTheta)*prod(*mA2x,*mxti)+mH*(1-mTheta)*prod(*mA2zs,*mzsti)+mH*(*mA2s);
    ACEprod(*mA,*mxti,*mxfree,false);
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
      *(mMLCP->mQ1)= *mD2s;
    //      scal(-1,*mD2s,*(mMLCP->mQ1));
    mMLCP->mQ2->setBlock(0,*mxfree);
    mMLCP->mQ2->setBlock(mDimx,*mB2s);
    //    scal(-1,*(mMLCP->mQ2),*(mMLCP->mQ2));
  }

}

bool linearSystemDAE::step(){
  //  cout<<"*******************begin step "<<mTcurrent<<" to ";
  ACE_times[ACE_TIMER_LS_STEP].start();
  mStepCmp++;
  mAllStepCmp++;
  if (mTcurrent+mH > mTstop){
    ACE_times[ACE_TIMER_LS_STEP].stop();
    return false;
  }
  //  cout<<mTcurrent<<endl;
  //  cout<<"before mzsti"<<*mzsti;
  if (mStepCmp%mLogFrequency==0 ){
    printf("-->%d : %e to %e\n",mPourMille, mTcurrent,mTcurrent + mH);
    mPourMille ++;
  }
  mTcurrent += mH;
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
  }else{
    ACE_WARNING("linearSystemDAE::step, simulation stoped because mlcp solver failed!\n");
  }
  ACE_times[ACE_TIMER_COMPUTE_VAR].stop();
  ACE_times[ACE_TIMER_LS_STEP].stop();
  //  cout<<"after mzsti "<<*mzsti;
  //  cout<<"*******************end step : "<<res<<"\n";
  return res;
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////PRINT FONCTIONS
 //
  //R=A1zns*C1l
  //Ax'=A2x*x + A2zs*Zs + R*lambda+A2s
  //0=B2x*x + B2zs*Zs + B2l*lambda + B2s
  //Y=D2x*x + D2zs*Zs + D2l*lambda + D2s
void linearSystemDAE::printSystem2(ostream& os){
  if (ACE_MUET_LEVEL == ACE_MUET)
    return;
    
  os<<"R=A1zns*C1l\nAx'=A2x*x + A2zs*Zs + R*lambda+A2s\n0=B2x*x + B2zs*Zs + B2l*lambda + B2s\nY=D2x*x + D2zs*Zs + D2l*lambda + D2s\n";
  os<<"A:\n";
  if (mA)
    os<<(*mA);

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
void linearSystemDAE::printA1(ostream& os ){
  if (ACE_MUET_LEVEL == ACE_MUET)
    return;

  os <<"system Ax'=A1x+A1zs+A1zns+s\n";
  os <<"A:\n";
  if (mA)
    os <<(*mA);
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
void linearSystemDAE::extractDynamicSystemSource(){
  if (mDimx)
    extractDynBockInVect(ms);// NB : mA1s=ms
  

}
