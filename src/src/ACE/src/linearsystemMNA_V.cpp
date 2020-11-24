/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2019 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/

/*! \file linearsystemMNA_V.cpp

*/
/************************************************************************
  			linearsystemDAE.cpp 
**************************************************************************/


#include "linearsystemMNA_V.h"
#include "SimpleMatrix.h"

using namespace std;
//|--------x'---------|----------x--------|-------Zs------|------Zns------|

linearSystemMNA_V::linearSystemMNA_V(){
  mAdaptiveStepEvaluation = false;
  mV0x=1;
  mV0zs=0;

}

linearSystemMNA_V::~linearSystemMNA_V(){
}
void linearSystemMNA_V::allocMemory(){
  ACE_CHECK_IERROR(!mReAlloc,"linearSystem::allocMemory again");
  mReAlloc = true;
  mNbUnknowns = 2*mx.size()+mZs.size()+mZns.size();
  mRS = mNbUnknowns;
  int index=-1;
  int i;
  int n=mKCL.size();
  for (i=0;i<n;i++){
    mKCL[i]->allocMemory(mNbUnknowns+1);
    mKCL[i]->mIsDyn=true;
    mKCL[i]->mLine=index;
    index++;
  }

  n=mIND.size();
  for (i=0;i<n;i++){
    mIND[i]->allocMemory(mNbUnknowns+1);
    mIND[i]->mIsDyn=true;
    mIND[i]->mLine=index;
    index++;
  }
  n=mCAP.size();
  for (i=0;i<n;i++){
    mCAP[i]->allocMemory(mNbUnknowns+1);
    mCAP[i]->mIsDyn=true;
    mCAP[i]->mLine=index;
    index++;
  }

  
  n=mVD.size();
  for (i=0;i<n;i++){
    mVD[i]->allocMemory(mNbUnknowns+1);
    mVD[i]->mIsDyn=true;
    mVD[i]->mLine=index;
    index++;
  }
  n=mTEN.size();
  for (i=0;i<n;i++){
    mTEN[i]->allocMemory(mNbUnknowns+1);
    mTEN[i]->mIsDyn=true;
    mTEN[i]->mLine=index;
    index++;
  }
  mNbDynEquations = mNbEquations;
  allocA1Matrix();
  allocC1Matrix();
  allocD1Matrix();



}

void  linearSystemMNA_V::addVUnknowns(){
  ACE_CHECK_IERROR(mNbNodes,"linearSystem::initAddVUnknowns no nodes");
 for(int i=0;i<mNbNodes;i++){
   mx.push_back(new unknown(ACE_TYPE_V,i));
 }
}
int linearSystemMNA_V::getIndexUnknown (int type,int node) {
  if (type == ACE_TYPE_V){
      return node + mx.size();
  }else{
    ACE_INTERNAL_WARNING("linearSystemMNA_V::getIndexUnknown not implemented");
    return -1;
  }
  
}

int linearSystemMNA_V::getDynIndexUnknown (int type,int node) {
  if (type == ACE_TYPE_V){
    return node ;
  }else{
    ACE_INTERNAL_WARNING("linearSystemMNA_V::getIndexUnknown not implemented");
    return -1;
  }
  
}



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////DISCRETISATION
void linearSystemMNA_V::buildMLCP(){

  if (!mMLCP)
    mMLCP = new mlcp(mDimLambda,mNbDynEquations,ACE_SOLVER_TYPE);
  
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
    if (mDimzs)
      scal(mH*mTheta,*mA2zs,*mHThetaA2zs[i]);
  }
  if (mDimLambda){
    *(mMLCP->mM11) = *mD2l;
    mMLCP->mM12->setBlock(0,0,*mD2x);
    if (mDimzs)
      mMLCP->mM12->setBlock(0,mDimx-mV0x,*mD2zs);
  }
  
}
void linearSystemMNA_V::fillMLCP(){
   if(mDimLambda)
     mMLCP->mM21->setBlock(0,0,*mhR[ACE_CUR_STEP]);
   mMLCP->mM22->setBlock(0,0,*mW[ACE_CUR_STEP]);
   if (mDimzs)
     mMLCP->mM22->setBlock(0,mDimx-mV0x,*mHThetaA2zs[ACE_CUR_STEP]);
   mMLCP->update();
}
void linearSystemMNA_V::preparMLCP(){
  if (mDimx ){//both
    //(*mxfree) = mA*(*mxti) + mH*((1-mTheta)*(prod(*mA2x,*mxti)+prod(*mA2zs,*mzsti) + (*mA2sti))+mThetap*(*mA2s));
    //mxfree->display();
    if (mDimzs){
      ACEprod(*mA2zs,*mzsti,*mxfree,true);
      ACEprod(*mA2x,*mxti,*mxfree,false);
      *mxfree+=*mA2sti;
    }else
      ACEprod(*mA2x,*mxti,*mxfree,true);
    scal((1-mTheta),*mxfree,*mxfree);
    if (mDimzs)
      *mxfree+=mThetap*(*mA2s);
    scal(mH,*mxfree,*mxfree);
    ACEprod(*mA,*mxti,*mxfree,false);
    //    cout<<"mxfree"<<*mxfree<<endl;
 
  } /*else if(mDimx){
      (*mxfree) =  mH*(1-mTheta)*(prod(*mA2x,*mxti)+prod(*mA2zs,*mzsti)+(*mA2sti))+mH*(*mA2s);
    ACEprod(*mA,*mxti,*mxfree,false);
    }*/else if (mDimLambda){//only lambda
    ;
  }else{
    ;
  }
//   cout<<"mQfree\n";
//   (mQfree)->display();
  if(ACE_MAT_TYPE == SPARSE){
    ACE_ERROR("linearSystemMNA_V::preparStep not implemented\n");
//     (mMLCP->mQ1)->set(*(mPfree));
//     scal(-1,*(mMLCP->mQ1),*(mMLCP->mQ1));
//     (mMLCP->mQ2)->set(*(mQfree));
//     scal(-1,*(mMLCP->mQ2),*(mMLCP->mQ2));
  }else{
    if (mDimLambda)
      *(mMLCP->mQ1)= *mD2s;
    //      scal(-1,*mD2s,*(mMLCP->mQ1));
    mMLCP->mQ2->setBlock(0,*mxfree);
    //    scal(-1,*(mMLCP->mQ2),*(mMLCP->mQ2));
  }

}

bool linearSystemMNA_V::step(){
  //  cout<<"*******************begin step "<<mTcurrent<<" to ";
  //     cout<<"before\n"<< *(mMLCP->mM22)<<endl;

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
  if(!mDimLambda)
    fillMLCP();
  ACE_times[ACE_TIMER_COMPUTE_VAR].start();
  if (res){
    //*mzsti=*(mMLCP->mZ2);
    if (ACE_MAT_TYPE == SPARSE){
      ACE_ERROR("linearSystemMNA_V::preparStep not implemented\n");
    }else{
      for (int i=0; i< mDimzs-1; i++){
	mzsti->setValue(i,mMLCP->mZ2->getValue(mDimx-mV0x+i));
      }
    }
    if (mDimx){
      for (int i=0; i< mDimx-mV0x; i++){
	mxti->setValue(i,mMLCP->mZ2->getValue(i));
      }
    }
    computeZnstiFromX_Zs();
  }else{
    ACE_WARNING("linearSystemMNA_V::step, simulation stoped because mlcp solver failed!\n");
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
  //Y=D2x*x + D2zs*Zs + D2l*lambda + D2s
void linearSystemMNA_V::printSystem2(ostream& os){
  if (ACE_MUET_LEVEL == ACE_MUET)
    return;
    
  os<<"R=A1zns*C1l\nAx'=A2x*x + A2zs*Zs + R*lambda+A2s\nY=D2x*x + D2zs*Zs + D2l*lambda + D2s\n";
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
void linearSystemMNA_V::printA1(ostream& os ){
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
void linearSystemMNA_V::extractDynamicSystemSource(){
  if (mDimx)
    extractDynBockInVect(ms);// NB : mA1s=ms
  //  cout<<"source\n"<<*ms<<endl;
  //  cout<<"ms\n"<<*ms<<endl;
  

}


void linearSystemMNA_V::printStep(ostream& os,aceVector * pVx,aceVector *pVzs){
//   int i;
   ACE_DOUBLE curDate = getCurrentTime();

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
    aux = pVx->getValue(pPrint->node1-1);
    if (pPrint->node2 >0)
      aux -=  pVx->getValue(pPrint->node2-1);
    // mSommeError += fabs(v-aux);
    
      os << aux ;
      //      os << aux << "\t"<<v;
  }
  os<<endl;
}
