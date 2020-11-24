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

/*! \file componentmos_nl2.cpp

*/
/************************************************************************
  			componentmos_nl2.cpp 

Convention:

A Non linear model without non-smooth laws.
\
**************************************************************************/

#include "componentmos_nl2.h"
#include "algo.h"



componentMOS_NL2::componentMOS_NL2(dataMOS1 *d)
:component_NONE_LINEAR_NS(){

  ACE_CHECK_IERROR(d,"componentMOS_NL2::componentMOS_NL2 : Mos data null");
  mType = ACE_TYPE_MOS_NL2;
  mData =(*d);
  mNodeD=mData.drain;
  mNodeG=mData.gate;
  mNodeS=mData.source;
  mB=mData.k/2.0;
  mK=mData.k;
  mMode = mData.mode;
  ACE_CHECK_IWARNING(mMode == 1 || mMode == -1,"componentMOS_NL2 mode value not 1 or -1.");
  mVt=mMode*mData.vt;
  cout<<"mMode=\n"<<mMode;
  cout<<"mVt=\n"<<mVt;
  cout<<"mK=\n"<<mK;
  cout<<"mB=\n"<<mB;
  

  mNodeNeg=mNodeD;
  mNodePos=mNodeS;

  mName = mData.name;
  ACE_CHECK_ERROR(mNodeD>=0 && mNodeG>=0 && mNodeS >=0,"componentMOS_NL2::componentMOS_NL2");

  //process parameters
  double VI = ACE_MOS_POWER_SUPPLY;	// Power supply
  double Vt0 = mData.vt;
  double Kval = mData.k;
  double HalfK = Kval/2.0;
  double ractolpwl;
  double widthhyp;

  mDimlambda=0;
  mDimZns=0;
  mIndiceStartZns=-1;
  mIndiceStartLambda=-1;

  
}

void componentMOS_NL2::addUnknowns(){
  mI=algo::spls->addinZs(ACE_TYPE_I,this);
}
void  componentMOS_NL2::addEquations (){
  mILaw = algo::spls->addNonLinearEquation();
}

void componentMOS_NL2::stamp(){
  if (ACE_FORMULATION ==  ACE_FORMULATION_MNA_V || ACE_FORMULATION ==  ACE_FORMULATION_STAMP_ONLY){
    stampMNA_V();
    return;
  }

  int ind=0;
  int i=mI->mIndex;
  //stamp equations.
  algo::spls->KCL(mNodeNeg)->mCoefs[i]-=1;
  algo::spls->KCL(mNodePos)->mCoefs[i]+=1;


  
}
void componentMOS_NL2::stampMNA_V(){

  ACE_ERROR("componentMOS_NL2::stampMNA_V not yet implemented.\n");
//   int ind=0;
//   int i=mI->mIndex;
//   //stamp equations.
//   algo::spls->KCL(mNodeNeg)->mCoefs[i]-=mMode*1;
//   algo::spls->KCL(mNodePos)->mCoefs[i]+=mMode*1;

//   //Zns = B*lamdba
//   for (ind=0;ind < mNbHyp;ind++){
//     algo::spls->mC1l->setValue(mIndiceStartZns,mIndiceStartLambda+ind,mCoefs[ind]);
//     algo::spls->mC1l->setValue(mIndiceStartZns,mIndiceStartLambda+mNbHyp+ind,-mCoefs[ind]);
//   }
  
//   //Y=C*zs+I*lambda+hyp

//   //C*zs
//   if (mNodeD){
//     for(ind=0;ind < mNbHyp;ind++)
//       algo::spls->mD1x->setValue(mIndiceStartLambda+mNbHyp+ind,mNodeD-1,mMode);
//   }
//   if (mNodeG){
//     for(ind=0;ind < mDimlambda;ind++)
//       algo::spls->mD1x->setValue(mIndiceStartLambda+ind,mNodeG-1,-mMode);
//   }
//   if (mNodeS){
//     for(ind=0;ind < mNbHyp;ind++)
//       algo::spls->mD1x->setValue(mIndiceStartLambda+ind,mNodeS-1,mMode);
//   }

//   //I*lambda
//   for(ind=0;ind < mDimlambda;ind++)
//     algo::spls->mD1l->setValue(mIndiceStartLambda+ind,mIndiceStartLambda+ind,1);
  
//   //hyp
//   for(ind=0;ind < mNbHyp;ind++){
//       algo::spls->mD1s->setValue(mIndiceStartLambda+ind,mHyp[ind]);
//       algo::spls->mD1s->setValue(mIndiceStartLambda+mNbHyp+ind,mHyp[ind]);
//     }

  
}

void componentMOS_NL2::computeNL(SiconosVector& SICONOS_X,SiconosVector& SICONOS_Lambda,SiconosVector& SICONOS_H){
  ACE_DOUBLE VGS=0;
  ACE_DOUBLE VGD=0;
  if (mNodeG){
    VGS+=SICONOS_Lambda.getValue(mNodeG-1);
    VGD+=SICONOS_Lambda.getValue(mNodeG-1);
    //    cout<<"VG "<<VGS<<endl;
  }
  if (mNodeS){
    VGS-=SICONOS_Lambda.getValue(mNodeS-1);
    //   cout<<"VS "<<SICONOS_Lambda.getValue(mNodeS-1)<<endl;
  }
  if (mNodeD){
    VGD-=SICONOS_Lambda.getValue(mNodeD-1);
    //cout<<"VD "<<SICONOS_Lambda.getValue(mNodeD-1)<<endl;
  }
  ACE_DOUBLE IDS = SICONOS_Lambda.getValue(mI->mIndex-2*algo::spls->mDimx-1);
//      ACE_DOUBLE aux=-IDS;
//    ACE_DOUBLE val;
//    cout<<"-IDS: "<<-IDS<<endl;
   
//    if (mMode*(VGS-mMode*mVt)>0)
//      aux+=mMode*mB*(VGS-mMode*mVt)*(VGS-mMode*mVt);
//    if (mMode*(VGD-mMode*mVt)>0)
//      aux-=mMode*mB*(VGD-mMode*mVt)*(VGD-mMode*mVt);
  
//    printf("L2=%e, L4=%e\n",L2,L4);
//    val = -IDS+mMode*mB*(L4*(VGS-mMode*mVt)*(VGS-mMode*mVt) - L2*(VGD-mMode*mVt)*(VGD-mMode*mVt));
//    printf("LCPval=%e, f=%e\n",val,aux);

  ACE_DOUBLE AUX = -IDS;
  if (mMode>0){
    if (VGS-mVt>0)
      AUX+=mB*(VGS-mVt)*(VGS-mVt);
    if (VGD-mVt>0)
      AUX-=mB*(VGD-mVt)*(VGD-mVt);
  }else{
    if (VGS+mVt<0)
      AUX-=mB*(VGS+mVt)*(VGS+mVt);
    if (VGD+mVt<0)
      AUX+=mB*(VGD+mVt)*(VGD+mVt);    
  }
  //  SICONOS_H.setValue(mILaw->mLine,-IDS+mMode*mB*(L4*(VGS-mMode*mVt)*(VGS-mMode*mVt) - L2*(VGD-mMode*mVt)*(VGD-mMode*mVt)));
  SICONOS_H.setValue(mILaw->mLine,AUX);
  

  
}
void componentMOS_NL2::computeJacNL(SiconosVector& SICONOS_X,SiconosVector& SICONOS_Lambda,SiconosMatrix& SICONOS_D){
  ACE_DOUBLE VGS=0;
  ACE_DOUBLE VGD=0;
  if (mNodeG){
    VGS+=SICONOS_Lambda.getValue(mNodeG-1);
    VGD+=SICONOS_Lambda.getValue(mNodeG-1);
  }
  if (mNodeS)
    VGS-=SICONOS_Lambda.getValue(mNodeS-1);
  if (mNodeD)
    VGD-=SICONOS_Lambda.getValue(mNodeD-1);
  int dimzs  = algo::spls->mDimzs - algo::spls->mV0zs;
  //-I_{DS}
  SICONOS_D.setValue(mILaw->mLine,mI->mIndex-2*algo::spls->mDimx-1,-1);
  if (mNodeG){
    ACE_DOUBLE aux=0;
    if (mMode>0){
      if (VGS-mVt>0)
	aux+=mK*(VGS-mVt);
      if (VGD-mVt>0)
	aux-=mK*(VGD-mVt);
    }else{
      if (VGS+mVt<0)
	aux-=mK*(VGS+mVt);
      if (VGD+mVt<0)
	aux+=mK*(VGD+mVt);
    }
    SICONOS_D.setValue(mILaw->mLine,mNodeG-1 ,aux);
  }
  if (mNodeS){
    ACE_DOUBLE aux=0;
    if (mMode>0){
      if (VGS-mVt > 0)
	aux=-mK*(VGS-mVt);
    }else{
      if (VGS+mVt < 0)
	aux=mK*(VGS+mVt);
    }
   SICONOS_D.setValue(mILaw->mLine,mNodeS-1 , aux);
  }
  if (mNodeD){
    ACE_DOUBLE aux=0;
    if (mMode>0){
      if (VGD-mVt>0)
	aux = mK*(VGD-mVt);
    }else{
      if (VGD+mVt<0)
	aux = -mK*(VGD+mVt);
    }
    SICONOS_D.setValue(mILaw->mLine,mNodeD-1 ,aux);// mMode*mK*(L2*(VGD-mMode*mVt)));
  }
}


componentMOS_NL2::~componentMOS_NL2(){
}

void componentMOS_NL2::print(){
  char name[ACE_CHAR_LENGTH];
  ACE_TYPE_TO_CHAR(mType,name);
  printf("component type : %s \n",name);
  if (mName)
    printf("\t Name : %s\n",mName);
  printf("\tnodePos, nodeNeg : %d %d\n",mNodePos,mNodeNeg);
  printf("\tdrain, source, gate, mode : %d %d %d %d\n",mNodeD,mNodeS,mNodeG,mMode);
  printf("\tw, vt : %f %f\n",mB,mVt);
  
}
