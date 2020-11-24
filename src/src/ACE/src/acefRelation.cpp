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

/*! \file acefRelation.cpp

*/
#ifndef ACEFRELATION_CPP
#define ACEFRELATION_CPP

#include "acefRelation.h"
//#include "aceMatrix.h"
#include "aceMatrix.h"
#include "algo.h"
#include "linearsystem.h"


#ifdef WITH_KERNEL_RELATION
//#define SICONOS_DEBUG

static int sOnlyOnce=1;
acefRelation::acefRelation():
  FirstOrderType2R()
{
}

void acefRelation::initJac(SimpleMatrix* C,SimpleMatrix* D,SimpleMatrix* B){
  mB=B;
  mC=C;
  mD=D;
  
}
void acefRelation::initialize(SP::Interaction inter)
{
  FirstOrderType2R::initialize(inter);
  unsigned int sizeY = interaction()->getSizeOfY();
  unsigned int sizeDS = interaction()->getSizeOfDS();
  SP::SiconosVector y = interaction()->y(0);
  SP::SiconosVector lambda = interaction()->lambda(0);

  
  
  workL.reset(new SimpleVector(interaction()->getSizeOfY()));
  JacXH->resize(sizeY,sizeDS);
  JacLH->resize(sizeY,sizeY);

  JacXG->resize(sizeDS,sizeDS);
  JacLG->resize(sizeDS,sizeY);

  
  *(JacXH)=*mC;
  *(JacLH)=*mD;
  *(JacLG)=*mB;

  double t0=0;
  mLastComputedSource =-1;
  for (int i=0;i <algo::spls->mzsti->size();i++)
    workL->setValue(i,algo::spls->mzsti->getValue(i));
  *lambda = *workL;

#ifdef SICONOS_DEBUG
    computeH(t0);
    cout<<"zsti\n";
    algo::spls->mzsti->display();
#endif



  
  computeG(t0);
  computeJacH(t0);
  computeJacG(t0);
  *data[r]=*data[g_alpha];
#ifdef SICONOS_DEBUG
  std::cout<<"data[r (g_alpha)] init\n";
  data[r]->display();
#endif


  
}



/*********************************/
/*y = h(X,lambda)                */
/*y = CX+D Lambda                */
/*and the NL relations           */
/*********************************/
void acefRelation::computeH(double t){

  *workX = *data[x];
  SP::SiconosVector lambda = interaction()->lambda(0);
  *workL = *lambda;
  SP::SiconosVector Heval = interaction()->relation()->Halpha();

#ifdef SICONOS_DEBUG
  std::cout<<"********         computeH at "<<t<<std::endl;
  //  CX+DL+NonLinear(X,L);
  if (sOnlyOnce){
    cout<<"acefRelation mC\n";
    mC->display();
    cout<<"acefRelation mD\n";
    mD->display();
    sOnlyOnce=0;
  }
  cout<<"acefRelation X\n";
  workX->display();
  cout<<"acefRelation L\n";
  workL->display();
#endif
  
  prod(*mC,*workX,*Heval);
  prod(*mD,*workL,*Heval,false);


  
  algo::sAlgo->computeNonLinearEquations(*workX,*workL,*Heval);

  if (t > mLastComputedSource){
#ifdef SICONOS_DEBUG
  std::cout<<"********         computeH(source) at "<<t<<std::endl;
#endif
    algo::sAlgo->preparStep(t);
    algo::spls->extractInteractionSource();
    algo::spls->updateInteractionSource();
    mLastComputedSource = t;
#ifdef SICONOS_DEBUG
    std::cout<<"B2s D2s:\n";
    algo::spls->mB2s->display();
    algo::spls->mD2s->display();
#endif

  }
  int s = algo::spls->mB2s->dimRow;
  int m = algo::spls->mD2s->dimRow;
  
  for (int i=0;i<s ;i++){
    Heval->setValue(i,Heval->getValue(i)+ algo::spls->mB2s->getValue(i));
  }
  for (int i=0;i<m ;i++){
    Heval->setValue(i+s,Heval->getValue(i+s)+ algo::spls->mD2s->getValue(i));
  }
  
#ifdef SICONOS_DEBUG
  std::cout<<"modif heval : \n";
  Heval->display();
#endif
  
}



void acefRelation::computeG(double t){
  SP::SiconosVector lambda = interaction()->lambda(0);
  *workL = *lambda;
  prod(*mB,*workL,*data[g_alpha],true);
#ifdef SICONOS_DEBUG
  std::cout<<"modif g_alpha : \n";
  data[g_alpha]->display();
#endif
  
}

  /** default function to compute jacobianH
   *  \param double : current time
   *  \param index for jacobian (0: jacobian according to x, 1 according to lambda)
   */
void acefRelation::computeJacXH(double t){
  /*Because jacobian according to x is null (in current version, may be not in futur)*/

}
void acefRelation::computeJacLH(double t){
  /*Because jacobian according to x is null (in current version, may be not in futur)*/

  SP::SiconosVector lambda = interaction()->lambda(0);
  *workX = *data[x];
  *workL = *lambda;
  algo::sAlgo->computeNonLinearJacL_H(*workX,*workL,*(JacLH));
#ifdef SICONOS_DEBUG
  std::cout<<"computeJacLH  at " <<" "<<t<<std::endl;
  (JacLH)->display();
#endif

}
	
  /** default function to compute jacobianG according to lambda
   *  \param double : current time
   *  \param index for jacobian: at the time only one possible jacobian => i = 0 is the default value .
   */
void acefRelation::computeJacLG(double t){

}
void acefRelation::computeJacXG(double t){

}


#endif
#endif
