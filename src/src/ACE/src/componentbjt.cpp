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

/*! \file componentbjt.cpp

*/
/************************************************************************
  			componentbjt.cpp 
Bipolar Junction Transistor

      /c
 b  |/
----|
    |
    |\
      \e
       >

compementarity model:
---------------------

*)Amplification:

Ic=0+U-V
j=-alpha*(Vc-Ve) + U
k=V+betha*Ib-alpha*(Vc-Ve)

*)Condition Amplification:

Ib=W
l=-(Vb-Ve-Vbiais)

*)voltage repartition ???

Ub=Vb = Vc  + Vbiais -Z
m = Z - (Vc)

ace formulation:
----------------

lambda=(U,V,W)
Y=(j,k,l)


Zns = (Ic,Ib)

      |1 -1 0  |          |0|  | 0 0 0  |
Zns = |0 0  1  |*lambda + |0| +| 0 0 0  |(Vc,Vb,Ve)

   |1 0 0 |         |-alpha  0   alpha|           |0  0      |
Y= |0 1 0 |*lambda+ |-alpha  0   alpha|(Vc,Vb,Ve)+|0  betha  | Zns +CST
   |0 0 0 |         |0       -1  1    |           |0  0      |
   
CST=(0,0,Vbiais,0)

0<lambda per Y>0
Convention:
\
**************************************************************************/

#include "componentbjt.h"
#include "algo.h"
#include "unknown.h"



componentBJT::componentBJT(dataBJT *d)
:component_LINEAR_NS(){

  ACE_CHECK_IERROR(d,"componentBJT::componentBJT : Bjt data null");
  mType = ACE_TYPE_BJT;

  mData =(*d);
  mNodeC=mData.collector;
  mNodeB=mData.base;
  mNodeE=mData.emitor;
  mMode = mData.mode;
  ACE_CHECK_IWARNING(mMode == 1 || mMode == -1,"componentBJT mode value not 1 or -1.");
  mBetha = 100;
  mAlpha = 10000;
  mVBiais = 0.7;

  mName = mData.name;
  
  mDimlambda=3;
  mDimZns=2;
  mIndiceStartZns=-1;
  mIndiceStartLambda=-1;

}

void componentBJT::addUnknowns(){
  mIb=algo::spls->addinZns(ACE_TYPE_I,this);
  sprintf(mIb->mName,"Ib_%d_%s",mNodeB,mName);
  mIc=algo::spls->addinZns(ACE_TYPE_I,this);
  sprintf(mIc->mName,"Ic_%d_%s",mNodeC,mName);
  mIndiceStartZns= mIb->mIndexInVector;
  mIndiceStartLambda= algo::spls->mDimLambda ;
  algo::spls->mDimLambda = algo::spls->mDimLambda + mDimlambda;
}
void componentBJT::stamp(){
  int ib=mIb->mIndex;
  int ic=mIc->mIndex;
  //stamp equations.
  algo::spls->KCL(mNodeB)->mCoefs[ib]-=1;
  algo::spls->KCL(mNodeC)->mCoefs[ic]-=1;
  algo::spls->KCL(mNodeE)->mCoefs[ib]+=1;
  algo::spls->KCL(mNodeE)->mCoefs[ic]+=1;

  //      |1 -1 0  |          |0|  | 0 0 0  |
  //Zns = |0 0  1  |*lambda + |0| +| 0 0 0  |(Vc,Vb,Ve)
    algo::spls->mC1l->setValue(mIndiceStartZns,mIndiceStartLambda,1);
    algo::spls->mC1l->setValue(mIndiceStartZns,mIndiceStartLambda+1,-1);
    algo::spls->mC1l->setValue(mIndiceStartZns,mIndiceStartLambda+2,0);
    algo::spls->mC1l->setValue(mIndiceStartZns+1,mIndiceStartLambda,0);
    algo::spls->mC1l->setValue(mIndiceStartZns+1,mIndiceStartLambda+1,0);
    algo::spls->mC1l->setValue(mIndiceStartZns+1,mIndiceStartLambda+2,1);
  
    //   |1 0 0 |         |-alpha  0   alpha|           |0  0      |
    //Y= |0 1 0 |*lambda+ |-alpha  0   alpha|(Vc,Vb,Ve)+|0  betha  | Zns +CST
    //   |0 0 0 |         |0       -1  1    |           |0  0      |

    // |-alpha  0   alpha|           
    // |-alpha  0   alpha|(Vc,Vb,Ve)
    // |0       -1  1    |           
  if (mNodeC){
    algo::spls->mD1zs->setValue(mIndiceStartLambda,mNodeC-1,-mAlpha);
    algo::spls->mD1zs->setValue(mIndiceStartLambda+1,mNodeC-1,-mAlpha);
  }
  if (mNodeE){
    algo::spls->mD1zs->setValue(mIndiceStartLambda,mNodeE-1,mAlpha);
    algo::spls->mD1zs->setValue(mIndiceStartLambda+1,mNodeE-1,mAlpha);
    algo::spls->mD1zs->setValue(mIndiceStartLambda+2,mNodeE-1,1);
  }
  if (mNodeB){
    algo::spls->mD1zs->setValue(mIndiceStartLambda+2,mNodeB-1,-1);
  }
  //   |1 0 0 |        
  //Y= |0 1 0 |*lambda
  //   |0 0 0 |       
  algo::spls->mD1l->setValue(mIndiceStartLambda,mIndiceStartLambda,1);
  algo::spls->mD1l->setValue(mIndiceStartLambda+1,mIndiceStartLambda+1,1);

  //|0  0      |
  //|0  betha  | Zns
  //|0  0      |
  algo::spls->mD1zns->setValue(mIndiceStartLambda+1,mIndiceStartZns+1,mBetha);

  //+CST = (0,0,Vbiais)
  algo::spls->mD1s->setValue(mIndiceStartLambda+2,mVBiais);

  
}

componentBJT::~componentBJT(){
}

void componentBJT::print(){
  char name[ACE_CHAR_LENGTH];
  ACE_TYPE_TO_CHAR(mType,name);
  printf("component type : %s \n",name);
  if (mName)
    printf("\t Name : %s\n",mName);
  printf("\tBase, Collector, Emitor , mode: %d %d %d %d \n",mNodeB,mNodeC,mNodeE,mMode);
  printf("\tbetha, vbiais : %f %f\n",mBetha,mVBiais);
  
}
