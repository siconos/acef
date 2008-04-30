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
:componentNLINEAR(){

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
  mIb=algo::sls.addinZns(ACE_TYPE_I,this);
  sprintf(mIb->mName,"Ib_%d_%s",mNodeB,mName);
  mIc=algo::sls.addinZns(ACE_TYPE_I,this);
  sprintf(mIc->mName,"Ic_%d_%s",mNodeC,mName);
  mIndiceStartZns= mIb->mIndexInVector;
  mIndiceStartLambda= algo::sls.mDimLambda ;
  algo::sls.mDimLambda = algo::sls.mDimLambda + mDimlambda;
}
void componentBJT::stamp(){
  int ib=mIb->mIndex;
  int ic=mIc->mIndex;
  //stamp equations.
  algo::sls.KCL(mNodeB)->mCoefs[ib]-=1;
  algo::sls.KCL(mNodeC)->mCoefs[ic]-=1;
  algo::sls.KCL(mNodeE)->mCoefs[ib]+=1;
  algo::sls.KCL(mNodeE)->mCoefs[ic]+=1;

  //      |1 -1 0  |          |0|  | 0 0 0  |
  //Zns = |0 0  1  |*lambda + |0| +| 0 0 0  |(Vc,Vb,Ve)
    algo::sls.mC1l->setValue(mIndiceStartZns,mIndiceStartLambda,1);
    algo::sls.mC1l->setValue(mIndiceStartZns,mIndiceStartLambda+1,-1);
    algo::sls.mC1l->setValue(mIndiceStartZns,mIndiceStartLambda+2,0);
    algo::sls.mC1l->setValue(mIndiceStartZns+1,mIndiceStartLambda,0);
    algo::sls.mC1l->setValue(mIndiceStartZns+1,mIndiceStartLambda+1,0);
    algo::sls.mC1l->setValue(mIndiceStartZns+1,mIndiceStartLambda+2,1);
  
    //   |1 0 0 |         |-alpha  0   alpha|           |0  0      |
    //Y= |0 1 0 |*lambda+ |-alpha  0   alpha|(Vc,Vb,Ve)+|0  betha  | Zns +CST
    //   |0 0 0 |         |0       -1  1    |           |0  0      |

    // |-alpha  0   alpha|           
    // |-alpha  0   alpha|(Vc,Vb,Ve)
    // |0       -1  1    |           
  if (mNodeC){
    algo::sls.mD1zs->setValue(mIndiceStartLambda,mNodeC-1,-mAlpha);
    algo::sls.mD1zs->setValue(mIndiceStartLambda+1,mNodeC-1,-mAlpha);
  }
  if (mNodeE){
    algo::sls.mD1zs->setValue(mIndiceStartLambda,mNodeE-1,mAlpha);
    algo::sls.mD1zs->setValue(mIndiceStartLambda+1,mNodeE-1,mAlpha);
    algo::sls.mD1zs->setValue(mIndiceStartLambda+2,mNodeE-1,1);
  }
  if (mNodeB){
    algo::sls.mD1zs->setValue(mIndiceStartLambda+2,mNodeB-1,-1);
  }
  //   |1 0 0 |        
  //Y= |0 1 0 |*lambda
  //   |0 0 0 |       
  algo::sls.mD1l->setValue(mIndiceStartLambda,mIndiceStartLambda,1);
  algo::sls.mD1l->setValue(mIndiceStartLambda+1,mIndiceStartLambda+1,1);

  //|0  0      |
  //|0  betha  | Zns
  //|0  0      |
  algo::sls.mD1zns->setValue(mIndiceStartLambda+1,mIndiceStartZns+1,mBetha);

  //+CST = (0,0,Vbiais)
  algo::sls.mD1s->setValue(mIndiceStartLambda+2,mVBiais);

  
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
