/************************************************************************
  			componentmos.cpp 

Convention:
current go from Neg node to Pos node ==> current =<0
tension = V_neg - V_pos ==> tension >=0
stamp() to write : Ax'=Bx+CZs+DZns+s
                    0 =Ex+FZs+DZns+t
**************************************************************************/

#include "componentmos.h"
#include "algo.h"



componentMOS::componentMOS(dataMOS1 *d)
:componentNLINEAR(){
  ACE_CHECK_IERROR(d,"componentMOS::componentMOS : Mos data null");
  mType = ACE_TYPE_MOS;

  mData =(*d);
  mNodeD=mData.drain;
  mNodeG=mData.gate;
  mNodeS=mData.source;
  //process parameters
  mB=mData.k/2.0;
  mMode = mData.mode;
  ACE_CHECK_IWARNING(mMode == 1 || mMode == -1,"componentMOS mode value not 1 or -1.");
  mVt=mData.vt;
  

  mNodeNeg=mNodeD;
  mNodePos=mNodeS;

  mName = mData.name;
  ACE_CHECK_ERROR(mNodeD>=0 && mNodeG>=0 && mNodeS >=0,"componentMOS::componentMOS");
  mN=5;
  mDimlambda=mN*2;
  mDimZns=1;
  mIndiceStartZns=-1;
  mIndiceStartLambda=-1;

  mCoefs = (ACE_DOUBLE *) calloc(mN,sizeof(ACE_DOUBLE));
  mHyp = (ACE_DOUBLE *) calloc(mN,sizeof(ACE_DOUBLE));
  mCoefs[0]=0.09;
  mCoefs[1]=0.2238;
  mCoefs[2]=0.4666;
  mCoefs[3]=1.1605;
  mCoefs[4]=2.8863;
  mHyp[0]=mB*mMode*mVt;
  mHyp[1]=mB*(mMode*mVt+0.1);
  mHyp[2]=mB*(mMode*mVt+0.2487);
  mHyp[3]=mB*(mMode*mVt+0.6185);
  mHyp[4]=mB*(mMode*mVt+1.5383);
  
}

void componentMOS::addUnknowns(){
  mI=algo::sls.addinZns(ACE_TYPE_I,this);
  //mIs=algo::sls.addinZns(ACE_TYPE_I,this);
  mIndiceStartZns= mI->mIndexInVector;
  mIndiceStartLambda= algo::sls.mDimLambda ;
  algo::sls.mDimLambda = algo::sls.mDimLambda + mDimlambda;
}
void componentMOS::stamp(){
  int ind=0;
  int i=mI->mIndex;
  //stamp equations.
  algo::sls.KCL(mNodeNeg)->mCoefs[i]-=mMode*1;
  algo::sls.KCL(mNodePos)->mCoefs[i]+=mMode*1;

  //Zns = B*lamdba
  for (ind=0;ind < mN;ind++){
    algo::sls.mC1l->setValue(mIndiceStartZns,mIndiceStartLambda+ind,mCoefs[ind]);
    algo::sls.mC1l->setValue(mIndiceStartZns,mIndiceStartLambda+mN+ind,-mCoefs[ind]);
  }
  
  //Y=C*zs+I*lambda+hyp

  //C*zs
  if (mNodeD){
    for(ind=0;ind < mN;ind++)
      algo::sls.mD1zs->setValue(mIndiceStartLambda+mN+ind,mNodeD-1,mMode*mB);
  }
  if (mNodeG){
    for(ind=0;ind < mDimlambda;ind++)
      algo::sls.mD1zs->setValue(mIndiceStartLambda+ind,mNodeG-1,-mMode*mB);
  }
  if (mNodeS){
    for(ind=0;ind < mN;ind++)
      algo::sls.mD1zs->setValue(mIndiceStartLambda+ind,mNodeS-1,mMode*mB);
  }

  //I*lambda
  for(ind=0;ind < mDimlambda;ind++)
    algo::sls.mD1l->setValue(mIndiceStartLambda+ind,mIndiceStartLambda+ind,1);
  
  //hyp
  for(ind=0;ind < mN;ind++){
      algo::sls.mD1s->setValue(mIndiceStartLambda+ind,0,mHyp[ind]);
      algo::sls.mD1s->setValue(mIndiceStartLambda+mN+ind,0,mHyp[ind]);
    }

  
}

componentMOS::~componentMOS(){
  free(mCoefs);
  free(mHyp);
}

void componentMOS::print(){
  char name[ACE_CHAR_LENGTH];
  ACE_TYPE_TO_CHAR(mType,name);
  printf("component type : %s \n",name);
  if (mName)
    printf("\t Name : %s\n",mName);
  printf("\tnodePos, nodeNeg : %d %d\n",mNodePos,mNodeNeg);
  printf("\tdrain, source, gate, mode : %d %d %d %d\n",mNodeD,mNodeS,mNodeG,mMode);
  printf("\tw, vt : %f %f\n",mB,mVt);
  
}
