/************************************************************************
  			componentmos.cpp 

Convention:
\
**************************************************************************/

#include "componentmos.h"
#include "algo.h"



componentMOS::componentMOS(dataMOS1 *d,int NbHyp)
:componentNLINEAR(){

  ACE_CHECK_IERROR(d,"componentMOS::componentMOS : Mos data null");
  ACE_CHECK_IERROR(NbHyp>=0,"componentMOS::componentMOS : NbHyp <0");
  mType = ACE_TYPE_MOS;

  mData =(*d);
  mNodeD=mData.drain;
  mNodeG=mData.gate;
  mNodeS=mData.source;
  if (NbHyp)
    mNbHyp=NbHyp;
  else
    mNbHyp=5;
  mB=mData.k/2.0;
  mMode = mData.mode;
  ACE_CHECK_IWARNING(mMode == 1 || mMode == -1,"componentMOS mode value not 1 or -1.");
  mVt=mData.vt;
  

  mNodeNeg=mNodeD;
  mNodePos=mNodeS;

  mName = mData.name;
  ACE_CHECK_ERROR(mNodeD>=0 && mNodeG>=0 && mNodeS >=0,"componentMOS::componentMOS");

  //process parameters
  double VI = 3.0;	// Power supply
  double Vt0 = mData.vt;
  double Kval = mData.k;
  double HalfK = Kval/2.0;
  double ractolpwl;
  double widthhyp;
  mHyp = (ACE_DOUBLE *) calloc(mNbHyp,sizeof(ACE_DOUBLE));
  mCoefs = (ACE_DOUBLE *) calloc(mNbHyp,sizeof(ACE_DOUBLE));


  if (mNbHyp > 1 && NbHyp)
    {
      ractolpwl = (VI-Vt0)/(2.0 * (1.0 + (sqrt(2.0)*(NbHyp - 1))));
      widthhyp = 2.0*sqrt(2.0)*ractolpwl;

      mHyp[0] = mMode*Vt0;
      mHyp[1] = mMode*Vt0+(1.0 + sqrt(2.0))*ractolpwl;

      mCoefs[0] = 2.0*ractolpwl*HalfK;
      mCoefs[1] = 2.0*widthhyp*HalfK;
      for (unsigned int i = 2;i < mNbHyp;i++) 
        {
	  mHyp[i] = mHyp[i-1] + widthhyp;
	  mCoefs[i] = 2.0*widthhyp*HalfK;
        }

    } else {
    mHyp[0] = mMode*Vt0;
    mCoefs[0] = 2.0*HalfK*(VI-Vt0)/(1.0 + sqrt(2.0));
  }
  
  mDimlambda=mNbHyp*2;
  mDimZns=1;
  mIndiceStartZns=-1;
  mIndiceStartLambda=-1;

  if (!NbHyp){
    mCoefs[0]=0.09*mB;
    mCoefs[1]=0.2238*mB;
    mCoefs[2]=0.4666*mB;
    mCoefs[3]=1.1605*mB;
    mCoefs[4]=2.8863*mB;
    mHyp[0]=mMode*mVt;
    mHyp[1]=(mMode*mVt+0.1);
    mHyp[2]=(mMode*mVt+0.2487);
    mHyp[3]=(mMode*mVt+0.6185);
    mHyp[4]=(mMode*mVt+1.5383);
  }
  
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
  for (ind=0;ind < mNbHyp;ind++){
    algo::sls.mC1l->setValue(mIndiceStartZns,mIndiceStartLambda+ind,mCoefs[ind]);
    algo::sls.mC1l->setValue(mIndiceStartZns,mIndiceStartLambda+mNbHyp+ind,-mCoefs[ind]);
  }
  
  //Y=C*zs+I*lambda+hyp

  //C*zs
  if (mNodeD){
    for(ind=0;ind < mNbHyp;ind++)
      algo::sls.mD1zs->setValue(mIndiceStartLambda+mNbHyp+ind,mNodeD-1,mMode);
  }
  if (mNodeG){
    for(ind=0;ind < mDimlambda;ind++)
      algo::sls.mD1zs->setValue(mIndiceStartLambda+ind,mNodeG-1,-mMode);
  }
  if (mNodeS){
    for(ind=0;ind < mNbHyp;ind++)
      algo::sls.mD1zs->setValue(mIndiceStartLambda+ind,mNodeS-1,mMode);
  }

  //I*lambda
  for(ind=0;ind < mDimlambda;ind++)
    algo::sls.mD1l->setValue(mIndiceStartLambda+ind,mIndiceStartLambda+ind,1);
  
  //hyp
  for(ind=0;ind < mNbHyp;ind++){
      algo::sls.mD1s->setValue(mIndiceStartLambda+ind,mHyp[ind]);
      algo::sls.mD1s->setValue(mIndiceStartLambda+mNbHyp+ind,mHyp[ind]);
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