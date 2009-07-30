/************************************************************************
  			componentmos_nl.cpp 

Convention:
\
**************************************************************************/

#include "componentmos_nl.h"
#include "algo.h"



componentMOS_NL::componentMOS_NL(dataMOS1 *d)
:component_NONE_LINEAR_NS(){

  ACE_CHECK_IERROR(d,"componentMOS_NL::componentMOS_NL : Mos data null");
  mType = ACE_TYPE_MOS_NL;
  mData =(*d);
  mNodeD=mData.drain;
  mNodeG=mData.gate;
  mNodeS=mData.source;
  mB=mData.k/2.0;
  mK=mData.k;
  mMode = mData.mode;
  ACE_CHECK_IWARNING(mMode == 1 || mMode == -1,"componentMOS_NL mode value not 1 or -1.");
  mVt=mMode*mData.vt;
  cout<<"mMode=\n"<<mMode;
  cout<<"mVt=\n"<<mVt;
  cout<<"mK=\n"<<mK;
  cout<<"mB=\n"<<mB;
  

  mNodeNeg=mNodeD;
  mNodePos=mNodeS;

  mName = mData.name;
  ACE_CHECK_ERROR(mNodeD>=0 && mNodeG>=0 && mNodeS >=0,"componentMOS_NL::componentMOS_NL");

  //process parameters
  double VI = ACE_MOS_POWER_SUPPLY;	// Power supply
  double Vt0 = mData.vt;
  double Kval = mData.k;
  double HalfK = Kval/2.0;
  double ractolpwl;
  double widthhyp;

  mDimlambda=4;
  mDimZns=0;
  mIndiceStartZns=-1;
  mIndiceStartLambda=-1;

  
}

void componentMOS_NL::addUnknowns(){
  mI=algo::spls->addinZs(ACE_TYPE_I,this);
  mIndiceStartLambda= algo::spls->mDimLambda ;
  algo::spls->mDimLambda = algo::spls->mDimLambda + mDimlambda;
}
void  componentMOS_NL::addEquations (){
  mILaw = algo::spls->addNonLinearEquation();
}

void componentMOS_NL::stamp(){
  if (ACE_FORMULATION ==  ACE_FORMULATION_MNA_V || ACE_FORMULATION ==  ACE_FORMULATION_STAMP_ONLY){
    stampMNA_V();
    return;
  }

  int ind=0;
  int i=mI->mIndex;
  //stamp equations.
  algo::spls->KCL(mNodeNeg)->mCoefs[i]-=1;
  algo::spls->KCL(mNodePos)->mCoefs[i]+=1;

  //Y=C*zs+D*lambda+cst

  //C*zs
  if (mNodeD){
    algo::spls->mD1zs->setValue(mIndiceStartLambda+1,mNodeD-1,mMode);
  }
  if (mNodeS){
    algo::spls->mD1zs->setValue(mIndiceStartLambda+3,mNodeS-1,mMode);
  } 
  if (mNodeG){
    algo::spls->mD1zs->setValue(mIndiceStartLambda+1,mNodeG-1,-mMode);
    algo::spls->mD1zs->setValue(mIndiceStartLambda+3,mNodeG-1,-mMode);
  }

  //D*lambda
  algo::spls->mD1l->setValue(mIndiceStartLambda,mIndiceStartLambda+1,-1);
  algo::spls->mD1l->setValue(mIndiceStartLambda+1,mIndiceStartLambda,1);
  algo::spls->mD1l->setValue(mIndiceStartLambda+2,mIndiceStartLambda+3,-1);
  algo::spls->mD1l->setValue(mIndiceStartLambda+3,mIndiceStartLambda+2,1);

  //cst
  algo::spls->mD1s->setValue(mIndiceStartLambda,1);
  algo::spls->mD1s->setValue(mIndiceStartLambda+1,mVt);
  algo::spls->mD1s->setValue(mIndiceStartLambda+2,1);
  algo::spls->mD1s->setValue(mIndiceStartLambda+3,mVt);
  
}
void componentMOS_NL::stampMNA_V(){

  ACE_ERROR("componentMOS_NL::stampMNA_V not yet implemented.\n");
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

void componentMOS_NL::computeNL(SiconosVector& SICONOS_X,SiconosVector& SICONOS_Lambda,SiconosVector& SICONOS_H){
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
  ACE_DOUBLE L2=SICONOS_Lambda.getValue(algo::spls->mB2zs->getDimRow()+mIndiceStartLambda+1);
  ACE_DOUBLE L4=SICONOS_Lambda.getValue(algo::spls->mB2zs->getDimRow()+mIndiceStartLambda+3);
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
  
  SICONOS_H.setValue(mILaw->mLine,-IDS+mMode*mB*(L4*(VGS-mMode*mVt)*(VGS-mMode*mVt) - L2*(VGD-mMode*mVt)*(VGD-mMode*mVt)));
  

  
}
void componentMOS_NL::computeJacNL(SiconosVector& SICONOS_X,SiconosVector& SICONOS_Lambda,SiconosMatrix& SICONOS_D){
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
 //dh/dl2
  SICONOS_D.setValue(mILaw->mLine,dimzs +mIndiceStartLambda+1,-mMode*mB*(VGD-mMode*mVt)*(VGD-mMode*mVt));
 //dh/dl4
  SICONOS_D.setValue(mILaw->mLine,dimzs+mIndiceStartLambda+3,mMode*mB*(VGS-mMode*mVt)*(VGS-mMode*mVt));
  ACE_DOUBLE L2=SICONOS_Lambda.getValue(algo::spls->mB2zs->getDimRow()+mIndiceStartLambda+1);
  ACE_DOUBLE L4=SICONOS_Lambda.getValue(algo::spls->mB2zs->getDimRow()+mIndiceStartLambda+3);
 if (mNodeG)
   SICONOS_D.setValue(mILaw->mLine,mNodeG-1 ,
		     mMode*mK*(L4*(VGS-mMode*mVt)
			   -L2*(VGD-mMode*mVt)));
 if (mNodeS)
   SICONOS_D.setValue(mILaw->mLine,mNodeS-1 ,
		      -mMode*mK*(L4*(VGS-mMode*mVt)));
 if (mNodeD)
   SICONOS_D.setValue(mILaw->mLine,mNodeD-1 ,
		      mMode*mK*(L2*(VGD-mMode*mVt)));
}


componentMOS_NL::~componentMOS_NL(){
}

void componentMOS_NL::print(){
  char name[ACE_CHAR_LENGTH];
  ACE_TYPE_TO_CHAR(mType,name);
  printf("component type : %s \n",name);
  if (mName)
    printf("\t Name : %s\n",mName);
  printf("\tnodePos, nodeNeg : %d %d\n",mNodePos,mNodeNeg);
  printf("\tdrain, source, gate, mode : %d %d %d %d\n",mNodeD,mNodeS,mNodeG,mMode);
  printf("\tw, vt : %f %f\n",mB,mVt);
  
}
