/************************************************************************
  			node.cpp 
**************************************************************************/


#include "node.h"
#include "algo.h"
#include "linearsystem.h"
#include "componentcap.h"


node::node(int i){
//   mXCoef=0;
//   mXBuffer=0;
//   mZsCoef=0;
//   mZsBuffer=0;
  mId = i;
  mC=0;
  mIndexInZs=-1;

}
node::~node(){
//   if (mXCoef)
//     free(mXCoef);
//   if (mXBuffer)
//     free(mXBuffer);
//   if (mZsCoef)
//     free(mZsCoef);
//   if (mZsBuffer)
//     free(mZsBuffer);
//   mXCoef=0;
//   mXBuffer=0;
//   mZsCoef=0;
//   mZsBuffer=0;
}

// void node::alloc(int dimX,int dimZs){
//   ACE_CHECK_IERROR(mXCoef!=0,"node::alloc");
//   ACE_CHECK_IERROR(mZsCoef!=0,"node::alloc");
//   mDimX=dimX;
//   mDimZs=dimZs;
//   mXCoef=(ACE_DOUBLE* )calloc(dimX*sizeof(ACE_DOUBLE));
//   mXBuffer=(ACE_DOUBLE* )calloc(dimX*sizeof(ACE_DOUBLE));
//   mZsCoef=(ACE_DOUBLE* )calloc(dimZs*sizeof(ACE_DOUBLE));
//   mZsBuffer=(ACE_DOUBLE* )calloc(dimZs*sizeof(ACE_DOUBLE));
// }
// void node::VProd(ACE_DOUBLE v){
//   memcpy(mXBuffer,mXCoef,mDimX*sizeof(ACE_DOUBLE));
//   memcpy(mZsBuffer,mZsCoef,mDimZs*sizeof(ACE_DOUBLE));
//   for(int i=0;i<mDimX;i++)
//     mXBuffer[i]*=v;
//   for(int i=0;i<mDimZs;i++)
//     mZsBuffer[i]*=v;
// }
bool node::isVUnknown(){
  if (mC)
    return false;
  return true;
}
void node::stampV(ACE_DOUBLE v,int lin,aceMatrix * X, aceMatrix * Zs){
  if (!mC){
    Zs->incValue(lin,mIndexInZs-1,v);
    return;
  }
  //In capacitor branche :
  //mU=nodeNeg-nodePos
  if (mId == mC->mNodePos){
    //nodePos=-mU+nodeNeg
    X->incValue(lin,mC->mU->mDynIndex,-v);
    if (mC->mNodeNeg)
      algo::spls->mNodes[mC->mNodeNeg]->stampV(v,lin,X,Zs);
  }else if (mId == mC->mNodeNeg){
    //nodeNeg=nodePos+mU;
    //X[mC->mU->mDynIndex]+=v;
    X->incValue(lin,mC->mU->mDynIndex,v);
    if (mC->mNodePos)
      algo::spls->mNodes[mC->mNodePos]->stampV(v,lin,X,Zs);
  }else
    ACE_ERROR("node::preparStamp, id out of capa!");
}
void node::stampV(ACE_DOUBLE v,ACE_DOUBLE * X, ACE_DOUBLE * Zs){
  if (!mC){
    Zs[mIndexInZs]+=v;
    return;
  }
  //In capacitor branche :
  //mU=nodeNeg-nodePos
  if (mId == mC->mNodePos){
    //nodePos=-mU+nodeNeg
    X[mC->mU->mDynIndex]-=v;
    algo::spls->mNodes[mC->mNodeNeg]->stampV(v,X,Zs);
  }else if (mId == mC->mNodeNeg){
    //nodeNeg=nodePos+mU;
    X[mC->mU->mDynIndex]+=v;
    algo::spls->mNodes[mC->mNodePos]->stampV(v,X,Zs);
  }else
    ACE_ERROR("node::preparStamp, id out of capa!");
}
void node::print () {
  
  
}


double node::getValue(aceVector * X, aceVector * Zs){
  if (!mId)
    return 0;
  if (!mC){
    return Zs->getValue(mIndexInZs-1);
  }
  //In capacitor branche :
  //mU=nodeNeg-nodePos
  ACE_DOUBLE aux=0;
  if (mId == mC->mNodePos){
    //nodePos=-mU+nodeNeg
    aux = -X->getValue(mC->mU->mDynIndex);
    aux += algo::spls->mNodes[mC->mNodeNeg]->getValue(X,Zs);
  }else if (mId == mC->mNodeNeg){
    //nodeNeg=nodePos+mU;
    aux = X->getValue(mC->mU->mDynIndex);
    aux += algo::spls->mNodes[mC->mNodePos]->getValue(X,Zs);
  }else
    ACE_ERROR("node::preparStamp, id out of capa!");
  return aux;

  
}
