/************************************************************************
  			linearsystem.cpp 
**************************************************************************/


#include "linearsystem.h"

//|--------x'---------|----------x--------|-------Zs------|------Zns------|

linearSystem::linearSystem(){
  mReAlloc = false;
  mDynIndex=0;
}
linearSystem::~linearSystem(){
  int i =0;
  //free unknows
  int n = mx.size();
  for (i=0; i<n; i++)
    delete mx[i];
  n = mZs.size();
  for (i=0; i<n; i++)
    delete mZs[i];
  n = mZns.size();
  for (i=0; i<n; i++)
    delete mZns[i];

  //free equations
  n=mKCL.size();
  for (i=0;i<n;i++)
    delete mKCL[i];
  n=mVD.size();
  for (i=0;i<n;i++)
    delete mVD[i];
  n=mTEN.size();
  for (i=0;i<n;i++)
    delete mTEN[i];
  n=mIND.size();
  for (i=0;i<n;i++)
    delete mIND[i];
  n=mCAP.size();
  for (i=0;i<n;i++)
    delete mCAP[i];
}

unknown* linearSystem::addinx(int type, component* c)
{
  unknown* res = new unknown(type,c);
  mx.push_back(res);
  return res;
}
unknown* linearSystem::addinZs(int type, component* c)
{
  unknown* res = new unknown(type,c);
  mZs.push_back(res);
  return res;
}
unknown* linearSystem::addinZns(int type, component* c)
{
  unknown* res = new unknown(type,c);
  mZns.push_back(res);
  return res;
}

int linearSystem::getIndexUnknown (int type,int node) {
  if (type == ACE_TYPE_V){
    return node + 2*mx.size();
  }else{
    ACE_WARNING("linearSystem::getIndexUnknown not implemented");
    return -1;
  }
  
}
equationKCL* linearSystem::KCL(int i){
  if (i >= (int)mKCL.size())
    ACE_ERROR("linearSystem::KCL");
  else
    return (equationKCL*) mKCL[i];
  return 0;
}
void linearSystem::initKCL(){
  if (mKCL.size()!=0)
    ACE_INTERNAL_ERROR("linearSystem::initKCL with mKCL not empty!");
  for (int i=0;i<mNbNodes;i++){
    mKCL.push_back(new equationKCL(i));
  }
  mKCL[0]->mAvailable=false;
}
void  linearSystem::addVUnknowns(){
 if (mZs.size()!=0)
   ACE_INTERNAL_ERROR("linearSystem::initAddVUnknowns mZs not empty");
 if (mNbNodes==0)
   ACE_INTERNAL_ERROR("linearSystem::initAddVUnknowns no nodes");
 for(int i=0;i<mNbNodes;i++){
   mZs.push_back(new unknown(ACE_TYPE_V,i));
 }
}

void linearSystem::preparForStamp(){
  int i =0;
  int nx = mx.size();
  for (i=0; i<nx; i++){
    mx[i]->mDynIndex = i;
    mx[i]->mIndex = nx+i;
    
  }
  int ns=mZs.size();
  for (i=0; i<ns; i++){
    mZs[i]->mIndex = 2*nx+i;
  }
  int Zns=mZns.size();
  for (i=0; i<Zns; i++){
    mZns[i]->mIndex = 2*nx+ns+i;
  }

  allocMemory();
}
void linearSystem::allocMemory(){
  if (mReAlloc)
    ACE_WARNING("linearSystem::allocMemory again");
  mReAlloc = true;
  mNbUnknowns = 2*mx.size()+mZs.size()+mZns.size();
  mRS = mNbUnknowns;
  int i;
  int n=mKCL.size();
  for (i=0;i<n;i++)
    mKCL[i]->allocMemory(mNbUnknowns);
  n=mVD.size();
  for (i=0;i<n;i++)
    mVD[i]->allocMemory(mNbUnknowns);
  n=mTEN.size();
  for (i=0;i<n;i++)
    mTEN[i]->allocMemory(mNbUnknowns);
  n=mIND.size();
  for (i=0;i<n;i++)
    mIND[i]->allocMemory(mNbUnknowns);
  n=mCAP.size();
  for (i=0;i<n;i++)
    mCAP[i]->allocMemory(mNbUnknowns);
}

/**
 * 
 */
bool linearSystem::isUnknown (int type, component* c,unknown **uout=0) {
  if(type == ACE_TYPE_U){
      int ii;
      for(ii=0; ii <(int) mx.size(); ii++)
	{
	  unknown* u = mx[ii];
	  
	  if (u->mType == ACE_TYPE_U && (
	      (u->mComponent->mNodePos == c->mNodePos && u->mComponent->mNodeNeg == c->mNodeNeg ) ||
	      (u->mComponent->mNodePos == c->mNodeNeg && u->mComponent->mNodeNeg == c->mNodePos ))){
	    if (uout)
	      (*uout)=u;
	    return true;
	  }
	}
      return false;
  }else{
    ACE_INTERNAL_ERROR("linearSystem::isUnknown not implemented");
    return false;
  }
  
}

equationCAP* linearSystem::addCapEquation(){
  equationCAP* res = new equationCAP();
  res->mLine = mDynIndex;
  mDynIndex++;
  mCAP.push_back(res);
  return res;
}
equationIND* linearSystem::addIndEquation(){
 equationIND* res = new equationIND();
 res->mLine = mDynIndex;
 mDynIndex++;
 mIND.push_back(res);
 return res;
}
equationVD* linearSystem::addVdEquation(){
  equationVD* res = new equationVD();
  mVD.push_back(res);
  return res;

}
equationTEN* linearSystem::addTenEquation(){
   equationTEN* eq = new equationTEN();
   mTEN.push_back(eq);
   return eq;
}

void linearSystem::addKCLinDyn(int j){
  if (j > mNbNodes){
    ACE_INTERNAL_ERROR("linearSystem::addKCLinDyn");
    return ;
  }
  mKCL[j]->mLine = mDynIndex;
  mKCL[j]->mIsDyn = true;
  mKCL[j]->mAvailable = false;
  mDynIndex++;

}
void linearSystem::print(){
  int i =0;
  int nx = mx.size();
  printf("x\n");
  for (i=0; i<nx; i++)
    mx[i]->print();
  printf("\nZs\n");
  int nzs = mZs.size();
  for (i=0; i<nzs; i++)
    mZs[i]->print();
  printf("\nZns\n");
  int nzns = mZns.size();
  for (i=0; i<nzns; i++)
    mZns[i]->print();
  printf("\n---------------------------------------\nequation\t");
  for (i=0; i<nx; i++)
    mx[i]->print();
  for (i=0; i<nx; i++)
    mx[i]->print();
  for (i=0; i<nzs; i++)
    mZs[i]->print();
  for (i=0; i<nzns; i++)
    mZns[i]->print();
  printf("\n");
  
  int n=mKCL.size();
  
  for (i=0;i<n;i++)
    mKCL[i]->print();
  n=mVD.size();
  for (i=0;i<n;i++)
    mVD[i]->print();
  n=mTEN.size();
  for (i=0;i<n;i++)
    mTEN[i]->print();
  n=mIND.size();
  for (i=0;i<n;i++)
    mIND[i]->print();
  n=mCAP.size();
  for (i=0;i<n;i++)
    mCAP[i]->print();

}

