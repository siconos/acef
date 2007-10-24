/************************************************************************
  			linearsystem.cpp 
**************************************************************************/


#include "linearsystem.h"

using namespace std;
//|--------x'---------|----------x--------|-------Zs------|------Zns------|

linearSystem::linearSystem(){
  mReAlloc = false;
  mDynIndex=0;
  mNbEquations=0;
  mNbDynEquations=0;
  mA=0;
  mB=0;
  mC=0;
  mD=0;
  ms=0;

  mA1x=0;
  mA1zs=0;
  mA1zns=0;
  mA1s=0;

}
void linearSystem::allocMatrix(){
  ACE_CHECK_IERROR(!mA,"linearSystem::buildABCDs : mA not NULL");
  ACE_CHECK_IERROR(!mA,"linearSystem::buildABCDs : mB not NULL");

  int nx = mx.size();
  if (nx == 0 ){
    ACE_MESSAGE("No Dynamique\n");
  }else{
    mA = new aceMatrix(nx,nx);
    mB= new aceMatrix(nx,nx);
    int nzs=mZs.size();
    int nzns=mZns.size();
    if (nzs)
      mC=new aceMatrix(nx,nzs);
    if (nzns)
      mD=new aceMatrix(nx,nzns);
    ms=new aceMatrix(nx,1);
    
    mA1x = mB;
    mA1zs= mC;
    mA1zns= mD;
    mA1s=ms;
  }
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
  if (mA)
    delete mA;
  if (mB)
    delete mB;
  if (mC)
    delete mC;
  if (mD)
    delete mD;
  if (ms)
    delete ms;

  
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
  ACE_CHECK_IERROR(mKCL.size()==0,"linearSystem::initKCL with mKCL not empty!");
  for (int i=0;i<mNbNodes;i++){
    mKCL.push_back(new equationKCL(i));
  }
  mKCL[0]->mAvailable=false;
  mNbEquations+=mNbNodes-1;
}

void  linearSystem::addVUnknowns(){
  ACE_CHECK_IERROR(mZs.size()==0,"linearSystem::initAddVUnknowns mZs not empty");
  ACE_CHECK_IERROR(mNbNodes,"linearSystem::initAddVUnknowns no nodes");
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
  ACE_CHECK_IERROR(!mReAlloc,"linearSystem::allocMemory again");
  mReAlloc = true;
  mNbUnknowns = 2*mx.size()+mZs.size()+mZns.size();
  mRS = mNbUnknowns;
  int i;
  int n=mKCL.size();
  for (i=0;i<n;i++)
    mKCL[i]->allocMemory(mNbUnknowns+1);
  n=mVD.size();
  for (i=0;i<n;i++)
    mVD[i]->allocMemory(mNbUnknowns+1);
  n=mTEN.size();
  for (i=0;i<n;i++)
    mTEN[i]->allocMemory(mNbUnknowns+1);
  n=mIND.size();
  for (i=0;i<n;i++)
    mIND[i]->allocMemory(mNbUnknowns+1);
  n=mCAP.size();
  for (i=0;i<n;i++)
    mCAP[i]->allocMemory(mNbUnknowns+1);
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
  mNbEquations++;
  mNbDynEquations++;
  res->mLine = mDynIndex;
  mDynIndex++;
  mCAP.push_back(res);
  return res;
}
equationIND* linearSystem::addIndEquation(){
 equationIND* res = new equationIND();
 mNbEquations++;
 mNbDynEquations++;
 res->mLine = mDynIndex;
 mDynIndex++;
 mIND.push_back(res);
 return res;
}
equationVD* linearSystem::addVdEquation(){
  equationVD* res = new equationVD();
  mNbEquations++;
  mVD.push_back(res);
  return res;

}
equationTEN* linearSystem::addTenEquation(){
   equationTEN* eq = new equationTEN();
   mNbEquations++;
   mTEN.push_back(eq);
   return eq;
}

void linearSystem::addKCLinDyn(int j){
  ACE_CHECK_IERROR(j <= mNbNodes,"linearSystem::addKCLinDyn");
  if (!mKCL[j]->mIsDyn)
    mNbDynEquations++;
  else
    ACE_INTERNAL_WARNING("linearSystem::addKCLinDyn kcl in dyn.");
  
  mKCL[j]->mLine = mDynIndex;
  mKCL[j]->mIsDyn = true;
  mKCL[j]->mAvailable = false;
  mKCL[j]->mLine=mDynIndex;
  mDynIndex++;

}

//fill matrix A1
void linearSystem::computedxdt(){
  allocMatrix();
  buildABCDs();
  if(!mA)
    return;
  printf("A:\n");
  mA->display();
  try{
    mA->PLUInverseInPlace();
  }
  catch(SiconosException e)
    {
      std::cout << e.report() << endl;
      ACE_INTERNAL_ERROR("linearSystem::computedxdt");
    }
  catch(...)
    {
      std::cout << "Exception caught." << endl;
      ACE_INTERNAL_ERROR("linearSystem::computedxdt");
    }
  printABCDs();
  try{
    if (mA)
      (*mA1x)=prod(*mA,*mB);
    if (mC)
      (*mA1zs)=prod(*mA,*mC);
    if (mD)
      (*mA1zns)=prod(*mA,*mD);
    
    (*mA1s)=prod(*mA,*ms);
  }
  catch(SiconosException e)
    {
      std::cout << e.report() << endl;
      ACE_INTERNAL_ERROR("linearSystem::computedxdt");
    }
  catch(...)
    {
      std::cout << "Exception caught." << endl;
      ACE_INTERNAL_ERROR("linearSystem::computedxdt");
    }
  printA1();

}

void linearSystem::getlinefromdxdt(int line, ACE_DOUBLE * coefs){
  if (!mA)
    return;
  ACE_CHECK_IERROR(coefs && line < mNbDynEquations && line >=0,"linearSystem::getlinefromdxdt");
  int idStart =0;
  int nx = mx.size();
  int nzs = mZs.size();
  int nzns = mZns.size();
  int i=0;
  for (i =idStart; i < nx;i++){
    coefs[i]= (*mA1x)(line,i-idStart);
  }
  idStart+=nx;
  for (i =idStart; i < idStart+nzs;i++){
    coefs[i]= (*mA1zs)(line,i-idStart);
  }
  idStart+=nzs;
  for (i =idStart; i < idStart+nzns;i++){
    coefs[i]= (*mA1zns)(line,i-idStart);
  }
  idStart+=nzns;
  ACE_CHECK_IERROR(idStart+nx == mNbUnknowns,"linearSystem::getlinefromdxdt idStart == mNbUnknowns");
  coefs[idStart]=(*mA1s)(line,0);

}


void linearSystem::buildABCDs(){
  int nx = mx.size();
  int nzs = mZs.size();
  int nzns = mZns.size();
  int istart=0;
  ACE_CHECK_IERROR( nx == mNbDynEquations,"linearSystem::buildABCDs : mx dim != nbDynEquations");

  //BUILD A
  if (mA)
    extractDynBockInMat(mA,istart,istart+nx);
     
  istart+=nx;
  //BUILD B
  if (mB)
    extractDynBockInMat(mB,istart,istart+nx);
  istart+=nx;

  //BUILD C
  if (mC)
    extractDynBockInMat(mC,istart,istart+nzs);    
  istart+=nzs;

  //BUILD D
  if (mD)
    extractDynBockInMat(mD,istart,istart+nzns);    
  istart+=nzns;

  //BUILD s
  ACE_CHECK_IERROR(istart == mNbUnknowns,"linearSystem::buildABCDs istart == mNbUnknowns");
  extractDynBockInMat(ms,istart,istart+1);
}



//copy in m the double contains in dynamique equation from IndexBegin to IndexEnd
void linearSystem::extractDynBockInMat(aceMatrix * m, int IndexBegin, int IndexEnd){
  ACE_CHECK_IERROR(m,"linearSystem::extractDynBockInMat m null");
  ACE_CHECK_IERROR(IndexBegin>=0,"linearSystem::extractDynBockInMat IndexBegin");
  ACE_CHECK_IERROR(IndexEnd>=IndexBegin && IndexEnd < mNbUnknowns+2 ,"linearSystem::extractDynBockInMat IndexEnd");
  
  
    int nbInd = mIND.size();
    int i,j=0;
    int line=0;
    for (i=0;i<nbInd;i++){
      ACE_DOUBLE * coefs=mIND[i]->mCoefs;
      line = mIND[i]->mLine;
      ACE_CHECK_IERROR(line>=0 && coefs,"linearSystem::extractDynBockInMat Ind");
      for(j=IndexBegin;j<IndexEnd;j++){
	m->setValue(line,j-IndexBegin,coefs[j]);
      }
    }
    int nbCap = mCAP.size();
    for (i=0;i<nbCap;i++){
      ACE_DOUBLE * coefs=mCAP[i]->mCoefs;    
      line = mCAP[i]->mLine;
      ACE_CHECK_IERROR(line>=0 && coefs,"linearSystem::extractDynBockInMat  Cap");
      for(j=IndexBegin;j<IndexEnd;j++){
	m->setValue(line,j-IndexBegin,coefs[j]);
      }
    }
    int nbKcl = mKCL.size();
    for (i=0;i<nbKcl;i++){
      if (mKCL[i]->mIsDyn){
	ACE_DOUBLE * coefs=mKCL[i]->mCoefs;    
	line = mKCL[i]->mLine;
	ACE_CHECK_IERROR(line>=0 && coefs,"linearSystem::extractDynBockInMat  Kcl");
	for(j=IndexBegin;j<IndexEnd;j++){
	  m->setValue(line,j-IndexBegin,coefs[j]);
	}
      }
    }
    
}









void linearSystem::printEquations(){
  int i,n =0;
  int nx = mx.size();
  printf("--->linearSystem with %d equations whose %d dynamic equations.\n",mNbEquations,mNbDynEquations);
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
  printf("\n---------------------------------------\nequation");
  for (i=0; i<nx; i++)
    mx[i]->print();
  for (i=0; i<nx; i++)
    mx[i]->print();
  for (i=0; i<nzs; i++)
    mZs[i]->print();
  for (i=0; i<nzns; i++)
    mZns[i]->print();
  printf("\n");
  
  //print dyn equation
  n=mKCL.size();
  for (i=0;i<n;i++)
    if (mKCL[i]->mIsDyn)
      mKCL[i]->print();
  n=mCAP.size();
  for (i=0;i<n;i++)
    mCAP[i]->print();  
  n=mIND.size();
  for (i=0;i<n;i++)
    mIND[i]->print();
  
  for (i=0;i<n;i++)
    if (!mKCL[i]->mIsDyn)
      mKCL[i]->print();
  n=mVD.size();
  for (i=0;i<n;i++)
    mVD[i]->print();
  n=mTEN.size();
  for (i=0;i<n;i++)
    mTEN[i]->print();
}

void linearSystem::printABCDs(){
  printf("inv A:\n");
  if (mA)
    mA->display();
  printf("B:\n");
  if (mB)
    mB->display();
  printf("C:\n");
  if (mC)
    mC->display();
  printf("D:\n");
  if (mD)
    mD->display();
  printf("s:\n");
  if (ms)
    ms->display();
 }
void linearSystem::printA1(){
  printf("system x'=A1x+A1zs+A1zns+s\n");
  printf("A1x:\n");
  if (mA1x)
    mA1x->display();
  printf("A1zs:\n");
  if (mA1zs)
    mA1zs->display();
  printf("A1zns:\n");
  if (mA1zns)
    mA1zns->display();
  printf("A1s:\n");
  if (mA1s)
    mA1s->display();
 }
