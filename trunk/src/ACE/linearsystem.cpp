/************************************************************************
  			linearsystem.cpp 
**************************************************************************/


#include "linearsystem.h"
#include <fstream>

using namespace std;
//|--------x'---------|----------x--------|-------Zs------|------Zns------|

linearSystem::linearSystem(){


  
  mReAlloc = false;
  mDynIndex=0;
  mNbNodes=0;
  mNbUnknowns=0;
  mNbEquations=0;
  mNbDynEquations=0;
  mNbNonDynEquations=0;
  mA=0;
  mB=0;
  mC=0;
  mD=0;
  ms=0;

  mA1x=0;
  mA1zs=0;
  mA1zns=0;
  mA1s=0;

  mB1x=0;
  mB1zs=0;
  mB1zns=0;
  mB1s=0;

  mC1x=0;
  mC1zs=0;
  mC1l=0;
  mC1s=0;



}
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////MEMORY

void linearSystem::allocA1Matrix(){
  ACE_CHECK_IERROR(!mA,"linearSystem::buildABCDs : mA not NULL");
  ACE_CHECK_IERROR(!mA,"linearSystem::buildABCDs : mB not NULL");

  int nx = mx.size();
  if (nx != 0 ){
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
void linearSystem::freeA1Matrix(){

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

//0 = B1x*x + B1zs*Zs + B1zns*Zns + B1s
void linearSystem::allocB1Matrix(){
  if (mNbNonDynEquations>0){
    int nx = mx.size();
    int nzs= mZs.size();
    int nzns=mZns.size();
    if (nx)
      mB1x=new aceMatrix(mNbNonDynEquations,nx);
    if (nzs)
      mB1zs=new aceMatrix(mNbNonDynEquations,nzs);
    if (nzns)
      mB1zns=new aceMatrix(mNbNonDynEquations,nzns);
    mB1s=new aceMatrix(mNbNonDynEquations,1);
  }
}
void linearSystem::freeB1Matrix(){
  if (mB1x)
    delete mB1x;
  if (mB1zs)
    delete mB1zs;
  if (mB1zns)
    delete mB1zns;
  if (mB1s)
    delete mB1s;
  
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
void linearSystem::allocC1Matrix(){
  int nx = mx.size();
  int nzs= mZs.size();
  int nzns=mZns.size();
  if (!nZns)
    return;
  if (nx)
    mC1x=new aceMatrix(nzns,nx);
  if (nzs)
    mC1zs=new aceMatrix(nzns,nzs);
  if(??)
    mC1l = new aceMatrix(nzns,???);
  mC1s = new aceMatrix(nzns,1);

}
void linearSystem::freeC1Matrix(){
  if (mC1x)
    delete mC1x;
  if (mC1zs)
    delete mC1zs;
  if (mC1l)
    delete mC1l;
  if (mC1s)
    delete mC1s;
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
  freeA1Matrix();
  freeB1Matrix();
  freeC1Matrix();
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////UNKNOWN
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
void  linearSystem::addVUnknowns(){
  ACE_CHECK_IERROR(mZs.size()==0,"linearSystem::initAddVUnknowns mZs not empty");
  ACE_CHECK_IERROR(mNbNodes,"linearSystem::initAddVUnknowns no nodes");
 for(int i=0;i<mNbNodes;i++){
   mZs.push_back(new unknown(ACE_TYPE_V,i));
 }
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
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////EQUATIONS
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

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////STAMP

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





///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////DYNAMICS EQUATIONS  : x'=A1x*x + A1zs*Zs + A1zns*Zns + A1s
//fill matrix: x'=A1x * mx + A1zs * mZs + A1zns * mZns; 
void linearSystem::computedxdt(){
  if (mx.size() == 0){
    ACE_MESSAGE("No Dynamic\n");
    return;
  }
  allocA1Matrix();
  buildABCDs();
  if(!mA)
    return;
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
  cout<<"inv A:\n-----\n";
  cout<<(*mA);

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
  ACE_CHECK_IERROR(idStart+nx == mNbUnknowns,"linearSystem::getlinefromdxdt idStart +nx == mNbUnknowns");
  coefs[idStart]=(*mA1s)(line,0);

}
//Ax'=Bx+CZs+DZns+s
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
  printABCDs();

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

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////LINEAR SYSTEM : 0 = B1x*x + B1zs*Zs + B1zns*Zns + B1s
void linearSystem::buildLinearSystem(){
  mNbNonDynEquations = mNbEquations - mNbDynEquations;
  ACE_CHECK_IERROR(mNbNonDynEquations >=0,"linearSystem::allocB1Matrix, mNbNonDynEquations <0.");
  allocB1Matrix();

  int nx = mx.size();
  int nzs = mZs.size();
  int nzns = mZns.size();
  int istart = nx;
  //BUILD B1x
  if (mB1x)
    extractNonDynBockInMat(mB1x,istart,istart+nx);
  istart+=nx;
  
  //BUILD B1zs
  if (mB1zs)
    extractNonDynBockInMat(mB1zs,istart,istart+nzs);
  istart+=nzs;

  //BUILD B1zns
  if (mB1zns)
    extractNonDynBockInMat(mB1zns,istart,istart+nzns);    
  istart+=nzns;

  //BUILD s
  ACE_CHECK_IERROR(istart == mNbUnknowns,"linearSystem::buildLinearSystem istart == mNbUnknowns");
  extractNonDynBockInMat(mB1s,istart,istart+1);
}


//copy in m the double contains in none dynamique equation from IndexBegin to IndexEnd
void linearSystem::extractNonDynBockInMat(aceMatrix * m, int IndexBegin, int IndexEnd){
  ACE_CHECK_IERROR(m,"linearSystem::extractNonDynBockInMat m null");
  int nx=mx.size();
  ACE_CHECK_IERROR(IndexBegin>=nx,"linearSystem::extractNonDynBockInMat IndexBegin");
  ACE_CHECK_IERROR(IndexEnd>=IndexBegin && IndexEnd < mNbUnknowns+2 ,"linearSystem::extractNonDynBockInMat IndexEnd");
  
  int line=0;
  int i,j;

  int nbKcl = mKCL.size();
  for (i=1;i<nbKcl;i++){
    if (!mKCL[i]->mIsDyn){
      ACE_DOUBLE * coefs=mKCL[i]->mCoefs;    
      ACE_CHECK_IERROR(coefs,"linearSystem::extractNonDynBockInMat  Kcl");
      for(j=IndexBegin;j<IndexEnd;j++){
	m->setValue(line,j-IndexBegin,coefs[j]);
      }
      line++;
    }
  }
  int nVD = mVD.size();
  for (i=0;i<nVD;i++){
      ACE_DOUBLE * coefs=mVD[i]->mCoefs;    
      ACE_CHECK_IERROR(coefs,"linearSystem::extractNonDynBockInMat  VD");
      for(j=IndexBegin;j<IndexEnd;j++){
	m->setValue(line,j-IndexBegin,coefs[j]);
      }
      line++;
  }

  int nTEN = mTEN.size();
  for (i=0;i<nTEN;i++){
      ACE_DOUBLE * coefs=mTEN[i]->mCoefs;    
      ACE_CHECK_IERROR(coefs,"linearSystem::extractNonDynBockInMat  TEN");
      for(j=IndexBegin;j<IndexEnd;j++){
	m->setValue(line,j-IndexBegin,coefs[j]);
      }
      line++;
  }

  ACE_CHECK_IERROR(mNbNonDynEquations == line,"linearSystem::extractNonDynBockInMat number equation!");
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////PRINT FONCTIONS
void linearSystem::printEquations(ostream& os ){
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
    mx[i]->printdev();
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

  //non dyn equation
  n=mKCL.size();
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

void linearSystem::printABCDs(ostream& os ){
  os <<"Ax'=Bx+CZs+DZns+s\n";
  os <<"-----------------\n";
  os <<"A:\n";
  if (mA)
    os << (*mA);
  os <<"B:\n";
  if (mB)
    os<< (*mB);
  os<<"C:\n";
  if (mC)
    os <<(*mC);
  os<<"D:\n";
  if (mD)
    os <<(*mD);
  os <<"s:\n";
  if (ms)
    os <<(*ms);
 }
void linearSystem::printA1(ostream& os ){
  os <<"system x'=A1x+A1zs+A1zns+s\n";
  os <<"A1x:\n";
  if (mA1x)
    os <<(*mA1x);
  os <<"A1zs:\n";
  if (mA1zs)
    os <<(*mA1zs);
  os <<"A1zns:\n";
  if (mA1zns)
    os <<(*mA1zns);
  os <<"A1s:\n";
  if (mA1s)
    os << (*mA1s);
 }
void linearSystem::printB1(ostream& os ){
  os <<"0 = B1x*x + B1zs*Zs + B1zns*Zns + B1s\n";
  os <<"B1x:\n";
  if (mB1x)
    os <<(*mB1x);
  os <<"B1zs:\n";
  if (mB1zs)
    os <<(*mB1zs);
  os <<"B1zns:\n";
  if (mB1zns)
    os <<(*mB1zns);
  os <<"B1s:\n";
  if (mB1s)
    os << (*mB1s);
 }

void linearSystem::printSystemInTabFile(char * file){
  ofstream pout(file);
  printA1(pout);
  printB1(pout);
  pout.close();
}
