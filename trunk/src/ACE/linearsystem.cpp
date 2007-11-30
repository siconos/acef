/************************************************************************
  			linearsystem.cpp 
**************************************************************************/


#include "linearsystem.h"
#include <fstream>
#include "SimpleMatrix.h"

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
  mDimLambda=0;
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

  mD1x=0;
  mD1zs=0;
  mD1zns=0;
  mD1l=0;
  mD1s=0;

  mR=0;
  mA2x=0;
  mA2zs=0;
  mA2s=0;
  mA2sti=0;
  
  mB2x=0;
  mB2zs=0;
  mB2l=0;
  mB2s=0;

  mD2x=0;
  mD2zs=0;
  mD2l=0;
  mD2s=0;

  mDimx=0;
  mDimzs=0;
  mDimzns=0;

  mxti=0;
  mzsti=0;
  mznsti=0;
  mxfree=0;

  mW=0;
  mD3l=0;
  mD3zs=0;
  mB3l=0;
  mB3zs=0;
  mPfree=0;
  mQfree=0;
  mMLCP=0;

  mTheta = 0.5;
  mThetap = 0.5;
  mH = 1;
}
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////MEMORY


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


  //Matrix allocation
  allocA1Matrix();
  mNbNonDynEquations = mNbEquations - mNbDynEquations;
  ACE_CHECK_IERROR(mNbNonDynEquations == mDimzs-1,"linearSystem::preparForStamp, mNbNonDynEquations and mDimzs-1 not coherent");
  ACE_CHECK_IERROR(mNbNonDynEquations >=0,"linearSystem::preparForStamp, mNbNonDynEquations <0.");
  allocB1Matrix();
  allocC1Matrix();
  allocD1Matrix();
  
}
//x'=A1x*x + A1zs*Zs + A1zns*Zns + A1s
void linearSystem::allocA1Matrix(){
  ACE_CHECK_IERROR(!mA,"linearSystem::buildABCDs : mA not NULL");
  ACE_CHECK_IERROR(!mB,"linearSystem::buildABCDs : mB not NULL");

  
  if (mDimx != 0 ){
    mA = new aceMatrix(mDimx,mDimx);
    mB= new aceMatrix(mDimx,mDimx);
    
    if (mDimzs-1>0)
      mC=new aceMatrix(mDimx,mDimzs-1);
    if (mDimzns)
      mD=new aceMatrix(mDimx,mDimzns);
    ms=new aceMatrix(mDimx,1);
    
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
  if (!mNbNonDynEquations)
    return;
  ACE_CHECK_IERROR(mB1zs ==0,"linearSystem::allocB1Matrix again");
  if (mDimx)
    mB1x=new aceMatrix(mNbNonDynEquations,mDimx);
  if (mDimzs-1>0)
    mB1zs=new aceMatrix(mNbNonDynEquations,mDimzs-1);
  if (mDimzns)
    mB1zns=new aceMatrix(mNbNonDynEquations,mDimzns);
  mB1s=new aceMatrix(mNbNonDynEquations,1);
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
void linearSystem::allocC1Matrix(){
  ACE_CHECK_IERROR(mC1zs == 0,"linearSystem::allocC1Matrix");
  if (!mDimzns)
    return;
  if (mDimx)
    mC1x=new aceMatrix(mDimzns,mDimx);
  if (mDimzs-1>0)
    mC1zs=new aceMatrix(mDimzns,mDimzs-1);
  if(mDimLambda)
    mC1l = new aceMatrix(mDimzns,mDimLambda);
  mC1s = new aceMatrix(mDimzns,1);

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
void linearSystem::allocD1Matrix(){
  ACE_CHECK_IERROR(mD1zs == 0,"linearSystem::allocD1Matrix");
  if (!mDimLambda)
    return;
  if (mDimx)
    mD1x=new aceMatrix(mDimLambda,mDimx);
  if (mDimzs-1>0)
    mD1zs=new aceMatrix(mDimLambda,mDimzs-1);
  if (mDimzns)
    mD1zns=new aceMatrix(mDimLambda,mDimzns);
  mD1l = new aceMatrix(mDimLambda,mDimLambda);
  mD1s = new aceMatrix(mDimLambda,1);
  if (mDimx)
    mR = new aceMatrix(mDimx,mDimLambda);

}
void linearSystem::set2matrix(){

  mA2x=mA1x;
  mA2zs=mA1zs;
  
  mB2x=mB1x;
  mB2zs=mB1zs;
  mB2l=new aceMatrix(mNbNonDynEquations,mDimLambda);

  mD2x=mD1x;
  mD2zs=mD1zs;
  mD2l=mD1l;
  mA2s=mA1s;
  if (mDimx)
    mA2sti = new aceMatrix(mDimx,1);
  mB2s=mB1s;
  mD2s=mD1s;

  if (mDimzns){
    if (mDimx){
      *mR=prod(*mA1zns,*mC1l);
      *mA2x = *mA1x + prod(*mA1zns,*mC1x);
      if (mDimzs)
	*mA2zs = *mA1zs + prod(*mA1zns,*mC1zs);
    }
    if (mNbNonDynEquations){
      if (mDimx)
	*mB2x = *mB1x + prod(*mB1zns,*mC1x);
      if (mDimzs)
	*mB2zs = *mB1zs + prod(*mB1zns,*mC1zs);
      if (mDimx)
	*mD2x=*mD1x + prod(*mD1zns,*mC1x);
      *mD2zs=*mD1zs+prod(*mD1zns,*mC1zs);
      *mD2l=*mD1l + prod(*mD1zns,*mC1l);

      *mB2l=prod(*mB1zns,*mC1l);
    }
  }

}
void linearSystem::freeD1Matrix(){
  if (mD1x)
    delete mD1x;
  if (mD1zs)
    delete mD1zs;
  if (mD1zns)
    delete mD1zns;
  if (mD1l)
    delete mD1l;
  if (mD1s)
    delete mD1s;
  if (mR)
    delete mR;
  if (mB2l)
    delete mB2l;
  if (mA2sti)
    delete mA2sti;
}


void linearSystem::allocDiscretisation(){
  if(mDimx){
    mxti = new aceMatrix(mDimx,1);
    mxfree = new aceMatrix(mDimx,1);
    mW = new aceMatrix(mDimx,mDimx);
  }
  
  mzsti = new aceMatrix(mDimzs-1,1);
  if(mDimzns)
    mznsti = new aceMatrix(mDimzns,1);
  if (mB2zs)
    mB3zs = new aceMatrix(mNbNonDynEquations,mDimzs-1);
  if (mNbNonDynEquations && mDimLambda)
    mB3l=new aceMatrix(mNbNonDynEquations,mDimLambda);
  if (mNbNonDynEquations)
    mQfree = new aceMatrix(mNbNonDynEquations,1);

  if (mDimLambda)
    mD3zs = new aceMatrix(mDimLambda,mDimzs-1);
  if (mDimLambda)
    mD3l = new aceMatrix(mDimLambda,mDimLambda);
  if (mDimLambda)
    mPfree = new aceMatrix(mDimLambda,1);
}
void linearSystem::freeDiscretisation(){
  if (mxti)
    delete mxti;
  if (mxfree)
    delete mxfree;
  if (mzsti)
    delete mzsti;
  if (mznsti)
    delete mznsti;
  if(mW)
    delete mW;
  if (mD3l)
    delete mD3l;
  if (mD3zs)
    delete mD3zs;
  if (mB3l)
    delete mB3l;
  if (mB3zs)
    delete mB3zs;
  if (mPfree)
    delete mPfree;
  if (mQfree)
    delete mQfree;
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
  freeD1Matrix();
  freeDiscretisation();
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////DISCRETISATION
void linearSystem::readInitialValue(){
  int i;
  double aux;
  int nbGuess;
  unsigned long guess;

  cout << "read init value\n";
  try{
    ifstream pin(mFile);
    ACE_CHECK_IERROR(pin,"linearSystem::readInitialValue need init file");

    pin >> mStepNumber ;
    pin >> mH;
    pin >> nbGuess;
    for (i=0;i<nbGuess;i++){
      pin >> guess;
      mMLCP->addGuess(guess);
    }
      
    initSimulation(PARSER_TSTEP,mH);
    initSimulation(PARSER_TSTOP,mH*mStepNumber);
    cout << "mStepNumber\t"<<mStepNumber<<"\tmH\t"<<mH<<endl;
    cout << "x\n";
    for (i=0;i<mDimx;i++){
      pin >> aux;
      mxti->setValue(i,0,aux);
      cout << aux<<"\t";
    }
    
    cout << "\nzs\n";
    for (i=0;i<mDimzs-1;i++){
      pin >> aux ;
      mzsti->setValue(i,0,aux);
      cout << aux<<"\t";
    }
    cout <<"\n";
  }
  catch(...)
    {
      std::cout << "Exception caught." << endl;
      ACE_INTERNAL_ERROR("linearSystem::readInitialValue");
    }

}
void linearSystem::initSimu(){
  allocDiscretisation();
  mMLCP = new mlcp(mDimLambda,mNbNonDynEquations);
  readInitialValue();
  mStepCmp=0;

  try{
    //ACE_CHECK_IERROR(mDimx >0 && mDimLambda >0,"linearSystem::initSimu case no x or no lambda not yet implemented");
    if (mDimx){
      *mW = SimpleMatrix(mDimx,mDimx,IDENTITY) - mH*mTheta*(*mA2x);
      mW->PLUInverseInPlace();
      *mB3zs = *mB2zs + mH*mTheta*prod(*mB2x,prod(*mW,*mA2zs)) ;
      *mD3zs = *mD2zs+mH*mTheta*prod(*mD2x,prod(*mW,*mA2zs));
      if (mDimLambda){
	*mB3l = mH*prod(*mB2x,prod(*mW,*mR))+*mB2l;
	*mD3l = mH*prod(*mD2x,prod(*mW,*mR))+*mD2l;
      }
    }else{
      *mD3zs = *mD2zs;
      *mB3zs = *mB2zs;
      if (mDimLambda){
	*mB3l = *mB2l;
	*mD3l = *mD2l;
      }
    }
    *(mMLCP->mM11) = *mD3l;
    *(mMLCP->mM12) = *mD3zs;
    *(mMLCP->mM21) = *mB3l;
    *(mMLCP->mM22) = *mB3zs;
  }
  catch(SiconosException e)
    {
      std::cout << e.report() << endl;
      ACE_INTERNAL_ERROR("linearSystem::initSimu");
    }
  catch(...)
    {
      std::cout << "Exception caught." << endl;
      ACE_INTERNAL_ERROR("linearSystem::initSimu");
    }

  
}
void linearSystem::preparStep(){
  
  ExtractAndCompute2Sources();
  if (mDimx && mDimLambda){//both
    (*mxfree) = (*mxti) + mH*((1-mTheta)*prod(*mA2x,*mxti)+(1-mTheta)*prod(*mA2zs,*mzsti)+mThetap*(*mA2s)+(1-mThetap)*(*mA2sti));
    (*mPfree) = prod(*mD2x,prod(*mW,*mxfree)) + (*mD2s);
    (*mQfree) = prod(*mB2x,prod(*mW,*mxfree)) + (*mB2s);
  } else if(mDimx){//only x
    (*mxfree) = (*mxti) + mH*(1-mTheta)*prod(*mA2x,*mxti)+mH*(1-mTheta)*prod(*mA2zs,*mzsti)+mH*(*mA2s);
    (*mQfree) = prod(*mB2x,prod(*mW,*mxfree)) + (*mB2s);
  }else if (mDimLambda){//only lambda
    (*mPfree) = (*mD2s);
    (*mQfree) = (*mB2s);
  }else{
    (*mQfree) = (*mB2s);
  }
  *(mMLCP->mQ1)= *mPfree;
  
  *(mMLCP->mQ2)= *mQfree;
}
void linearSystem::computeZnstiFromX_Zs(){
  if (mDimLambda==0)
    return;
  if (mDimx)
    (*mznsti)=prod(*mC1x,*mxti)+prod(*mC1zs,*mzsti)+prod(*mC1l,*(mMLCP->mZ1))+(*mC1s);
  else
    (*mznsti)=prod(*mC1zs,*mzsti)+prod(*mC1l,*(mMLCP->mZ1))+(*mC1s);
}
bool linearSystem::step(){
  mStepCmp++;
  if (mStepCmp >= mStepNumber)
    return false;
    bool res = mMLCP->solve();
    if (res){
      *mzsti=*(mMLCP->mZ2);
      if (mDimx)
	*mxti=prod(*mW,*mxfree)+mH*mTheta*prod(*mW,prod(*mA2zs,*mzsti))+mH*prod(*mW,prod(*mR,*(mMLCP->mZ1)));
      computeZnstiFromX_Zs();
    }
    return res;

}
void linearSystem::stopSimu(){
  mMLCP->printGuess();
  if (mMLCP)
    delete mMLCP;
  mMLCP=0;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////UNKNOWN
unknown* linearSystem::addinx(int type, component* c)
{
  unknown* res = new unknown(type,c);
  res->mIndexInVector = mx.size();
  mx.push_back(res);
  return res;
}
unknown* linearSystem::addinZs(int type, component* c)
{
  unknown* res = new unknown(type,c);
  res->mIndexInVector = mZs.size();
  mZs.push_back(res);
  return res;
}
unknown* linearSystem::addinZns(int type, component* c)
{
  unknown* res = new unknown(type,c);
  res->mIndexInVector = mZns.size();
  mZns.push_back(res);
  return res;
}

int linearSystem::getIndexUnknown (int type,int node) {
  if (type == ACE_TYPE_V){
    return node + 2*mx.size();
  }else{
    ACE_INTERNAL_WARNING("linearSystem::getIndexUnknown not implemented");
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
equationVD* linearSystem::addVdEquation(char* name){
  equationVD* res = new equationVD(name);
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
  mDimx = mx.size();
  for (i=0; i<mDimx; i++){
    mx[i]->mDynIndex = i;
    mx[i]->mIndex = mDimx+i;    
  }
  mDimzs=mZs.size();
  for (i=0; i<mDimzs; i++){
    mZs[i]->mIndex = 2*mDimx+i;
  }
  mDimzns=mZns.size();
  for (i=0; i<mDimzns; i++){
    mZns[i]->mIndex = 2*mDimx+mDimzs+i;
  }

  allocMemory();

}





///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////DYNAMICS EQUATIONS  : x'=A1x*x + A1zs*Zs + A1zns*Zns + (not A1s because time's dependent)
//fill matrix: x'=A1x * mx + A1zs * mZs + A1zns * mZns +A1s; 
void linearSystem::computedxdt(){
  if (mx.size() == 0){
    ACE_MESSAGE("No Dynamic\n");
    return;
  }
  
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
    
    // (*mA1s)=prod(*mA,*ms);
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
  int i=0;
  for (i =idStart; i < mDimx;i++){
    coefs[i]= (*mA1x)(line,i-idStart);
  }
  idStart+=mDimx;
  coefs[idStart]=0;//Because V0
  for (i =idStart+1; i < idStart+mDimzs;i++){
    coefs[i]= (*mA1zs)(line,i-idStart-1);
  }
  idStart+=mDimzs;
  for (i =idStart; i < idStart+mDimzns;i++){
    coefs[i]= (*mA1zns)(line,i-idStart);
  }
  idStart+=mDimzns;
  ACE_CHECK_IERROR(idStart+mDimx == mNbUnknowns,"linearSystem::getlinefromdxdt idStart +mDimx == mNbUnknowns");
  coefs[idStart]=(*mA1s)(line,0);

}
//Ax'=Bx+CZs+DZns+s
void linearSystem::buildABCDs(){
  int istart=0;
  if (mDimx==0)
    return;
  ACE_CHECK_IERROR( mDimx == mNbDynEquations,"linearSystem::buildABCDs : mDimx != nbDynEquations");
  //BUILD A
  if (mA)
    extractDynBockInMat(mA,istart,istart+mDimx);

  istart+=mDimx;
  //BUILD B
  if (mB)
    extractDynBockInMat(mB,istart,istart+mDimx);
  istart+=mDimx;

  //BUILD C
  if (mC)
    extractDynBockInMat(mC,istart+1,istart+mDimzs);//because V0
  istart+=mDimzs;

  //BUILD D
  if (mD)
    extractDynBockInMat(mD,istart,istart+mDimzns);    
  istart+=mDimzns;

  //BUILD s
  ACE_CHECK_IERROR(istart == mNbUnknowns,"linearSystem::buildABCDs istart == mNbUnknowns");

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
///////////////////////////////////////////////////////////////////////////SOURCES must be call each step
void linearSystem::extractSources(){
  //BUILD s
  extractNonDynBockInMat(mB1s,mNbUnknowns,mNbUnknowns+1);
  if (mDimx){
    extractDynBockInMat(ms,mNbUnknowns,mNbUnknowns+1);
    (*mA1s)=prod(*mA,*ms);
  }
  /*  cout<<"vecteur source s:";
  if (ms)
    cout<<(*ms);
  cout<<"\nvecteur source Bs1:";
  if (mB1s)
    cout<<(*mB1s);
  */
}
void linearSystem::ExtractAndCompute2Sources(){
  if (mDimx)
    *mA2sti=*mA2s;
  extractSources();
  if (mDimx && mDimLambda){
    (*mA2s)=(*mA1s)+prod(*mA1zns,*mC1s);
    (*mB2s)=(*mB1s)+prod(*mB1zns,*mC1s);
    (*mD2s)=(*mD1s)+prod(*mD1zns,*mC1s);
  }else if (mDimx){
    (*mA2s)=(*mA1s);
    (*mB2s)=(*mB1s);    
  }else if (mDimLambda){
    (*mB2s)=(*mB1s)+prod(*mB1zns,*mC1s);
    (*mD2s)=(*mD1s)+prod(*mD1zns,*mC1s);    
  }
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////LINEAR SYSTEM : 0 = B1x*x + B1zs*Zs + B1zns*Zns + (not B1s because time dependent)
void linearSystem::buildLinearSystem(){
  
  int istart = mDimx;
  //BUILD B1x
  if (mB1x)
    extractNonDynBockInMat(mB1x,istart,istart+mDimx);
  istart+=mDimx;
  
  //BUILD B1zs
  if (mB1zs)
    extractNonDynBockInMat(mB1zs,istart+1,istart+mDimzs);
  istart+=mDimzs;

  //BUILD B1zns
  if (mB1zns)
    extractNonDynBockInMat(mB1zns,istart,istart+mDimzns);    
  istart+=mDimzns;
  ACE_CHECK_IERROR(istart == mNbUnknowns,"linearSystem::buildLinearSystem istart == mNbUnknowns");
}


//copy in m the double contains in none dynamique equation from IndexBegin to IndexEnd
void linearSystem::extractNonDynBockInMat(aceMatrix * m, int IndexBegin, int IndexEnd){
  ACE_CHECK_IERROR(m,"linearSystem::extractNonDynBockInMat m null");
  ACE_CHECK_IERROR(IndexBegin>=mDimx,"linearSystem::extractNonDynBockInMat IndexBegin");
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
 //
  //R=A1zns*C1l
  //x'=A2x*x + A2zs*Zs + R*lambda+A2s
  //0=B2x*x + B2zs*Zs + B2l*lambda + B2s
  //Y=D2x*x + D2zs*Zs + D2l*lambda + D2s
void linearSystem::printSystem2(ostream& os){
  os<<"R=A1zns*C1l\nx'=A2x*x + A2zs*Zs + R*lambda+A2s\n0=B2x*x + B2zs*Zs + B2l*lambda + B2s\nY=D2x*x + D2zs*Zs + D2l*lambda + D2s\n";
  os<<"R:\n";
  if (mR)
    os<<(*mR);
  os<<"A2x:\n";
  if (mA2x)
    os<<(*mA2x);
  os<<"A2zs:\n";
  if (mA2zs)
    os<<(*mA2zs);
  os<<"A2s:\n";
  if (mA2s)
    os<<(*mA2s);
  os<<"B2x:\n";
  if (mB2x)
    os<<(*mB2x);
  os<<"B2zs:\n";
  if (mB2zs)
    os<<(*mB2zs);
  os<<"B2l:\n";
  if (mB2l)
    os<<(*mB2l);
  os<<"B2s:\n";
  if (mB2s)
    os<<(*mB2s);
  os<<"D2x:\n";
  if (mD2x)
    os<<(*mD2x);
  os<<"D2zs:\n";
  if (mD2zs)
    os<<(*mD2zs);
  os<<"D2l:\n";
  if (mD2l)
    os<<(*mD2l);
  os<<"D2s:\n";
  if (mD2s)
    os<<(*mD2s);
}
void linearSystem::printStep(ostream& os){
  int i;
  os << "xt("<<mStepCmp*mH<<")\t";
  if (mxti)
    for (i=0;i<mDimx;i++)
      os << mxti->getValue(i,0)<<"\t";
  os << "zs("<<mStepCmp*mH<<")\t";
  if (mzsti)
    for (i=0;i<mDimzs-1;i++)
      os << mzsti->getValue(i,0)<<"\t";
  
  os << "zns("<<mStepCmp*mH<<")\t";
  if (mznsti)
    for (i=0;i<mDimzns;i++)
      os << mznsti->getValue(i,0)<<"\t";
  os<<"\n";
}

void linearSystem::printEquations(ostream& os ){
  int i,n =0;
  
  printf("--->linearSystem with %d equations whose %d dynamic equations.\n",mNbEquations,mNbDynEquations);
  printf("x\n");
  for (i=0; i<mDimx; i++)
    mx[i]->print();
  printf("\nZs\n");
  for (i=0; i<mDimzs; i++)
    mZs[i]->print();
  printf("\nZns\n");
  for (i=0; i<mDimzns; i++)
    mZns[i]->print();
  printf("\n---------------------------------------\nequation");
  for (i=0; i<mDimx; i++)
    mx[i]->printdev();
  for (i=0; i<mDimx; i++)
    mx[i]->print();
  for (i=0; i<mDimzs; i++)
    mZs[i]->print();
  for (i=0; i<mDimzns; i++)
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
void linearSystem::printC1(ostream& os ){
  os <<"Zns = C1x*x + C1s*Zs + C1l*lamdba + C1s\n";
  os <<"C1x:\n";
  if (mC1x)
    os <<(*mC1x);
  os <<"C1zs:\n";
  if (mC1zs)
    os <<(*mC1zs);
  os <<"C1l:\n";
  if (mC1l)
    os <<(*mC1l);
  os <<"C1s:\n";
  if (mC1s)
    os << (*mC1s);
 }
void linearSystem::printD1(ostream& os ){
  os <<"Y = D1x*x + D1s*Zs + D1ns*Zns + D1l*lambda +D1s\n";
  os <<"D1x:\n";
  if (mD1x)
    os <<(*mD1x);
  os <<"D1zs:\n";
  if (mD1zs)
    os <<(*mD1zs);
  os <<"D1zns:\n";
  if (mD1zns)
    os <<(*mD1zns);
  os <<"D1l:\n";
  if (mD1l)
    os <<(*mD1l);
  os <<"D1s:\n";
  if (mD1s)
    os << (*mD1s);
 }
void linearSystem::printDiscretisation(ostream& os){
  os <<"W(I - h*ThetaA2x)"<<endl;
  if (mW)
    os << mW;
  os <<"xfree"<<endl;
    
  os <<"0=Qfree + B3zs * Zs + B3l * Lambda"<<endl;
  os <<"Qfree"<<endl;
  if (mQfree)
    os << (*mQfree);
  os <<"B3zs"<<endl;
  if (mB3zs)
    os << (*mB3zs);
  os <<"B3l"<<endl;
  if (mB3l)
    os<<(*mB3l);
  os <<"Y=pfree + D3zs*Zs + D3l*Lambda"<<endl;
  os << "pfree\n";
  if (mPfree)
    os<<(*mPfree);
  os <<"D3zs"<<endl;
  if (mD3zs)
    os <<(*mD3zs);
  os <<"D3l"<<endl;
  if (mD3l)
    os <<(*mD3l);
}

void linearSystem::printSystemInTabFile(char * file){
  ofstream pout(file);
  printA1(pout);
  printB1(pout);
  printC1(pout);
  printD1(pout);
  printSystem2(pout);
  
  pout.close();
}
