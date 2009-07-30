#ifndef ACEFRELATION_CPP
#define ACEFRELATION_CPP

#include "acefRelation.h"
//#include "aceMatrix.h"
#include "aceMatrix.h"
#include "algo.h"
#include "linearsystem.h"


#ifdef WITH_KERNEL_RELATION
//#define SICONOS_DEBUG

static int sOnlyOnce=1;
acefRelation::acefRelation():
  FirstOrderType2R()
{
  JacH.resize(2);
  JacG.resize(2);
  JacH[0].reset(new PluggedMatrix());
  JacH[1].reset(new PluggedMatrix());
  JacG[0].reset(new PluggedMatrix());
  JacG[1].reset(new PluggedMatrix());
}

void acefRelation::initJac(SimpleMatrix* C,SimpleMatrix* D,SimpleMatrix* B){
  mB=B;
  mC=C;
  mD=D;
  
}
void acefRelation::initialize(SP::Interaction inter)
{
  FirstOrderType2R::initialize(inter);
  unsigned int sizeY = getInteractionPtr()->getSizeOfY();
  unsigned int sizeDS = getInteractionPtr()->getSizeOfDS();
  SP::SiconosVector y = getInteractionPtr()->getYPtr(0);
  SP::SiconosVector lambda = getInteractionPtr()->getLambdaPtr(0);

  
  
  workL.reset(new SimpleVector(getInteractionPtr()->getSizeOfY()));
  JacH[0]->resize(sizeY,sizeDS);
  JacH[1]->resize(sizeY,sizeY);

  JacG[0]->resize(sizeDS,sizeDS);
  JacG[1]->resize(sizeDS,sizeY);

  
  *(JacH[0])=*mC;
  *(JacH[1])=*mD;
  *(JacG[1])=*mB;

  double t0=0;
  mLastComputedSource =-1;
  for (int i=0;i <algo::spls->mzsti->size();i++)
    workL->setValue(i,algo::spls->mzsti->getValue(i));
  *lambda = *workL;

#ifdef SICONOS_DEBUG
    computeH(t0);
    cout<<"zsti\n";
    algo::spls->mzsti->display();
#endif



  
  computeG(t0);
  computeJacH(t0,0);
  computeJacH(t0,1);
  computeJacG(t0,0);
  computeJacG(t0,1);
  *data[r]=*data[g_alpha];
#ifdef SICONOS_DEBUG
  std::cout<<"data[r (g_alpha)] init\n";
  data[r]->display();
#endif


  
}



/*********************************/
/*y = h(X,lambda)                */
/*y = CX+D Lambda                */
/*and the NL relations           */
/*********************************/
void acefRelation::computeH(double t){

  *workX = *data[x];
  SP::SiconosVector lambda = getInteractionPtr()->getLambdaPtr(0);
  *workL = *lambda;
  SP::SiconosVector Heval = getInteractionPtr()->getRelationPtr()->getHalphaPtr();

#ifdef SICONOS_DEBUG
  std::cout<<"********         computeH at "<<t<<std::endl;
  //  CX+DL+NonLinear(X,L);
  if (sOnlyOnce){
    cout<<"acefRelation mC\n";
    mC->display();
    cout<<"acefRelation mD\n";
    mD->display();
    sOnlyOnce=0;
  }
  cout<<"acefRelation X\n";
  workX->display();
  cout<<"acefRelation L\n";
  workL->display();
#endif
  
  prod(*mC,*workX,*Heval);
  prod(*mD,*workL,*Heval,false);


  
  algo::sAlgo->computeNonLinearEquations(*workX,*workL,*Heval);

  if (t > mLastComputedSource){
#ifdef SICONOS_DEBUG
  std::cout<<"********         computeH(source) at "<<t<<std::endl;
#endif
    algo::sAlgo->preparStep(t);
    algo::spls->extractInteractionSource();
    algo::spls->updateInteractionSource();
    mLastComputedSource = t;
#ifdef SICONOS_DEBUG
    std::cout<<"B2s D2s:\n";
    algo::spls->mB2s->display();
    algo::spls->mD2s->display();
#endif

  }
  int s = algo::spls->mB2s->dimRow;
  int m = algo::spls->mD2s->dimRow;
  
  for (int i=0;i<s ;i++){
    Heval->setValue(i,Heval->getValue(i)+ algo::spls->mB2s->getValue(i));
  }
  for (int i=0;i<m ;i++){
    Heval->setValue(i+s,Heval->getValue(i+s)+ algo::spls->mD2s->getValue(i));
  }
  
#ifdef SICONOS_DEBUG
  std::cout<<"modif heval : \n";
  Heval->display();
#endif
  
}



void acefRelation::computeG(double t){
  SP::SiconosVector lambda = getInteractionPtr()->getLambdaPtr(0);
  *workL = *lambda;
  prod(*mB,*workL,*data[g_alpha],true);
#ifdef SICONOS_DEBUG
  std::cout<<"modif g_alpha : \n";
  data[g_alpha]->display();
#endif
  
}

  /** default function to compute jacobianH
   *  \param double : current time
   *  \param index for jacobian (0: jacobian according to x, 1 according to lambda)
   */
void acefRelation::computeJacH(double t, unsigned int index){
  /*Because jacobian according to x is null (in current version, may be not in futur)*/
  if (!index)
    return;

  SP::SiconosVector lambda = getInteractionPtr()->getLambdaPtr(0);
  *workX = *data[x];
  *workL = *lambda;
  algo::sAlgo->computeNonLinearJacL_H(*workX,*workL,*(JacH[index]));
#ifdef SICONOS_DEBUG
  std::cout<<"computeJacH "<<index <<" at " <<" "<<t<<std::endl;
  (JacH[index])->display();
#endif

}
	
  /** default function to compute jacobianG according to lambda
   *  \param double : current time
   *  \param index for jacobian: at the time only one possible jacobian => i = 0 is the default value .
   */
void acefRelation::computeJacG(double t, unsigned int index){

}


#endif
#endif
