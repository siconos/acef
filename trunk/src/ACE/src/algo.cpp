/************************************************************************
  			algo.cpp

**************************************************************************/
#include "algo.h"
#include "linearsystemDAE.h"
#include "componentcap.h"
#include "componentind.h"
#include "componentres.h"
#include "componentdio.h"
#include "componentvsrc.h"
#include "componentvcvs.h"
#include "componentvccs.h"
#include "componentarb.h"
#include "componentcomp.h"
#include "componentisrc.h"
#include "componentmos.h"
#include "componentbjt.h"
#include "graph.h"
#include <fstream>

linearSystem * algo::spls;

algo::algo(char * file){
  ACE_CHECK_IERROR(strlen(file) < ACE_CHAR_LENGTH && strlen(file) >3,"algo::algo: file name length!");
  if (ACE_FORMULATION_WITH_INVERSION)
    spls = new linearSystem();
  else
    spls = new linearSystemDAE();
  strncpy(mFile,file,strlen(file)-3);
  strcat(mFile,"txt");
  
  strncpy(spls->mFile,file,strlen(file)-3);
  strcat(spls->mFile,"ini");
  strncpy(mSimuFile,file,strlen(file)-3);
  strcat(mSimuFile,"sim");
  

  if (!ParserInitLibrary()){
    ACE_INTERNAL_ERROR("initParserLibrary");
  }
  ParserReadFile(file);
  ParserPrintCircuit();


  
}
algo::~algo(){
  ParserStopLibrary();

  //destroy components
  int n=0;
  int i=0;
  n = mInds.size();
  for(i=0;i<n;i++)
    delete mInds[i];
  n = mCaps.size();
  for(i=0;i<n;i++)
    delete mCaps[i];
  n = mRess.size();
  for(i=0;i<n;i++)
    delete mRess[i];
  n = mIsrcs.size();
  for(i=0;i<n;i++)
    delete mIsrcs[i];
  n = mVsrcs.size();
  for(i=0;i<n;i++)
    delete mVsrcs[i];
  n = mDios.size();
  for(i=0;i<n;i++)
    delete mDios[i];
  n = mMos.size();
  for(i=0;i<n;i++)
    delete mMos[i];
  n = mVcvs.size();
  for(i=0;i<n;i++)
    delete mVcvs[i];
  n = mArbs.size();
  for(i=0;i<n;i++)
    delete mArbs[i];
  n = mVccs.size();
  for(i=0;i<n;i++)
    delete mVccs[i];
  n= mComps.size();
  for(i=0;i<n;i++)
    delete mComps[i];
  delete spls;
  
 }

////////////////////////////////////////////////////////////////////// ALGO
void algo::perform(){
  ACE_times[ACE_TIMER_EQUATION].start();
 spls->mNbNodes = ParserGetNbElementsOfType("Node");
 spls->initKCL();
 spls->addVUnknowns();
 //BUILD x VECTOR
 //get CAPACITOR from parser
 initGraph(spls->mNbNodes,ParserGetNbElementsOfType("Capacitor"));
 ParserInitComponentList("Capacitor");
 dataCAP dCap;
 while(ParserNextComponent(&dCap)){
   componentCAP *c=new componentCAP(&dCap);
   mCaps.push_back(c);
   
   int np = dCap.nodePos;
   int nn = dCap.nodeNeg;
   unknown *uout;
   if (!spls->isUnknown(ACE_TYPE_U,c,&uout)){
     c->addTensionUnknown();
     c->addTensionEquation();
     if (np == 0){
       equationKCL *eq=spls->KCL(nn);
       ACE_CHECK_IERROR(eq->mAvailable,"algo::perform : KCL not available");
       spls->addKCLinDyn(nn);
     }else if(nn ==0){
       equationKCL *eq=spls->KCL(np);
       ACE_CHECK_IERROR(eq->mAvailable,"algo::perform : KCL not available");
       spls->addKCLinDyn(np);
     }else{
       //add only edge not connected on node 0.
       graphAddEdge(np,nn,c);
     }
   }else{
     c->mU=uout;
   }
  }
 computeGraphMST();
 
 int n1;
 int n2;
 componentCAP *c=0;
 while (nextEdgeInMST(n1,n2,&c)){
      
      if (spls->KCL(n1)->mAvailable) {
	spls->addKCLinDyn(n1);
      }
      else if (spls->KCL(n2)->mAvailable){
	spls->addKCLinDyn(n2);
      }else{
	ACE_MESSAGE("algo::perform : Add unknown current in capacitor branche!!\n");
	c->addCurrentUnknown();
	c->addCurrentEquation();
	
      }
 }

 while (nextEdgeOutMST(n1,n2,&c)){   
   c->addCurrentUnknown();
   c->addCurrentEquation();
 }
 stopGraph();

 //get INDUCTOR from parser
 ParserInitComponentList("Inductor");
 dataIND dInd;
 while(ParserNextComponent(&dInd)){
   componentIND* c = new componentIND(&dInd);
   mInds.push_back(c);
   c->addUnknowns();
   c->addEquations();
 }
//BUILD Zns component
//get INDUCTOR from parser
 ParserInitComponentList("Diode");
 dataDIO dIo;
 while(ParserNextComponent(&dIo)){
   componentDIO *c=new componentDIO(&dIo);
   mDios.push_back(c);
   c->addUnknowns();
   c->addEquations();
   
 }
 ParserInitComponentList("Mos1");
 dataMOS1 mos;
 while(ParserNextComponent(&mos)){
   componentMOS *c=new componentMOS(&mos,1);
   mMos.push_back(c);
   c->addUnknowns();
   c->addEquations();
 }
  ParserInitComponentList("BJT");
  dataBJT bjt;
  while(ParserNextComponent(&bjt)){
    componentBJT *c=new componentBJT(&bjt);
    mBjt.push_back(c);
    c->addUnknowns();
    c->addEquations();
  }

 
//get RESISTOR from parser
 ParserInitComponentList("Resistor");
 dataRES dRes;
 while(ParserNextComponent(&dRes)){
   componentRES *c=new componentRES(&dRes);
   mRess.push_back(c);   
 }
 ParserComputeSourcesValues(0.0);
//get Vsource from parser
 ParserInitComponentList("Vsource");
 dataVSRC dVsrc;
 while(ParserNextComponent(&dVsrc)){
   componentVSRC *c=new componentVSRC(&dVsrc);
   c->addUnknowns();
   c->addEquations();
   mVsrcs.push_back(c);   
 }

 ParserInitComponentList("VCVS");
 dataVCVS dVcvs;
 while(ParserNextComponent(&dVcvs)){
   componentVCVS *c=new componentVCVS(&dVcvs);
   c->addUnknowns();
   c->addEquations();
   mVcvs.push_back(c);   
 }
 
//  initComponentList("ASRC");
//  dataARB dArb;
//  dataCOMP dCOMP;
//  while(nextComponent(&dArb)){
//    ACE_MESSAGE("replace arb by comp\n");
//    dCOMP.nodePos=7;
//    dCOMP.nodeNeg=12;
//    dCOMP.nodeOut=3;
//    dCOMP.vmax = 3;
//    dCOMP.vmin = 0;
//    dCOMP.vepsilon=0.1;
//    dCOMP.name=ACE_name;
//    componentCOMP *c=new componentCOMP(&dCOMP);
//    c->addUnknowns();
//    c->addEquations();
//    mArbs.push_back(c);   
//  }

 ParserInitComponentList("Comparator");
 dataCOMP dCOMP;
 while(ParserNextComponent(&dCOMP)){
   componentCOMP *c=new componentCOMP(&dCOMP);
   c->addUnknowns();
   c->addEquations();
   mComps.push_back(c);   
 }
 
 
 ParserInitComponentList("VCCS");
 dataVCCS dVCCS;
 while(ParserNextComponent(&dVCCS)){
   componentVCCS *c=new componentVCCS(&dVCCS);
   c->addUnknowns();
   c->addEquations();
   mVccs.push_back(c);   
 }
 
//get Isource from parser
 ParserInitComponentList("Isource");
 dataISRC dIsrc;
 while(ParserNextComponent(&dIsrc)){
   componentISRC *c=new componentISRC(&dIsrc);
   mIsrcs.push_back(c);   
 }
 printComponents();
 spls->preparForStamp();
 stamp();
 //  spls->printEquations();
 
 //compute matrix: x'=A1x * mx + A1zs * mZs + A1zns * mZns;
 ACE_times[ACE_TIMER_TEST_1].start();
 spls->computedxdt();
 ACE_times[ACE_TIMER_TEST_1].stop();
 stampAfterInvertion();
 ACE_MESSAGE("final equation ;\n");
 spls->printEquations();
 ACE_times[ACE_TIMER_TEST_2].start();
 spls->buildLinearSystem();
 ACE_times[ACE_TIMER_TEST_2].stop();
 spls->printA1();
 spls->printB1();
 spls->printC1();
 spls->printD1();
 spls->set2matrix();
 spls->printSystem2();
 ACE_times[ACE_TIMER_EQUATION].stop();
}
////////////////////////////////////////////////////////////////////// STAMP
//with x'=A1x * mx + A1zs * mZs + A1zns * mZns, compute curent in all capacitor branche, and fill KCL law
void algo::stampAfterInvertion(){
  int n = mCaps.size();
  for(int i=0;i<n;i++){
    ((componentCAP *)(mCaps[i]))->stampAfterInvertion();
  }
}
//write Ax'=Bx+CZs+DZns+s
void algo::stamp(){
  int n=0;
  int i=0;
  n = mInds.size();
  for(i=0;i<n;i++)
    mInds[i]->stamp();
  n = mCaps.size();
  for(i=0;i<n;i++)
    ((componentCAP *)(mCaps[i]))->stampBeforeInvertion();
  n = mRess.size();
  for(i=0;i<n;i++)
    mRess[i]->stamp();
  n = mVsrcs.size();
  for(i=0;i<n;i++)
    mVsrcs[i]->stamp();

  n = mIsrcs.size();
  for(i=0;i<n;i++)
    mIsrcs[i]->stamp();
  n = mDios.size();
  for(i=0;i<n;i++)
    mDios[i]->stamp();
  n = mMos.size();
  for(i=0;i<n;i++)
    mMos[i]->stamp();
  n = mVcvs.size();
  for(i=0;i<n;i++)
    mVcvs[i]->stamp();
  n = mArbs.size();
  for(i=0;i<n;i++)
    mArbs[i]->stamp();
  n = mComps.size();
  for(i=0;i<n;i++)
    mComps[i]->stamp();
  n = mVccs.size();
  for(i=0;i<n;i++)
    mVccs[i]->stamp();
  n = mBjt.size();
  for(i=0;i<n;i++)
    mBjt[i]->stamp();

}
////////////////////////////////////////////////////////////////////// SIMULATION
void algo::preparStep(double time){
  int n;
  int i;
  //ACE_DOUBLE time = spls->mStepCmp*spls->mH;
  ParserComputeSourcesValues(time);
  n = mVsrcs.size();
  for(i=0;i<n;i++)
    mVsrcs[i]->stampTimer();
  n = mIsrcs.size();
  for(i=0;i<n;i++)
    mIsrcs[i]->stampTimer();
}
void algo::simulate(){
  spls->initSimu();
  mSimuStream = new ofstream(mSimuFile);
  ACE_times[ACE_TIMER_SIMULATION].start();
  preparStep(spls->mStepCmp*spls->mH);
  spls->ExtractAndCompute2Sources();
  spls->preparStep();
  while(spls->step()){
    if (ACE_MUET_LEVEL != ACE_MUET)
      spls->printStep();
    spls->printStep(*mSimuStream);
    preparStep(spls->mStepCmp*spls->mH);
    spls->ExtractAndCompute2Sources();
    spls->preparStep();
  }
  ACE_times[ACE_TIMER_SIMULATION].stop();
  spls->stopSimu();
  mSimuStream->close();
  delete mSimuStream;
  mSimuStream=0;

}
////////////////////////////////////////////////////////////////////// PRINT
void algo::printComponents(){
  int n =0;
  int i;
  ACE_MESSAGE("ACE components:\n");
  n=mIsrcs.size();
  printf("-->%d I source:\n",n);
  for(i=0;i<n;i++)
    mIsrcs[i]->print();
  
  n=mVsrcs.size();
  printf("-->%d Vsource:\n",n);
  for(i=0;i<n;i++)
    mVsrcs[i]->print();
  
  n=mArbs.size();
  printf("-->%d Arb:\n",n);
  for(i=0;i<n;i++)
    mArbs[i]->print();

  n=mVccs.size();
  printf("-->%d VCCS:\n",n);
  for(i=0;i<n;i++)
    mVccs[i]->print();

  n=mVcvs.size();
  printf("-->%d VCVS:\n",n);
  for(i=0;i<n;i++)
    mVcvs[i]->print();
  

  n=mInds.size();
  printf("-->%d Inductors:\n",n);
  for(i=0;i<n;i++)
    mInds[i]->print();
  
  n=mCaps.size();
  printf("-->%d Capacitors:\n",n);
  for(i=0;i<n;i++)
    mCaps[i]->print();
  
  n=mRess.size();
  printf("-->%d Resistors:\n",n);
  for(i=0;i<n;i++)
    mRess[i]->print();
  
  n=mDios.size();
  printf("-->%d Diode:\n",n);
  for(i=0;i<n;i++)
    mDios[i]->print();

  n=mMos.size();
  printf("-->%d MOS:\n",n);
  for(i=0;i<n;i++)
    mMos[i]->print();

  n=mBjt.size();
  printf("-->%d BJT:\n",n);
  for(i=0;i<n;i++)
    mBjt[i]->print();

  n=mComps.size();
  printf("-->%d COMPARATOR:\n",n);
  for(i=0;i<n;i++)
    mComps[i]->print();
  

  
}
