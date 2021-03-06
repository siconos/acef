/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2019 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/

/*! \file algo.cpp

*/
/************************************************************************
  			algo.cpp

**************************************************************************/
#include "algo.h"
#include "linearsystemDAE.h"
#include "linearsystemMNA.h"
#include "linearsystemMNA_V.h"
#include "linearsystemSTAMP_ONLY.h"
#include "componentcap.h"
#include "componentcapdae.h"
#include "componentcapmna.h"
#include "componentcapmna_v.h"
#include "componentind.h"
#include "componentres.h"
#include "componentdio.h"
#include "componentvsrc.h"
#include "componentvcvs.h"
#include "componentvccs.h"
#include "componentarb.h"
#include "componentcomp.h"
#include "componentrelay.h"
#include "componentisrc.h"
#include "componentmos.h"
#include "componentmos_nl.h"
#include "componentmos_nl2.h"
#include "componentbjt.h"
#include "graph.h"
#include "node.h"
#include "componentcomp_SE1.h"
#include "componentrelay_SE1.h"
#include "componentmos_SE1.h"
#include "componentind_SE1.h"
#include "componentres_SE1.h"
#include "componentdio_SE1.h"
#include "componentvsrc_SE1.h"
#include "componentvcvs_SE1.h"
#include <fstream>

linearSystem * algo::spls;
algo * algo::sAlgo;

algo::algo(char * file){
  ACE_CHECK_IERROR(strlen(file) < ACE_CHAR_LENGTH && strlen(file) >3,"algo::algo: file name length!");
  if (ACE_FORMULATION == ACE_FORMULATION_SEMI_EXPLICT || ACE_FORMULATION == ACE_FORMULATION_SE1)
    spls = new linearSystem();
  else if (ACE_FORMULATION == ACE_FORMULATION_WITHOUT_INVERT)
    spls = new linearSystemDAE();
  else if (ACE_FORMULATION == ACE_FORMULATION_MNA)
    spls = new linearSystemMNA();
  else if (ACE_FORMULATION == ACE_FORMULATION_MNA_V)
    spls = new linearSystemMNA_V();
  else if (ACE_FORMULATION == ACE_FORMULATION_STAMP_ONLY)
    spls = new linearSystemSTAMP_ONLY();
  else
    ACE_INTERNAL_ERROR("algo::algo, ACE_FORMULATION wrong type");
    
    
  
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
  n = mMos_NL.size();
  for(i=0;i<n;i++)
    delete mMos_NL[i];
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
  n= mRelays.size();
  for(i=0;i<n;i++)
    delete mRelays[i];
  delete spls;
  
 }

void algo::parseComponents(){
 
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
 if (!ACE_USE_NL_MOS){
   ParserInitComponentList("Mos1");
   dataMOS1 mos;
   while(ParserNextComponent(&mos)){
     componentMOS *c=new componentMOS(&mos,ACE_MOS_NB_HYP);
     mMos.push_back(c);
     c->addUnknowns();
     c->addEquations();
   }
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
   if (dCOMP.vepsilon == 0){
     componentRELAY *c=new componentRELAY(&dCOMP);
     c->addUnknowns();
     c->addEquations();
     mRelays.push_back(c);   
   }else{
     componentCOMP *c=new componentCOMP(&dCOMP);
     c->addUnknowns();
     c->addEquations();
     mComps.push_back(c);
   }
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
 if (ACE_USE_NL_MOS && !ACE_USE_SMOOTH_MOS){
   ParserInitComponentList("Mos1");
   dataMOS1 mos;
   while(ParserNextComponent(&mos)){
     componentMOS_NL *c=new componentMOS_NL(&mos);
     mMos_NL.push_back(c);
     c->addUnknowns();
     c->addEquations();
   }
 }else if (ACE_USE_SMOOTH_MOS){
   ParserInitComponentList("Mos1");
   dataMOS1 mos;
   while(ParserNextComponent(&mos)){
     componentMOS_NL2 *c=new componentMOS_NL2(&mos);
     mMos_NL.push_back(c);
     c->addUnknowns();
     c->addEquations();
   }
 }

 
 printComponents();

}




void algo::parseComponents_SE(){
 
 //get INDUCTOR from parser
 ParserInitComponentList("Inductor");
 dataIND dInd;
 while(ParserNextComponent(&dInd)){
   componentIND_SE1* c = new componentIND_SE1(&dInd);
   mInds.push_back(c);
   c->addUnknowns();
   c->addEquations();
 }
//BUILD Zns component
//get INDUCTOR from parser
 ParserInitComponentList("Diode");
 dataDIO dIo;
 while(ParserNextComponent(&dIo)){
   componentDIO_SE1 *c=new componentDIO_SE1(&dIo);
   mDios.push_back(c);
   c->addUnknowns();
   c->addEquations();
   
 }
 if (!ACE_USE_NL_MOS){
   ParserInitComponentList("Mos1");
   dataMOS1 mos;
   while(ParserNextComponent(&mos)){
     componentMOS_SE1 *c=new componentMOS_SE1(&mos,ACE_MOS_NB_HYP);
     mMos.push_back(c);
     c->addUnknowns();
     c->addEquations();
   }
 }
 ParserInitComponentList("BJT");
 dataBJT bjt;
 while(ParserNextComponent(&bjt)){
   cout<<"WARNING : not implemented for SE1"<<endl;
   exit(0);
   componentBJT *c=new componentBJT(&bjt);
   mBjt.push_back(c);
   c->addUnknowns();
   c->addEquations();
 }

 
//get RESISTOR from parser
 ParserInitComponentList("Resistor");
 dataRES dRes;
 while(ParserNextComponent(&dRes)){
   componentRES_SE1 *c=new componentRES_SE1(&dRes);
   mRess.push_back(c);   
 }
 ParserComputeSourcesValues(0.0);
//get Vsource from parser
 ParserInitComponentList("Vsource");
 dataVSRC dVsrc;
 while(ParserNextComponent(&dVsrc)){
   componentVSRC_SE1 *c=new componentVSRC_SE1(&dVsrc);
   c->addUnknowns();
   c->addEquations();
   mVsrcs.push_back(c);   
 }

 ParserInitComponentList("VCVS");
 dataVCVS dVcvs;
 while(ParserNextComponent(&dVcvs)){
   componentVCVS_SE1 *c=new componentVCVS_SE1(&dVcvs);
   c->addUnknowns();
   c->addEquations();
   mVcvs.push_back(c);   
 }
 


 ParserInitComponentList("Comparator");
 dataCOMP dCOMP;
 while(ParserNextComponent(&dCOMP)){
   if (dCOMP.vepsilon == 0){
     componentRELAY_SE1 *c=new componentRELAY_SE1(&dCOMP);
     c->addUnknowns();
     c->addEquations();
     mRelays.push_back(c);   
   }else{
     componentCOMP_SE1 *c=new componentCOMP_SE1(&dCOMP);
     c->addUnknowns();
     c->addEquations();
     mComps.push_back(c);
   }
 }
 
 
 ParserInitComponentList("VCCS");
 dataVCCS dVCCS;
 while(ParserNextComponent(&dVCCS)){
   cout<<"WARNING : not implemented for SE1"<<endl;
   exit(0);
   componentVCCS *c=new componentVCCS(&dVCCS);
   c->addUnknowns();
   c->addEquations();
   mVccs.push_back(c);   
 }
 
//get Isource from parser
 ParserInitComponentList("Isource");
 dataISRC dIsrc;
 while(ParserNextComponent(&dIsrc)){
   cout<<"WARNING : not implemented for SE1"<<endl;
   exit(0);
   componentISRC *c=new componentISRC(&dIsrc);
   mIsrcs.push_back(c);   
 }
 if (ACE_USE_NL_MOS && !ACE_USE_SMOOTH_MOS){
   ParserInitComponentList("Mos1");
   dataMOS1 mos;
   while(ParserNextComponent(&mos)){
     cout<<"WARNING : not implemented for SE1"<<endl;
     exit(0);
     componentMOS_NL *c=new componentMOS_NL(&mos);
     mMos_NL.push_back(c);
     c->addUnknowns();
     c->addEquations();
   }
 }else if (ACE_USE_SMOOTH_MOS){
   ParserInitComponentList("Mos1");
   dataMOS1 mos;
   while(ParserNextComponent(&mos)){
     cout<<"WARNING : not implemented for SE1"<<endl;
     exit(0);
     componentMOS_NL2 *c=new componentMOS_NL2(&mos);
     mMos_NL.push_back(c);
     c->addUnknowns();
     c->addEquations();
   }
 }

 
 printComponents();

}




////////////////////////////////////////////////////////////////////// ALGO
void algo::perform(){
  ACE_times[ACE_TIMER_EQUATION].start();
  spls->mNbNodes = ParserGetNbElementsOfType("Node");
  spls->initKCL();
  if (ACE_FORMULATION != ACE_FORMULATION_SE1)
    spls->addVUnknowns();

  if (ACE_FORMULATION == ACE_FORMULATION_SEMI_EXPLICT)
    performSemiExplicit();
  else if (ACE_FORMULATION == ACE_FORMULATION_WITHOUT_INVERT)
    performWithOutInvert();
  else if (ACE_FORMULATION == ACE_FORMULATION_MNA)
    performMNA();
  else if (ACE_FORMULATION == ACE_FORMULATION_MNA_V || ACE_FORMULATION == ACE_FORMULATION_STAMP_ONLY )
    performMNA_V();
  else if (ACE_FORMULATION == ACE_FORMULATION_SE1)
    performSE();
  else{
    ACE_INTERNAL_ERROR("algo::perform, ACE_FORMULATION wrong value.");
    ACE_times[ACE_TIMER_EQUATION].stop();
  }
}


void algo::performSE(){
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
     //     c->addTensionEquation();
     if (np == 0){
       equationKCL *eq=spls->KCL(nn);
       ACE_CHECK_IERROR(eq->mAvailable,"algo::perform : KCL not available");
       spls->addKCLinDyn(nn);
       spls->mNodes[nn]->setCapa(c);
     }else if(nn ==0){
       equationKCL *eq=spls->KCL(np);
       ACE_CHECK_IERROR(eq->mAvailable,"algo::perform : KCL not available");
       spls->addKCLinDyn(np);
       spls->mNodes[np]->setCapa(c);
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
 components Caux;
 while (nextEdgeInMST(n1,n2,&c)){
      
      if (spls->KCL(n1)->mAvailable) {
	spls->addKCLinDyn(n1);
	spls->mNodes[n1]->setCapa(c);
      }
      else if (spls->KCL(n2)->mAvailable){
	spls->addKCLinDyn(n2);
	spls->mNodes[n2]->setCapa(c);
      }else{
	ACE_MESSAGE("algo::perform : Add unknown current in capacitor branche!!\n");
	Caux.push_back(c);
	//	c->addCurrentUnknown();
	//	c->addCurrentEquation();
	
      }
 }

 while (nextEdgeOutMST(n1,n2,&c)){
   Caux.push_back(c);
   //   c->addCurrentUnknown();
   //   c->addCurrentEquation();
 }
 stopGraph();
 spls->addVUnknowns();
 int n = Caux.size();
 for(int i=0;i<n;i++){
   ((componentCAP *) Caux[i])->addTensionEquation();
   ((componentCAP *) Caux[i])->addCurrentUnknown();
   ((componentCAP *) Caux[i])->addCurrentEquation();
 }
 
 parseComponents_SE();
 

 spls->preparForStamp();
 stamp();

 spls->printEquations();
 
 //compute matrix: x'=A1x * mx + A1zs * mZs + A1zns * mZns;
 spls->computedxdt();
 stampAfterInvertion();
  ACE_MESSAGE("final equation ;\n");
  spls->printEquations();
 spls->buildLinearSystem();
  spls->printA1();
  spls->printB1();
  spls->printC1();
  spls->printD1();
 spls->set2matrix();
  spls->printSystem2();
}


void algo::performSemiExplicit(){
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
 parseComponents();

 spls->preparForStamp();
 stamp();

 spls->printEquations();
 //compute matrix: x'=A1x * mx + A1zs * mZs + A1zns * mZns;
 spls->computedxdt();
 stampAfterInvertion();
  ACE_MESSAGE("final equation ;\n");
  spls->printEquations();
 spls->buildLinearSystem();
  spls->printA1();
  spls->printB1();
  spls->printC1();
  spls->printD1();
 spls->set2matrix();
  spls->printSystem2();
}

void algo::performWithOutInvert(){
  //ACE_FORMULATION == ACE_FORMULATION_WITHOUT_INVERT
  ParserInitComponentList("Capacitor");
  dataCAP dCap;
  while(ParserNextComponent(&dCap)){
   componentCAPDAE *c=new componentCAPDAE(&dCap);
   mCaps.push_back(c);
   int np = dCap.nodePos;
   int nn = dCap.nodeNeg;
   //Tension unknown management
   unknown *uout;
   if (!spls->isUnknown(ACE_TYPE_U,c,&uout) || 1){
     c->addTensionUnknown();
     c->addTensionEquation();
   }else{
     c->mU=uout;
   }
   //Current unknown management
   c->addCurrentUnknown();
   c->addCurrentEquation();
  }
  parseComponents();
  spls->preparForStamp();
  stamp();
  spls->printEquations();
  
  //compute Ax'=Bx+CZs+DZns+s
  spls->buildABCDs();


  spls->buildLinearSystem();
  spls->printA1();
  spls->printB1();
  spls->printC1();
  spls->printD1();

  spls->set2matrix();
  spls->printSystem2();
  
}

void algo::performMNA(){
  //ACE_FORMULATION == ACE_FORMULATION_MNA
   //get INDUCTOR from parser
  ParserInitComponentList("Capacitor");
  dataCAP dCap;
  while(ParserNextComponent(&dCap)){
    componentCAPMNA *c=new componentCAPMNA(&dCap); 
    mCaps.push_back(c);
    c->addUnknowns();
    c->addEquations();
  }
  parseComponents();
  spls->preparForStamp();
  stamp();
  spls->printEquations();

  //compute Ax'=Bx+CZs+DZns+s
  spls->buildABCDs();
  spls->printA1();
  spls->printC1();
  spls->printD1();
  
  spls->set2matrix();
  spls->printSystem2();
 
}
void algo::performMNA_V(){
  //ACE_FORMULATION == ACE_FORMULATION_MNA_V
   //get INDUCTOR from parser
  ParserInitComponentList("Capacitor");
  dataCAP dCap;
  while(ParserNextComponent(&dCap)){
    componentCAPMNA_V *c=new componentCAPMNA_V(&dCap); 
    mCaps.push_back(c);
    c->addUnknowns();
    c->addEquations();
  }
  parseComponents();
  spls->preparForStamp();
  stamp();
  spls->printEquations();

  //compute Ax'=Bx+CZs+DZns+s
  spls->buildABCDs();
  spls->printA1();
  spls->printC1();
  spls->printD1();
  if (ACE_FORMULATION == ACE_FORMULATION_MNA_V){
    spls->set2matrix();
    spls->printSystem2();
  }
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
    ((componentCAP *)(mCaps[i]))->stampBeforeInvertion();//even if ACE_FORMULATION_WITHOUT_INVERT.
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
  n = mMos_NL.size();
  for(i=0;i<n;i++)
    mMos_NL[i]->stamp();
  n = mVcvs.size();
  for(i=0;i<n;i++)
    mVcvs[i]->stamp();
  n = mArbs.size();
  for(i=0;i<n;i++)
    mArbs[i]->stamp();
  n = mComps.size();
  for(i=0;i<n;i++)
    mComps[i]->stamp();
  n = mRelays.size();
  for(i=0;i<n;i++)
    mRelays[i]->stamp();
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
  int cmp =0;
  mSimuStream = new ofstream(mSimuFile);
  spls->mCoef=1;
  for (int i=0;i<1;i++){
    spls->initSimu(mSimuStream);
    dataPrint * pPrint;
    ParserInitPrintElem();
    //write first line.
    //   (*mSimuStream)<<"Index\ttime";
    //   while(ParserGetPrintElem((void**)&pPrint)){
    //     (*mSimuStream)<<"\t\t";
    //     (*mSimuStream)<<pPrint->name;
    //   }
    //   (*mSimuStream)<<endl<<endl;
    
    
    preparStep(spls->getCurrentTime());
    spls->ExtractAndCompute2Sources();
    spls->preparStep();
    ACE_times[ACE_TIMER_SIMULATION].start();
    while(spls->step()){
      //    if (ACE_MUET_LEVEL != ACE_MUET)
      //      spls->printStep();
      //    (*mSimuStream)<<cmp<<"\t";
      //    cmp++;
      //    spls->printStep(*mSimuStream);
      preparStep(spls->getCurrentTime());
      spls->ExtractAndCompute2Sources();
      spls->preparStep();
    }
    spls->printLog();
    //    spls->mCoef=spls->mCoef*2;
    ACE_RTOL_LOCAL=ACE_RTOL_LOCAL/2.0;
    ACE_ATOL_LOCAL=ACE_ATOL_LOCAL/2.0;
  }
  
  ACE_times[ACE_TIMER_SIMULATION].stop();
  spls->stopSimu();
    
  mSimuStream->close();
  delete mSimuStream;
  mSimuStream=0;

}
////////////////////////////////////////////////////////////////////// PRINT
void algo::printComponents(){
  if (ACE_MUET_LEVEL == ACE_MUET)
    return;

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

  n=mMos_NL.size();
  printf("-->%d MOS:\n",n);
  for(i=0;i<n;i++)
    mMos_NL[i]->print();

  n=mBjt.size();
  printf("-->%d BJT:\n",n);
  for(i=0;i<n;i++)
    mBjt[i]->print();

  n=mComps.size();
  printf("-->%d COMPARATOR:\n",n);
  for(i=0;i<n;i++)
    mComps[i]->print();

  n=mRelays.size();
  printf("-->%d RELAY:\n",n);
  for(i=0;i<n;i++)
    mRelays[i]->print();
  

  
}
void algo::computeNonLinearEquations(SiconosVector& SICONOS_X,SiconosVector& SICONOS_Lambda,SiconosVector& SICONOS_H){
  int n=mMos_NL.size();
  for(int i=0;i<n;i++)
    mMos_NL[i]->computeNL(SICONOS_X,SICONOS_Lambda,SICONOS_H);

}
void algo::computeNonLinearJacL_H(SiconosVector& SICONOS_X,SiconosVector& SICONOS_Lambda,SiconosMatrix& SICONOS_D){
  int n=mMos_NL.size();
  for(int i=0;i<n;i++)
    mMos_NL[i]->computeJacNL(SICONOS_X,SICONOS_Lambda,SICONOS_D);
}
