/************************************************************************
  			algo.cpp

**************************************************************************/
#include "algo.h"
#include "componentcap.h"
#include "componentind.h"
#include "componentres.h"
#include "componentdio.h"
#include "componentvsrc.h"
#include "componentisrc.h"
#include "graph.h"

linearSystem algo::sls;

algo::algo(){;}
algo::~algo(){

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

}


void algo::perform(){
 sls.mNbNodes = getNbElementsOfType("Node");
 sls.initKCL();
 sls.addVUnknowns();
 //BUILD x VECTOR
 //get CAPACITOR from parser
 initGraph(sls.mNbNodes,getNbElementsOfType("Capacitor"));
 initComponentList("Capacitor");
 dataCAP dCap;
 while(nextComponent(&dCap)){
   componentCAP *c=new componentCAP(&dCap);
   mCaps.push_back(c);
   
   int np = dCap.nodePos;
   int nn = dCap.nodeNeg;
   unknown *uout;
   if (!sls.isUnknown(ACE_TYPE_U,c,&uout)){
     c->addTensionUnknown();
     c->addTensionEquation();
     if (np == 0){
       equationKCL *eq=sls.KCL(nn);
       if(!eq->mAvailable)
	 ACE_INTERNAL_ERROR("algo::perform : KCL not available");
       sls.addKCLinDyn(nn);
     }else if(nn ==0){
       equationKCL *eq=sls.KCL(np);
       if(!eq->mAvailable)
	 ACE_INTERNAL_ERROR("algo::perform : KCL not available");
       sls.addKCLinDyn(np);
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
      
      if (sls.KCL(n1)->mAvailable) {
	sls.addKCLinDyn(n1);
      }
      else if (sls.KCL(n2)->mAvailable){
	sls.addKCLinDyn(n2);
      }else{
	ACE_WARNING("algo::perform : Add unknown I could be avoid!!");
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
 initComponentList("Inductor");
 dataIND dInd;
 while(nextComponent(&dInd)){
   componentIND* c = new componentIND(&dInd);
   mInds.push_back(c);
   c->addUnknowns();
   c->addEquations();
 }
//BUILD Zns component
//get INDUCTOR from parser
 initComponentList("Diode");
 dataDIO dIo;
 while(nextComponent(&dIo)){
   componentDIO *c=new componentDIO(&dIo);
   mDios.push_back(c);
   c->addUnknowns();
   c->addEquations();
   
 }

//get RESISTOR from parser
 initComponentList("Resistor");
 dataRES dRes;
 while(nextComponent(&dRes)){
   componentRES *c=new componentRES(&dRes);
   mRess.push_back(c);   
 }
//get Vsource from parser
 initComponentList("Vsource");
 dataVSRC dVsrc;
 while(nextComponent(&dVsrc)){
   componentVSRC *c=new componentVSRC(&dVsrc);
   c->addUnknowns();
   c->addEquations();
   mVsrcs.push_back(c);   
 }
//get Isource from parser
 initComponentList("Isource");
 dataISRC dIsrc;
 while(nextComponent(&dIsrc)){
   componentISRC *c=new componentISRC(&dIsrc);
   mIsrcs.push_back(c);   
 }
 sls.preparForStamp();
 stamp();
 sls.print();
 
}
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
  n = mIsrcs.size();
  for(i=0;i<n;i++)
    mIsrcs[i]->stamp();
  n = mVsrcs.size();
  for(i=0;i<n;i++)
    mVsrcs[i]->stamp();
  n = mDios.size();
  for(i=0;i<n;i++)
    mDios[i]->stamp();

}

void algo::printComponents(){
  int n =0;
  int i;
  
  n=mInds.size();
  printf("ACE %d Inductors:\n",n);
  for(i=0;i<n;i++)
    mInds[i]->print();
  
  n=mCaps.size();
  printf("ACE %d Capacitors:\n",n);
  for(i=0;i<n;i++)
    mCaps[i]->print();
  
  n=mRess.size();
  printf("ACE %d Resistors:\n",n);
  for(i=0;i<n;i++)
    mRess[i]->print();
  
  n=mIsrcs.size();
  printf("ACE %d Isource:\n",n);
  for(i=0;i<n;i++)
    mIsrcs[i]->print();
  
  n=mVsrcs.size();
  printf("ACE %d Vsource:\n",n);
  for(i=0;i<n;i++)
    mVsrcs[i]->print();
  
  n=mDios.size();
  printf("ACE %d Diode:\n",n);
  for(i=0;i<n;i++)
    mDios[i]->print();

  
}
