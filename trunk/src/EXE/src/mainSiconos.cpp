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

/*! \file mainSiconos.cpp

*/
#include "ACEF.h"
#include "linearsystem.h"

#include "SiconosKernel.hpp"
#include <stdio.h>
#include <stdlib.h>

static  algo *sAlgo=0;
/************************************************************/
/************************************************************/
/************************************************************/
/*call back for the source*/
/*call back for the formulation with inversion*/
void (bLDS) (double t, unsigned int N, double* b, unsigned int z, double*zz){
  sAlgo->preparStep(t);
  sAlgo->spls->extractDynamicSystemSource();
  sAlgo->spls->updateDynamicSystemSource();
  for (unsigned int i=0;i<N;i++){
    b[i]= sAlgo->spls->mA2s->getValue(i);
  }
}

void (DAEbLDS) (double t, unsigned int N, double* b, unsigned int z, double*zz){
  sAlgo->preparStep(t);
  sAlgo->spls->extractDynamicSystemSource();
  sAlgo->spls->updateDynamicSystemSource();
  sAlgo->spls->extractInteractionSource();
  sAlgo->spls->updateInteractionSource();
  unsigned int s = sAlgo->spls->mB2s->dimRow;
  for (unsigned int i=0;i<N-s;i++){
    b[i]= sAlgo->spls->mA2s->getValue(i);
  }
  for (unsigned int i=N-s;i<N;i++){
    b[i]= sAlgo->spls->mB2s->getValue(i-(N-s));
  }
}
/*call back for the formulation without inversion*/
void (eLDS) (double t, unsigned int N, double* e, unsigned int z, double*zz){
  sAlgo->preparStep(t);
  sAlgo->spls->extractInteractionSource();
  sAlgo->spls->updateInteractionSource();
  int s = sAlgo->spls->mB2s->dimRow;
  int m = sAlgo->spls->mD2s->dimRow;
  for (int i=0;i<s ;i++){
    e[i]= sAlgo->spls->mB2s->getValue(i);
  }
  for (int i=0;i<m ;i++){
    e[i+s]= sAlgo->spls->mD2s->getValue(i);
  }
}
void (DAEeLDS) (double t, unsigned int N, double* e, unsigned int z, double*zz){
  sAlgo->preparStep(t);
  sAlgo->spls->extractInteractionSource();
  sAlgo->spls->updateInteractionSource();
  int m = sAlgo->spls->mD2s->dimRow;
  for (int i=0;i<m ;i++){
    e[i]= sAlgo->spls->mD2s->getValue(i);
  }
}
/************************************************************/
/************************************************************/
/************************************************************/
/************************************************************/
/*main program*/
int main(int argc, char **argv){
  if (argc<3){
    printf("usage : noselect file.cir ENUM|SIMPLEX|PATH|DIRECT_ENUM|DIRECT_SIMPLEX|DIRECT_PATH [INV/NOINV/MNA]\n");
    return 0;
  }
  ACE_FORMULATION=ACE_FORMULATION_SEMI_EXPLICT;
  ACE_USE_NL_MOS = 0;

  /*MOS Parameters*/
  /*formulation with inversion ?*/
  if (argc>3){
    if (!strcmp(argv[3],"NOINV"))
      ACE_FORMULATION=ACE_FORMULATION_WITHOUT_INVERT;
    if (!strcmp(argv[3],"MNA"))
      ACE_FORMULATION=ACE_FORMULATION_MNA;
  }


  string solverDirEnum = "DIRECT_ENUM" ;
  string solverDirPath = "DIRECT_PATH" ;
  string solverDirSimplex = "DIRECT_SIMPLEX" ;
  string solverEnum = "ENUM" ;
  string solverSimplex = "SIMPLEX" ;
  string solverPath = "PATH" ;
  string * solverName =0;
  // One Step non smooth problem
  IntParameters iparam(10);
  DoubleParameters dparam(10);
  
  double* floatWorkingMem = 0;
  int * intWorkingMem= 0;

  int freq = 1000;
  int Nfreq=0;
  int cmp=0;

  int NbDataMax = 10000;
  int NData=0;
  int len = strlen(argv[1]);
  char *simTraj;
  simTraj =(char*)malloc(sizeof(char)*len);
  strncpy(simTraj,argv[1],len-3);
  strcat(simTraj,"sim");
  cout<<"simTraj : "<<simTraj<<endl;
  /*Get the solver type*/
  if (!strcmp(argv[2],"ENUM")){
    ACE_SOLVER_TYPE = ACE_SOLVER_ENUM;
  }else if(!strcmp(argv[2],"SIMPLEX")){
    ACE_SOLVER_TYPE = ACE_SOLVER_SIMPLEX;
  }else if(!strcmp(argv[2],"PATH")){
    ACE_SOLVER_TYPE = ACE_SOLVER_PATH;
  }else if(!strcmp(argv[2],"DIRECT_ENUM")){
    ACE_SOLVER_TYPE=ACE_SOLVER_NUMERICS_DIRECT_ENUM;
  }else if(!strcmp(argv[2],"DIRECT_SIMPLEX")){
    ACE_SOLVER_TYPE=ACE_SOLVER_NUMERICS_DIRECT_SIMPLEX;
  }else if(!strcmp(argv[2],"DIRECT_PATH")){
    ACE_SOLVER_TYPE=ACE_SOLVER_NUMERICS_DIRECT_PATH;
  }else{
    printf("param2 must be : ENUM|SIMPLEX|PATH|DIRECT_ENUM|DIRECT_SIMPLEX|DIRECT_PATH \n");
    return 0;
  }
//   if (!strcmp(argv[3],"0"))
//     ACE_MUET_LEVEL=0;
//   else if (!strcmp(argv[3],"1"))
//     ACE_MUET_LEVEL=1;
//   else if (!strcmp(argv[3],"2"))
//     ACE_MUET_LEVEL=2;
//   else
  ACE_MUET_LEVEL=9;
  
//   if (!strcmp(argv[4],"SPARSE"))
//     ACE_MAT_TYPE=SPARSE;
//   else
    ACE_MAT_TYPE=DENSE;

    /************************************************************/
    /************************************************************/
    /*Solver options*/
    iparam[5] = 5; // nb config
    if (ACE_SOLVER_TYPE == ACE_SOLVER_ENUM){
    iparam[0] = 0; // verbose
    dparam[0] = 0.0000001; // Tolerance
    solverName = &solverEnum;
  }else if(ACE_SOLVER_TYPE == ACE_SOLVER_SIMPLEX){
    iparam[0]= 1000000;
    iparam[1]=1;
    dparam[0]=1e-12;
    dparam[1]=1e-12;
    dparam[2]=1e-9;
    solverName = &solverSimplex;
  }else if(ACE_SOLVER_TYPE == ACE_SOLVER_PATH){
    iparam[0]=0;//verbose
    dparam[0]=1e-12;
    solverName = &solverPath;
  }else if (ACE_SOLVER_TYPE==ACE_SOLVER_NUMERICS_DIRECT_ENUM){
    iparam[0] = 0; // enum verbose
    iparam[5] = 2; // nb config
    iparam[6] = 0; // direct verbose
    dparam[0] = 0.0000001; // Tolerance
    dparam[5] = 0.0000001; // Tolerance neg
    dparam[5] = 1e-20; // Tolerance pos
    solverName = &solverDirEnum;
  }else if (ACE_SOLVER_TYPE==ACE_SOLVER_NUMERICS_DIRECT_SIMPLEX){
    //simplex option
    iparam[0]= 1000000;
    iparam[1]=0;//verbose
    dparam[0]=1e-12;
    dparam[1]=1e-12;
    dparam[2]=1e-9;
    //direct options
    iparam[6] = 0; // direct verbose
    dparam[5] = 0.0000001; // Tolerance
    dparam[5] = 1e-20; // Tolerance pos
    solverName = &solverDirSimplex;
   }else if (ACE_SOLVER_TYPE==ACE_SOLVER_NUMERICS_DIRECT_PATH){
    //path option
    iparam[0]=0;//verbose
    //direct options
    iparam[6] = 0; // direct verbose
    solverName = &solverDirPath;
    dparam[0]=1e-12;
    dparam[5]=1e-12;
    dparam[6]=1e-12;
  }

  ACE_times[ACE_TIMER_MAIN].start();
  ACE_INIT();
  /*build object from ACEF*/
  sAlgo=new algo(argv[1]);
  sAlgo->perform();
  if (!sAlgo->spls->mDimLambda){
    sAlgo->simulate();
    ACE_times[ACE_TIMER_MAIN].stop();
    ACE_PRINT_TIME();
    ACE_STOP();
    return 0;
  }
  int ACEFdimX = sAlgo->spls->mDimx;
  int dimX=ACEFdimX;
  int s=sAlgo->spls->mB2zs->getDimRow();
  int m=sAlgo->spls->mD2l->getDimRow();
  aceMatrix * DAE_M =0;
  aceMatrix * DAE_A =0;
  aceVector* DAE_X0 =0;
  aceVector* DAE_As =0;
  
  if (ACE_FORMULATION==ACE_FORMULATION_WITHOUT_INVERT){
    dimX += s;
    DAE_M = new aceMatrix(dimX,dimX);
    DAE_A = new aceMatrix(dimX,dimX);
    DAE_X0 = new aceVector(dimX);
    DAE_As = new aceVector(dimX);
  }
  
  sAlgo->spls->allocForInitialValue();
  sAlgo->spls->readInitialValue();
  double h=0.00001;
  h=sAlgo->spls->mH;
  double finalTime = h*sAlgo->spls->mStepNumber ;

  //*****BUILD THE DYNAMIC SYSTEM
  SP::FirstOrderLinearDS aDS ;
  SP::SiconosMatrix A2x =createSPtrSiconosMatrix(*(sAlgo->spls->mA2x));
  SP::SiconosVector Axti =createSPtrSiconosVector(*(sAlgo->spls->mxti));
  SP::SiconosVector A2s =createSPtrSiconosVector(*(sAlgo->spls->mA2s));
  if (ACE_FORMULATION==ACE_FORMULATION_SEMI_EXPLICT){
    //    aDS.reset(new FirstOrderLinearDS(*(sAlgo->spls->mxti),*(sAlgo->spls->mA2x),*(sAlgo->spls->mA2s)));
    aDS.reset(new FirstOrderLinearDS(Axti,A2x,A2s));
    aDS->setComputebFunction(&bLDS);
  }else if (ACE_FORMULATION==ACE_FORMULATION_WITHOUT_INVERT){
    DAE_As->setBlock(0,*(sAlgo->spls->mA2s));
    DAE_As->setBlock(ACEFdimX,*(sAlgo->spls->mB2s));
    DAE_X0->setBlock(0,*(sAlgo->spls->mxti));
    DAE_X0->setBlock(ACEFdimX,*(sAlgo->spls->mzsti));
    DAE_A->setBlock(0,0,*(sAlgo->spls->mA2x));
    DAE_A->setBlock(ACEFdimX,0,*(sAlgo->spls->mB2x));
    DAE_A->setBlock(0,ACEFdimX,*(sAlgo->spls->mA2zs));
    DAE_A->setBlock(ACEFdimX,ACEFdimX,*(sAlgo->spls->mB2zs));
    for (int ii =0; ii < ACEFdimX; ii++)
      DAE_M->setValue(ii,ii,1);
    SP::SiconosMatrix SP_DAE_A = createSPtrSiconosMatrix(*DAE_A);
    SP::SiconosVector SP_DAE_X0 =createSPtrSiconosVector(*DAE_X0);
    SP::SiconosVector SP_DAE_As =createSPtrSiconosVector(*DAE_As);
    aDS.reset(new FirstOrderLinearDS(SP_DAE_X0,SP_DAE_A,SP_DAE_As));
    SP::SiconosMatrix SP_DAE_M = createSPtrSiconosMatrix(*DAE_M);
    aDS->setMPtr(SP_DAE_M);
    aDS->setComputebFunction(&DAEbLDS);
  }else{
    ACE_INTERNAL_ERROR("main, ACE_FORMULATION not managed");
  }
  cout<<"FirstOrderLinearTIDS with :"<<endl;
  aDS->display();
  DynamicalSystemsSet  Inter_DS ;
  Inter_DS.insert(aDS);

  //******BUILD THE RELATION
  aceMatrix* C= 0;
  aceMatrix* D= 0;
  aceMatrix* B= 0;
  if (ACE_FORMULATION==ACE_FORMULATION_SEMI_EXPLICT){
    C= new aceMatrix(s+m,dimX);
    D= new aceMatrix(s+m,s+m);
    B= new aceMatrix(dimX,s+m);
    C->setBlock(0,0,*(sAlgo->spls->mB2x));
    C->setBlock(s,0,*(sAlgo->spls->mD2x));
    D->setBlock(0,0,*(sAlgo->spls->mB2zs));
    D->setBlock(0,s,*(sAlgo->spls->mB2l));
    D->setBlock(s,0,*(sAlgo->spls->mD2zs));
    D->setBlock(s,s,*(sAlgo->spls->mD2l));
    B->setBlock(0,0,*(sAlgo->spls->mA2zs));
    B->setBlock(0,s,*(sAlgo->spls->mR));
  }else  if (ACE_FORMULATION==ACE_FORMULATION_WITHOUT_INVERT){
    C= new aceMatrix(m,dimX);
    D= new aceMatrix(m,m);
    B= new aceMatrix(dimX,m);
    C->setBlock(0,0,*(sAlgo->spls->mD2x));
    C->setBlock(0,ACEFdimX,*(sAlgo->spls->mD2zs));

    D->setBlock(0,0,*(sAlgo->spls->mD2l));

    B->setBlock(0,0,*(sAlgo->spls->mR));
    B->setBlock(ACEFdimX,0,*(sAlgo->spls->mB2l));
  }
  SP::SiconosMatrix SP_C =  createSPtrSiconosMatrix(*C);
  SP::SiconosMatrix SP_B =  createSPtrSiconosMatrix(*B);
  SP::FirstOrderLinearR aR( new FirstOrderLinearR(SP_C,SP_B));
  SP::SiconosMatrix SP_D =  createSPtrSiconosMatrix(*D);
  aR->setDPtr(SP_D);
  if (ACE_FORMULATION==ACE_FORMULATION_SEMI_EXPLICT){
    aR->setComputeEFunction(&eLDS);
  }else{
    aR->setComputeEFunction(&DAEeLDS);
  }
  cout<<"FirstOrderLinearR with :"<<endl;
  aR->display();

  //*****BUILD THE NSLAW
  SP::NonSmoothLaw aNSL;
  int NSLawSize=0;
  if (ACE_FORMULATION==ACE_FORMULATION_SEMI_EXPLICT){
    NSLawSize=m+s;
    aNSL.reset(new MixedComplementarityConditionNSL(m,s));
  }else{
    NSLawSize=m;
    aNSL.reset(new ComplementarityConditionNSL(m));
  }

  //****BUILD THE INTERACTION
  SP::Interaction aI( new Interaction("MLCP",Inter_DS,1,NSLawSize,aNSL,aR));
  //****BUILD THE SYSTEM
  SP::NonSmoothDynamicalSystem  aNSDS( new NonSmoothDynamicalSystem(aDS,aI));
  
  SP::Model  aM( new Model(0,finalTime));
  aM->setNonSmoothDynamicalSystemPtr(aNSDS);
  SP::TimeDiscretisation  aTD( new TimeDiscretisation(0,h));
  SP::TimeStepping aS ( new TimeStepping(aTD));
  
  //*****BUILD THE STEP INTEGRATOR
  SP::OneStepIntegrator  aMoreau ;
  if (ACE_FORMULATION==ACE_FORMULATION_SEMI_EXPLICT){
    aMoreau.reset( new Moreau(aDS,0.5));
  }else{
    aMoreau.reset( new Moreau2(aDS,0.5));
  }
  aS->recordIntegrator(aMoreau);
  SP::NonSmoothSolver  mySolver( new NonSmoothSolver((*solverName),iparam,dparam,floatWorkingMem,intWorkingMem));

  //**** BUILD THE STEP NS PROBLEM
  SP::MLCP  aMLCP ;
  if (ACE_FORMULATION==ACE_FORMULATION_SEMI_EXPLICT){
    aMLCP.reset(new MLCP(mySolver,"MLCP"));
  }else{
    //    aMLCP.reset(new MLCP2(mySolver,"MLCP2"));
  }
  aS->recordNonSmoothProblem(aMLCP);
  aM->initialize(aS);
  //  Alloc working mem for the solver
  if (ACE_SOLVER_TYPE==ACE_SOLVER_NUMERICS_DIRECT_ENUM ||
      ACE_SOLVER_TYPE==ACE_SOLVER_ENUM  ||
      ACE_SOLVER_TYPE==ACE_SOLVER_NUMERICS_DIRECT_SIMPLEX ||
      ACE_SOLVER_TYPE==ACE_SOLVER_NUMERICS_DIRECT_PATH){
    int aux =mlcp_driver_get_iwork(aMLCP->getNumericsMLCP().get(),mySolver->numericsSolverOptions().get());
    intWorkingMem = (int*)malloc(aux*sizeof(int));
    mySolver->numericsSolverOptions()->iWork = intWorkingMem;
    aux =mlcp_driver_get_dwork(aMLCP->getNumericsMLCP().get(),mySolver->numericsSolverOptions().get());
    floatWorkingMem = (double*)malloc(aux*sizeof(double));
    mySolver->numericsSolverOptions()->dWork = floatWorkingMem;
  }
  /*init solver*/
  mlcp_driver_init(aMLCP->getNumericsMLCP().get(),mySolver->numericsSolverOptions().get());

  
  SP::SiconosVector  x = aDS->x();
  SP::SiconosVector  y = aI->y(0);
  SP::SiconosVector  lambda = aI->lambda(0);


  unsigned int count = 0; // events counter. 
  // do simulation while events remains in the "future events" list of events manager. 
  int N = sAlgo->spls->mStepNumber;
  cout << " ==== Start of  simulation : "<<N<<" steps====" << endl;
  Nfreq = N/freq;
  if (N > NbDataMax){
    NData = N/NbDataMax;
  }
  
  ofstream pout(simTraj);

  dataPrint * pPrint;
  ParserInitPrintElem();
  pout<<"Index\ttime";
  while(ParserGetPrintElem((void**)&pPrint)){
    pout<<"\t\t";
    pout<<pPrint->name;
  }
  pout<<endl<<endl;
  /*print t0*/
  ParserInitPrintElem();
  pout <<0<<"\t"<<0;
  while(ParserGetPrintElem((void**)&pPrint)){
    pout<<"\t\t";
    double aux = sAlgo->spls->mzsti->getValue(pPrint->node1-1);
    if (pPrint->node2 >0)
      aux -=  sAlgo->spls->mzsti->getValue(pPrint->node2-1);
    pout << aux;
  }
  pout<<endl;
 
  //  unsigned int vin = ParserGetId("Vsource","vin");
  //  cout <<"vin"<<vin<<endl;
  //  ParserSetIputFromId("Vsource",vin,10.0);

  
  for(int k = 0 ; k < N ; k++){
    //  double new_v=10*cos(((double)k)*6.28/((double)N));
    //ParserSetIputFromId("Vsource",vin,new_v);
    if ((!Nfreq) || k % Nfreq == 0){
      cout<<"..."<<cmp<<endl;
      cmp++;
    }
    // solve ... 
    aS->computeOneStep();

    //    if ((!NData) || k % NData == 0){
      ParserInitPrintElem();
      //pout <<k+1<<"\t"<<(k+1)*h;
      pout <<(k+1)*h;
      while(ParserGetPrintElem((void**)&pPrint)){
	pout<<"\t\t";
	if (ACE_FORMULATION==ACE_FORMULATION_SEMI_EXPLICT){
	  double aux = (*lambda)(pPrint->node1-1);
	  if (pPrint->node2 >0)
	    aux -= (*lambda)(pPrint->node2-1);
	  pout << aux;
	}else{
	  double aux = (*x)(ACEFdimX+pPrint->node1-1);
	  if (pPrint->node2 >0)
	    aux -= (*x)(ACEFdimX+pPrint->node2-1);
	  pout << aux;
	}
      }
      //      pout<<"\t"<<(*x)(4)<<"\t";
      pout<<endl;
      //    }

    aS->nextStep();
	
  }
  aMLCP->reset();
  pout.close();
  cout << "===== End of simulation. ==== " << endl; 



  
  ACE_times[ACE_TIMER_MAIN].stop();
  ACE_PRINT_TIME();
  if (floatWorkingMem)
    free(floatWorkingMem);
  if (intWorkingMem)
    free(intWorkingMem);
  ACE_STOP();

  return 0;
 
}
