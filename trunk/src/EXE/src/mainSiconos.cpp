#include "ACEF.h"
#include "linearsystem.h"

#include "/usr/local/include/Siconos/Kernel/MLCP.h"
#include "SiconosKernel.h"
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
  sAlgo->spls->computeDynamicSystemSource();
  for (unsigned int i=0;i<N;i++){
    b[i]= sAlgo->spls->mA2s->getValue(i);
  }
}

void (DAEbLDS) (double t, unsigned int N, double* b, unsigned int z, double*zz){
  sAlgo->preparStep(t);
  sAlgo->spls->extractDynamicSystemSource();
  sAlgo->spls->computeDynamicSystemSource();
  sAlgo->spls->extractInteractionSource();
  sAlgo->spls->computeInteractionSource();
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
  sAlgo->spls->computeInteractionSource();
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
  sAlgo->spls->computeInteractionSource();
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
    printf("usage : noselect file.cir ENUM|SIMPLEX|PATH|DIRECT_ENUM|DIRECT_SIMPLEX|DIRECT_PATH [INV/NOINV]\n");
    return 0;
  }
  /*MOS Parameters*/
  ACE_MOS_NB_HYP=3;
  ACE_MOS_POWER_SUPPLY=3.0;
  /*formulation with inversion ?*/
  if (argc>3){
    if (!strcmp(argv[3],"NOINV"))
      ACE_FORMULATION_WITH_INVERSION=0;
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
  ACE_MUET_LEVEL=10;
  
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
  
  if (!ACE_FORMULATION_WITH_INVERSION){
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
  FirstOrderLinearTIDS * aDS = 0;
  if (ACE_FORMULATION_WITH_INVERSION){
    aDS = new FirstOrderLinearTIDS(1,*(sAlgo->spls->mxti),*(sAlgo->spls->mA2x),*(sAlgo->spls->mA2s));
    aDS->setComputeBFunction(&bLDS);
  }else{
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
    aDS = new FirstOrderLinearTIDS(1,*DAE_X0,*DAE_A,*DAE_As);
    aDS->setMPtr(DAE_M);
    aDS->setComputeBFunction(&DAEbLDS);
  }
  cout<<"FirstOrderLinearTIDS with :"<<endl;
  aDS->display();
  DynamicalSystemsSet Inter_DS;
  Inter_DS.insert(aDS);

  //******BUILD THE RELATION
  aceMatrix* C= 0;
  aceMatrix* D= 0;
  aceMatrix* B= 0;
  if (ACE_FORMULATION_WITH_INVERSION){
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
  }else{//FORMULATION DAE
    C= new aceMatrix(m,dimX);
    D= new aceMatrix(m,m);
    B= new aceMatrix(dimX,m);
    C->setBlock(0,0,*(sAlgo->spls->mD2x));
    C->setBlock(0,ACEFdimX,*(sAlgo->spls->mD2zs));

    D->setBlock(0,0,*(sAlgo->spls->mD2l));

    B->setBlock(0,0,*(sAlgo->spls->mR));
    B->setBlock(ACEFdimX,0,*(sAlgo->spls->mB2l));
  }
  FirstOrderLinearR * aR = new FirstOrderLinearR(C,B);
  aR->setD(*D);
  if (ACE_FORMULATION_WITH_INVERSION){
    aR->setComputeEFunction(&eLDS);
  }else{
    aR->setComputeEFunction(&DAEeLDS);
  }
  cout<<"FirstOrderLinearR with :"<<endl;
  aR->display();

  //*****BUILD THE NSLAW
  NonSmoothLaw * aNSL=0;
  int NSLawSize=0;
  if (ACE_FORMULATION_WITH_INVERSION){
    NSLawSize=m+s;
    aNSL = new MixedComplementarityConditionNSL(m,s);
  }else{
    NSLawSize=m;
    aNSL = new ComplementarityConditionNSL(m);
  }

  //****BUILD THE INTERACTION
  Interaction * aI = new Interaction("MLCP",Inter_DS,1,NSLawSize,aNSL,aR);
  //****BUILD THE SYSTEM
  NonSmoothDynamicalSystem * aNSDS = new NonSmoothDynamicalSystem(aDS,aI);
  
  Model * aM = new Model(0,finalTime);
  aM->setNonSmoothDynamicalSystemPtr(aNSDS);
  TimeDiscretisation * aTD = new TimeDiscretisation(h,aM);
  TimeStepping *aS = new TimeStepping(aTD);
  
  //*****BUILD THE STEP INTEGRATOR
  OneStepIntegrator * aMoreau =0;
  if (ACE_FORMULATION_WITH_INVERSION){
    aMoreau = new Moreau(aDS,0.5,aS);
  }else{
    aMoreau = new Moreau2(aDS,0.5,aS);
  }
  
  NonSmoothSolver * mySolver = new NonSmoothSolver((*solverName),iparam,dparam,floatWorkingMem,intWorkingMem);

  //**** BUILD THE STEP NS PROBLEM
  MLCP * aMLCP =0;
  if (ACE_FORMULATION_WITH_INVERSION){
    aMLCP = new MLCP(aS,mySolver,"MLCP");
  }else{
    aMLCP = new MLCP2(aS,mySolver,"MLCP2");
  }
  aS->initialize();
  //  Alloc working mem for the solver
  if (ACE_SOLVER_TYPE==ACE_SOLVER_NUMERICS_DIRECT_ENUM ||
      ACE_SOLVER_TYPE==ACE_SOLVER_ENUM  ||
      ACE_SOLVER_TYPE==ACE_SOLVER_NUMERICS_DIRECT_SIMPLEX ||
      ACE_SOLVER_TYPE==ACE_SOLVER_NUMERICS_DIRECT_PATH){
    int aux =mlcp_driver_get_iwork(aMLCP->getNumericsMLCP(),mySolver->getNumericsSolverOptionsPtr());
    intWorkingMem = (int*)malloc(aux*sizeof(int));
    mySolver->getNumericsSolverOptionsPtr()->iWork = intWorkingMem;
    aux =mlcp_driver_get_dwork(aMLCP->getNumericsMLCP(),mySolver->getNumericsSolverOptionsPtr());
    floatWorkingMem = (double*)malloc(aux*sizeof(double));
    mySolver->getNumericsSolverOptionsPtr()->dWork = floatWorkingMem;
  }
  /*init solver*/
  mlcp_driver_init(aMLCP->getNumericsMLCP(),mySolver->getNumericsSolverOptionsPtr());

  
  SiconosVector * x = aDS->getXPtr();
  SiconosVector * y = aI->getYPtr(0);
  SiconosVector * lambda = aI->getLambdaPtr(0);


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
  
  for(int k = 0 ; k < N ; k++){
    if ((!Nfreq) || k % Nfreq == 0){
      cout<<"..."<<cmp<<endl;
      cmp++;
    }
    // solve ... 
    aS->computeOneStep();

    if ((!NData) || k % NData == 0){
      ParserInitPrintElem();
      //pout <<k+1<<"\t"<<(k+1)*h;
      pout <<(k+1)*h;
      while(ParserGetPrintElem((void**)&pPrint)){
	pout<<"\t\t";
	if (ACE_FORMULATION_WITH_INVERSION){
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
    }

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
