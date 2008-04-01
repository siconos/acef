#include"ace.h"
#include "algo.h"

#include "/usr/local/include/Siconos/Kernel/MLCP.h"
#include "SiconosKernel.h"

static  algo *sAlgo=0;

void (bLDS) (double t, unsigned int N, double* b, unsigned int z, double*zz){
  sAlgo->preparStep(t);
  sAlgo->sls.extractDynamicSystemSource();
  sAlgo->sls.computeDynamicSystemSource();
//   printf("bLDS %lf :\n",t);
//   printf("mA2s");
//   cout<<(*(sAlgo->sls.mA2s));
  for (int i=0;i<N;i++){
    b[i]= sAlgo->sls.mA2s->getValue(i);
  }
}
void (eLDS) (double t, unsigned int N, double* e, unsigned int z, double*zz){
  sAlgo->preparStep(t);
  sAlgo->sls.extractInteractionSource();
  sAlgo->sls.computeInteractionSource();
  int s = sAlgo->sls.mB2s->dimRow;
  int m = sAlgo->sls.mD2s->dimRow;
//   printf("eLDS %lf :\n",t);
//   printf("mB2s");
//   cout<<(*(sAlgo->sls.mB2s));
//   printf("mD2s");
//   cout<<(*(sAlgo->sls.mD2s));
  for (int i=0;i<s ;i++){
    e[i]= sAlgo->sls.mB2s->getValue(i);
  }
  for (int i=0;i<m ;i++){
    e[i+s]= sAlgo->sls.mD2s->getValue(i);
  }
}


int main(int argc, char **argv){
  if (argc<5){
    printf("usage : toto file.cir ENUM|SIMPLEX|PATH|DIR_ENUM 10 DENSE|SPARSE\n");
    return 0;
  }
  string solverDirEnum = "DIRECT_ENUM" ;
  string solverEnum = "ENUM" ;
  string solverSimplex = "SIMPLEX" ;
  string solverPath = "PATH" ;
  string * solverName =0;
  // One Step non smooth problem
  IntParameters iparam(10);
  DoubleParameters dparam(10);
  
  double* floatWorkingMem = 0;
  int * intWorkingMem= 0;

  if (!strcmp(argv[2],"ENUM")){
    ACE_SOLVER_TYPE = ACE_SOLVER_ENUM;
  }else if(!strcmp(argv[2],"SIMPLEX")){
    ACE_SOLVER_TYPE = ACE_SOLVER_SIMPLEX;
  }else if(!strcmp(argv[2],"PATH")){
    ACE_SOLVER_TYPE = ACE_SOLVER_PATH;
  }else if(!strcmp(argv[2],"DIR_ENUM")){
    ACE_SOLVER_TYPE=ACE_SOLVER_NUMERICS_DIRECT_ENUM;
  }else{
    printf("param2 must be : ENUM|SIMPLEX|PATH|DIR_ENUM \n");
    return 0;
  }
  if (!strcmp(argv[3],"0"))
    ACE_MUET_LEVEL=0;
  else if (!strcmp(argv[3],"1"))
    ACE_MUET_LEVEL=1;
  else if (!strcmp(argv[3],"2"))
    ACE_MUET_LEVEL=2;
  else
    ACE_MUET_LEVEL=10;
  
  if (!strcmp(argv[4],"SPARSE"))
    ACE_MAT_TYPE=SPARSE;
  else
    ACE_MAT_TYPE=DENSE;
    
  ACE_times[ACE_TIMER_MAIN].start();
  ACE_INIT();
  sAlgo=new algo(argv[1]);
  sAlgo->perform();
  int dimX = sAlgo->sls.mDimx;
  int s=sAlgo->sls.mB2zs->getDimRow();
  int m=sAlgo->sls.mD2l->getDimRow();
  aceVector* X0 = new aceVector(dimX);
  
  sAlgo->sls.allocForInitialValue();
  sAlgo->sls.readInitialValue();
  double h=0.00001;
  h=sAlgo->sls.mH;
  double finalTime = h*sAlgo->sls.mStepNumber ;
  FirstOrderLinearTIDS * aDS = new FirstOrderLinearTIDS(1,*(sAlgo->sls.mxti),*(sAlgo->sls.mA2x),*(sAlgo->sls.mA2s));
  aDS->setComputeBFunction(&bLDS);
  cout<<"FirstOrderLinearTIDS with :"<<endl;
  aDS->display();
  DynamicalSystemsSet Inter_DS;
  Inter_DS.insert(aDS);

  aceMatrix* C= new aceMatrix(s+m,dimX);
  aceMatrix* D= new aceMatrix(s+m,s+m);
  aceMatrix* B= new aceMatrix(dimX,s+m);
  C->setBlock(0,0,*(sAlgo->sls.mB2x));
  C->setBlock(s,0,*(sAlgo->sls.mD2x));
  D->setBlock(0,0,*(sAlgo->sls.mB2zs));
  D->setBlock(0,s,*(sAlgo->sls.mB2l));
  D->setBlock(s,0,*(sAlgo->sls.mD2zs));
  D->setBlock(s,s,*(sAlgo->sls.mD2l));
  B->setBlock(0,0,*(sAlgo->sls.mA2zs));
  B->setBlock(0,s,*(sAlgo->sls.mR));
  aceVector *q= new aceVector(s+m);
  q->setBlock(0,*(sAlgo->sls.mB2s));
  q->setBlock(s,*(sAlgo->sls.mD2s));
  FirstOrderLinearR * aR = new FirstOrderLinearR(C,B);
  aR->setD(*D);
  aR->setComputeEFunction(&eLDS);
  cout<<"FirstOrderLinearR with :"<<endl;
  aR->display();
  MixedComplementarityConditionNSL * aNSL = new MixedComplementarityConditionNSL(m,s);
  Interaction * aI = new Interaction("MLCP",Inter_DS,1,m+s,aNSL,aR);
  NonSmoothDynamicalSystem * aNSDS = new NonSmoothDynamicalSystem(aDS,aI);
  Model * aM = new Model(0,finalTime);
  aM->setNonSmoothDynamicalSystemPtr(aNSDS);
  TimeDiscretisation * aTD = new TimeDiscretisation(h,aM);
  TimeStepping *aS = new TimeStepping(aTD);
  Moreau * aMoreau = new Moreau(aDS,0.5,aS);
  
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
    iparam[0]=101;
    dparam[0]=1e-9;
    solverName = &solverPath;
  }else if (ACE_SOLVER_TYPE==ACE_SOLVER_NUMERICS_DIRECT_ENUM){
    iparam[0] = 0; // enum verbose
    iparam[5] = 2; // nb config
    iparam[6] = 0; // direct verbose
    dparam[0] = 0.0000001; // Tolerance
    dparam[5] = 0.0000001; // Tolerance
    solverName = &solverDirEnum;
  }

    
  NonSmoothSolver * mySolver = new NonSmoothSolver((*solverName),iparam,dparam,floatWorkingMem,intWorkingMem);
  
  MLCP * aMLCP = new MLCP(aS,mySolver,"MLCP");
  aS->initialize();
  //  Alloc working mem
  if (ACE_SOLVER_TYPE==ACE_SOLVER_NUMERICS_DIRECT_ENUM || ACE_SOLVER_TYPE==ACE_SOLVER_ENUM  ){
    int aux =mlcp_driver_get_iwork(aMLCP->getNumericsMLCP(),mySolver->getNumericsSolverOptionsPtr());
    intWorkingMem = (int*)malloc(aux*sizeof(int));
    mySolver->getNumericsSolverOptionsPtr()->iWork = intWorkingMem;
    aux =mlcp_driver_get_dwork(aMLCP->getNumericsMLCP(),mySolver->getNumericsSolverOptionsPtr());
    floatWorkingMem = (double*)malloc(aux*sizeof(double));
    mySolver->getNumericsSolverOptionsPtr()->dWork = floatWorkingMem;
  }
    
  mlcp_driver_init(aMLCP->getNumericsMLCP(),mySolver->getNumericsSolverOptionsPtr());

  
  SiconosVector * x = aDS->getXPtr();
  SiconosVector * y = aI->getYPtr(0);
  SiconosVector * lambda = aI->getLambdaPtr(0);


  unsigned int count = 0; // events counter. 
  // do simulation while events remains in the "future events" list of events manager. 
  cout << " ==== Start of  simulation ====" << endl;
  ofstream pout("acefSimu.dat");

  int N = sAlgo->sls.mStepNumber;
  dataPrint * pPrint;
  initPrintElem();
  pout<<"Index\ttime";
  while(getPrintElem((void**)&pPrint)){
    pout<<"\t\t";
    pout<<pPrint->name;
  }
  pout<<endl<<endl;
  /*print t0*/
  initPrintElem();
  pout <<0<<"\t"<<0;
  while(getPrintElem((void**)&pPrint)){
    pout<<"\t\t";
    double aux = sAlgo->sls.mzsti->getValue(pPrint->node1-1);
    if (pPrint->node2 >0)
      aux -=  sAlgo->sls.mzsti->getValue(pPrint->node2-1);
    pout << aux;
  }
  pout<<endl;
  
  for(int k = 0 ; k < N ; k++){
    // solve ... 
    aS->computeOneStep();
    

    initPrintElem();
    pout <<k+1<<"\t"<<(k+1)*h;
    while(getPrintElem((void**)&pPrint)){
      pout<<"\t\t";
      double aux = (*lambda)(pPrint->node1-1);
      if (pPrint->node2 >0)
	aux -= (*lambda)(pPrint->node2-1);
      pout << aux;
    }
    pout<<endl;

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
  return 0;
 
}
