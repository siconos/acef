#include "ACEF.h"

int main(int argc, char **argv){
  /*MOS Parameters*/
  /*  ACE_MOS_NB_HYP=5;*/
  ACE_MOS_POWER_SUPPLY=3.0;
  ACE_FORMULATION=ACE_FORMULATION_SEMI_EXPLICT;
  

  if (argc<6){
    printf("usage : toto file.cir ENUM|SIMPLEX|PATH|FB 10 DENSE|SPARSE INV/NOINV/MNA/STAMP/SE1 [FIX/ADAPT]\n");
    return 0;
  }
  if (!strcmp(argv[2],"ENUM")){
    ACE_SOLVER_TYPE = ACE_SOLVER_ENUM;
  }else if(!strcmp(argv[2],"SIMPLEX")){
    ACE_SOLVER_TYPE = ACE_SOLVER_SIMPLEX;
  }else if(!strcmp(argv[2],"PATH")){
    ACE_SOLVER_TYPE = ACE_SOLVER_PATH;
  }else if(!strcmp(argv[2],"FB")){
    ACE_SOLVER_TYPE = ACE_SOLVER_FB;
  }else{
    printf("usage : toto file.cir ENUM|SIMPLEX|PATH|FB 10 DENSE|SPARSE INV/NOINV/MNA/STAMP/SE1 [FIX/ADAPT]\n");
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
  
  if (!strcmp(argv[5],"NOINV"))
    ACE_FORMULATION=ACE_FORMULATION_WITHOUT_INVERT;
  if (!strcmp(argv[5],"MNA"))
    ACE_FORMULATION=ACE_FORMULATION_MNA;
  if (!strcmp(argv[5],"MNA_V"))
    ACE_FORMULATION=ACE_FORMULATION_MNA_V;
  if (!strcmp(argv[5],"STAMP"))
    ACE_FORMULATION=ACE_FORMULATION_STAMP_ONLY;
  if (!strcmp(argv[5],"SE1")){
        printf("SE1 not managed\n");
        return 0;
	//    ACE_FORMULATION=ACE_FORMULATION_SE1;
  }
  
  if (argc > 6 && !strcmp(argv[6],"ADAPT")){
    ACE_RTOL_LOCAL=1e-7;
    ACE_ATOL_LOCAL=1e-7;
    ACE_WITH_ADAPTATIVE_TIME_STEPPING=1;
  }
    
  ACE_times[ACE_TIMER_MAIN].start();
  ACE_INIT();
  algo *a=new algo(argv[1]);
  a->perform();
  a->simulate();
  ACE_times[ACE_TIMER_MAIN].stop();
  ACE_PRINT_TIME();
  ACE_STOP();
  return 0;
 
}
