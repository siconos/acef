#include "ACEF.h"

int main(int argc, char **argv){
  if (argc<6){
    printf("usage : toto file.cir ENUM|SIMPLEX|PATH 10 DENSE|SPARSE INV/NOINV [FIX/ADAPT]\n");
    return 0;
  }
  if (!strcmp(argv[2],"ENUM")){
    ACE_SOLVER_TYPE = ACE_SOLVER_ENUM;
  }else if(!strcmp(argv[2],"SIMPLEX")){
    ACE_SOLVER_TYPE = ACE_SOLVER_SIMPLEX;
  }else if(!strcmp(argv[2],"PATH")){
    ACE_SOLVER_TYPE = ACE_SOLVER_PATH;
  }else{
    printf("usage : toto file.cir ENUM|SIMPLEX|PATH 10\n");
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
    ACE_FORMULATION_WITH_INVERSION=0;
  if (argc > 6 && !strcmp(argv[6],"ADAPT")){
    ACE_MAX_LOCAL_ERROR=1e-4;
    ACE_WITH_ADAPTATIVE_TIME_STEPPING=1;
  }
    
  ACE_times[ACE_TIMER_MAIN].start();
  ACE_INIT();
  algo *a=new algo(argv[1]);
  a->perform();
  a-> simulate();
  ACE_times[ACE_TIMER_MAIN].stop();
  ACE_PRINT_TIME();
  return 0;
 
}
