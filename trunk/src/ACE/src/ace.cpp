#include "extern.h"
#include "ace.h"
#include <math.h>
#include <cstdlib>

char ACE_name[12]="ace";
aceTime ACE_times[ACE_TIMER_LAST];
int ACE_SOLVER_TYPE=0;
int ACE_MUET_LEVEL=0;
UBLAS_TYPE ACE_MAT_TYPE=DENSE;
int ACE_FORMULATION_WITH_INVERSION=1;
int ACE_WITH_ADAPTATIVE_TIME_STEPPING=0;
ACE_DOUBLE ACE_MAX_LOCAL_ERROR=1e-3;
static ofstream* ACE_LOG_FILE=0;
int ACE_CUR_STEP=0;
int ACE_CMP_ADAT[ACE_NB_ADAPT_STEP+1];
int  ACE_MOS_NB_HYP=3;
ACE_DOUBLE ACE_MOS_POWER_SUPPLY=3.0;

using namespace std;


ofstream& ACE_GET_LOG_STREAM(){
  return *ACE_LOG_FILE;
}

void ACE_INIT(){
  ACE_INIT_TIME();
  ACE_LOG_FILE = new ofstream("ACE.log");
}
void ACE_STOP(){
  delete ACE_LOG_FILE;
}
void ACE_INIT_TIME(){
  ACE_times[ACE_TIMER_MAIN].setName("main ");
  ACE_times[ACE_TIMER_SOLVE_LU].setName("compute LU ");
  ACE_times[ACE_TIMER_SOLVE_ENUM].setName("solver enum ");
  ACE_times[ACE_TIMER_SOLVE_SIMPLEX].setName("solver simplex ");
  ACE_times[ACE_TIMER_SOLVE_PATH].setName("solver path ");
  ACE_times[ACE_TIMER_EQUATION].setName("equation formulation ");
  ACE_times[ACE_TIMER_SIMPLEX_FIRST].setName("simplex first getsolution ");
  ACE_times[ACE_TIMER_SIMPLEX_GUESS].setName("simplex try previous guess ");
  ACE_times[ACE_TIMER_SIMPLEX_TREE].setName("simplex tree exploration ");
  ACE_times[ACE_TIMER_SIMPLEX_TRY_NODE].setName("simplex try node ");
  ACE_times[ACE_TIMER_SOLVE_GUESS].setName("solver guess ");
  ACE_times[ACE_TIMER_SOLVER].setName("mlcp solver  ");
  ACE_times[ACE_TIMER_DIRECT].setName("solver direct  ");
  ACE_times[ACE_TIMER_LU_DIRECT].setName(" lu direct  ");
  ACE_times[ACE_TIMER_LS_STEP].setName("ls step  ");
  ACE_times[ACE_TIMER_COMPUTE_VAR].setName("lcp --> netlist  ");
  ACE_times[ACE_TIMER_SIMULATION].setName("Simulation  ");
  ACE_times[ACE_TIMER_PROD_MAT].setName("Prod mat  ");
  ACE_times[ACE_TIMER_TEST].setName("test  ");
  ACE_times[ACE_TIMER_TEST_1].setName("computedxdt  ");
  ACE_times[ACE_TIMER_TEST_2].setName("test_2  ");
  ACE_times[ACE_TIMER_TEST_3].setName("test_3  ");
  ACE_times[ACE_TIMER_TEST_4].setName("test_4  ");
  ACE_times[ACE_TIMER_TEST_5].setName("test_5  ");
  ACE_times[ACE_TIMER_TEST_6].setName("test_6  ");
  ACE_times[ACE_TIMER_TEST_7].setName("test_7  ");
  ACE_times[ACE_TIMER_TEST_8].setName("test_8  ");
  ACE_times[ACE_TIMER_TEST_9].setName("test_9  ");
  ACE_times[ACE_TIMER_TEST_10].setName("test_10  ");
  ACE_times[ACE_TIMER_TEST_11].setName("test_11  ");
  
  
  
  
}
void ACE_STOP_SOLVER_TIME(){
  ACE_times[ACE_TIMER_SOLVE_ENUM].stop();
  ACE_times[ACE_TIMER_SOLVE_SIMPLEX].stop();
  ACE_times[ACE_TIMER_SOLVE_PATH].stop();
  ACE_times[ACE_TIMER_SOLVE_GUESS].stop();
  ACE_times[ACE_TIMER_SOLVER].stop();
}
void ACE_PRINT_TIME(){
  for (int i=0;i <ACE_TIMER_LAST;i++)
    ACE_times[i].print();
}
bool ACE_IS_NULL(ACE_DOUBLE d){
  return (fabs(d)<ACE_INF);
}
void ACE_MESSAGE( char * mess,int level){
  if (ACE_MUET_LEVEL>level)
    return;
  (*ACE_LOG_FILE)<<"ACE MESSAGE :"<<mess<<endl;
  ACE_LOG_FILE->flush();
}

void ACE_ERROR(char * mess){
  (*ACE_LOG_FILE)<<"ACE ERROR: "<<mess<<endl;
  ACE_LOG_FILE->flush();
  ACE_STOP();
  exit(1);

}
void ACE_WARNING(char * mess){
  (*ACE_LOG_FILE)<<"ACE WARNING: "<<mess<<endl;
  ACE_LOG_FILE->flush();
}
void ACE_INTERNAL_ERROR(char *mess){
  (*ACE_LOG_FILE)<<"ACE INTERNAL ERROR: "<<mess<<endl;
  ACE_LOG_FILE->flush();
  ACE_STOP();
  exit(1);

}
void ACE_INTERNAL_WARNING(char *mess){
  (*ACE_LOG_FILE)<<"ACE INTERNAL WARNING: "<<mess<<endl;
  ACE_LOG_FILE->flush();

}
void ACE_TYPE_TO_CHAR(int type,char* name){
  switch(type){
    case ACE_TYPE_NO:
      strcpy(name,"NO");
      break;
    case ACE_TYPE_RES:
      strcpy(name,"Resistor");
      break;
    case ACE_TYPE_CAP:
      strcpy(name,"Capacitor");
      break;
    case ACE_TYPE_IND:
      strcpy(name,"Inductor");
      break;
    case ACE_TYPE_VSRC:
      strcpy(name,"Vsource");
      break;
    case ACE_TYPE_ISRC:
      strcpy(name,"Isource");
      break;
    case ACE_TYPE_DIO:
      strcpy(name,"Diode");
      break;
    case ACE_TYPE_COMP:
      strcpy(name,"Comp");
      break;
    case ACE_TYPE_U:
      strcpy(name,"U");
      break;
    case ACE_TYPE_V:
      strcpy(name,"V");
      break;
    case ACE_TYPE_I:
      strcpy(name,"I");
      break;
    case ACE_TYPE_MOS:
      strcpy(name,"MOS");
      break;
    case ACE_TYPE_VCVS:
      strcpy(name,"VCVS");
      break;
    case ACE_TYPE_ARB:
      strcpy(name,"ARB");
      break;
    case ACE_TYPE_VCCS:
      strcpy(name,"VCCS");
    case ACE_TYPE_BJT:
      strcpy(name,"BJT");
      break;
  default:
    strcpy(name,"????");
  }
}
void ACE_CHECK_IERROR(bool b,char* mess){
  if (!b)
    ACE_INTERNAL_ERROR(mess);
}
void ACE_CHECK_IWARNING(bool b,char* mess){
  if (!b)
    ACE_INTERNAL_WARNING(mess);
}
void ACE_CHECK_WARNING(bool b,char* mess){
  if (!b)
    ACE_WARNING(mess);
}
void ACE_CHECK_ERROR(bool b,char* mess){
  if (!b)
    ACE_ERROR(mess);
}
void startTime(int ID){
#ifdef ACE_WITH_TIMER
  ;
#endif
}
double stopTime(int ID){
#ifdef ACE_WITH_TIMER
  ;
#endif
  return 0;
}
