#include "extern.h"
#include "ace.h"
#include <math.h>
#include <cstdlib>

char ACE_name[12]="ace";
aceTime ACE_times[ACE_TIMER_LAST];
int ACE_SOLVER_TYPE=0;
int ACE_MUET_LEVEL=0;
UBLAS_TYPE ACE_MAT_TYPE=DENSE;
int ACE_FORMULATION=ACE_FORMULATION_SEMI_EXPLICT;
int ACE_WITH_ADAPTATIVE_TIME_STEPPING=0;
ACE_DOUBLE ACE_ATOL_LOCAL=1e-3;
ACE_DOUBLE ACE_RTOL_LOCAL=1e-2;
static ofstream* ACE_LOG_FILE=0;
static ofstream* ACE_LOG1_FILE=0;
int ACE_CUR_STEP=0;
int ACE_CMP_ADAT[ACE_NB_ADAPT_STEP+1];
int ACE_USE_NL_MOS = 0;
int  ACE_MOS_NB_HYP=5;
ACE_DOUBLE ACE_MOS_POWER_SUPPLY=3.0;
ACE_DOUBLE ACE_DIODE_THRESHOLD=0;

double ACE_THETA_X=0.5;
double ACE_THETA_ZS=1.0;
double ACE_THETA_P=1.0;


using namespace std;


ofstream& ACE_GET_LOG_STREAM(){
  return *ACE_LOG_FILE;
}
ofstream& ACE_GET_LOG1_STREAM(){
  return *ACE_LOG1_FILE;
}

void ACE_INIT(){
  ifstream foption("./ace.opt");
  char opt_name[128];
  ACE_INIT_TIME();
  ACE_LOG_FILE = new ofstream("ACE.log");
  ACE_LOG1_FILE = new ofstream("ACE1.log");
  if (foption.is_open()){
    while(! foption.eof()){
      foption >> opt_name ;
      if (!strcmp("DIO_TH",opt_name))
	foption >> ACE_DIODE_THRESHOLD;
      else if (!strcmp("MOS_NB_HYP",opt_name))
	foption>>ACE_MOS_NB_HYP;
      else if (!strcmp("MOS_POWER",opt_name))
	foption>>ACE_MOS_POWER_SUPPLY;
      else if (!strcmp("ACE_ATOL",opt_name))
	foption>>ACE_ATOL_LOCAL;
      else if (!strcmp("ACE_RTOL",opt_name))
	foption>>ACE_RTOL_LOCAL;
      else if (!strcmp("ACE_USE_NL_MOS",opt_name))
	foption>>ACE_USE_NL_MOS;
      else if (!strcmp("ACE_THETA_X",opt_name))
	foption>>ACE_THETA_X;
      else if (!strcmp("ACE_THETA_ZS",opt_name))
	foption>>ACE_THETA_ZS;
      else if (!strcmp("ACE_THETA_P",opt_name))
	foption>>ACE_THETA_P;
      else
	ACE_WARNING("UNKNOWN OPTION\n");
    
    }
    foption.close();
  }else
      ACE_WARNING("NONE OPTION FILE\n");
  
}
void ACE_STOP(){
  delete ACE_LOG_FILE;
  delete ACE_LOG1_FILE;
}
void ACE_INIT_TIME(){
  ACE_times[ACE_TIMER_MAIN].setName("main ");
  ACE_times[ACE_TIMER_SOLVE_NUMERICS].setName("solver numerics ");
  ACE_times[ACE_TIMER_EQUATION].setName("equation formulation ");
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
  ACE_times[ACE_TIMER_LS_STEP].setName("LS step  ");
  
  
  
  
}
void ACE_STOP_SOLVER_TIME(){
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
    case ACE_TYPE_RELAY:
      strcpy(name,"Relay");
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
    case ACE_TYPE_MOS_NL:
      strcpy(name,"MOS_NL");
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
unsigned int ACE_GET_ID(char* type,char* name){
  return ParserGetId(type,name);
}
void ACE_SET_INPUT(char* type,unsigned int id, double v){
  ParserSetIputFromId(type,id,v);
}

