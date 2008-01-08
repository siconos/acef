#include "extern.h"
#include "ace.h"
#include <math.h>
#include <cstdlib>

char ACE_name[12]="ace";
aceTime ACE_times[ACE_TIMER_LAST];
int ACE_SOLVER_TYPE=0;
int ACE_MUET_LEVEL=0;

void ACE_INIT(){
  ACE_INIT_TIME();
}
void ACE_INIT_TIME(){
  ACE_times[ACE_TIMER_MAIN].setName("main ");
  ACE_times[ACE_TIMER_SOLVE_GUESS].setName("solver guess ");
  ACE_times[ACE_TIMER_SOLVE_ENUM].setName("solver enum ");
  ACE_times[ACE_TIMER_SOLVE_SIMPLEX].setName("solver simplex ");
  ACE_times[ACE_TIMER_SOLVE_PATH].setName("solver path ");
  ACE_times[ACE_TIMER_EQUATION].setName("equation formulation ");
  ACE_times[ACE_TIMER_DIRECT].setName("solver direct ");
}
void ACE_STOP_SOLVER_TIME(){
  ACE_times[ACE_TIMER_SOLVE_GUESS].stop();
  ACE_times[ACE_TIMER_SOLVE_ENUM].stop();
  ACE_times[ACE_TIMER_SOLVE_SIMPLEX].stop();
  ACE_times[ACE_TIMER_SOLVE_PATH].stop();
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
  printf("ACE MESSAGE :");
  printf(mess);
}

void ACE_ERROR(char * mess){
  printf("ACE ERROR: ");
  printf(mess);
  printf("\n");
  exit(1);

}
void ACE_WARNING(char * mess){
  printf("ACE WARNING: ");
  printf(mess);
  printf("\n");
}
void ACE_INTERNAL_ERROR(char *mess){
  printf("ACE INTERNAL ERROR: ");
  printf(mess);
  printf("\n");
  exit(1);

}
void ACE_INTERNAL_WARNING(char *mess){
  printf("ACE INTERNAL WARNING: ");
  printf(mess);
  printf("\n");

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
