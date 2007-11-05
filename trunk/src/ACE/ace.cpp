#include "extern.h"
#include "ace.h"
#include <math.h>
#include <cstdlib>

char ACE_name[12]="ace_name";

bool ACE_IS_NULL(ACE_DOUBLE d){
  return (fabs(d)<ACE_INF);
}
void ACE_MESSAGE(char * mess){
#ifndef NO_MESSAGE
  printf("ACE MESSAGE :");
  printf(mess);
#endif
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
      strcpy(name,"Comparator?");
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
