/************************************************************************
  			equationvd.cpp 
**************************************************************************/

#include "equationvd.h"

equationVD::equationVD(char * name){
  mName = name;
}
void equationVD::print(){
  printf("VD");
  if (mName)
    printf(mName);
  printf("\t");
  equation::print();
}
