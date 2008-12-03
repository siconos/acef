/************************************************************************
  			equationvd.cpp 
**************************************************************************/

#include "equationvd.h"

equationVD::equationVD(char * name){
  mName = name;
}
void equationVD::print(){
  if (mLine > -1)
    printf("VD%d",mLine);
  else
    printf("VD");
  if (mName)
    printf(mName);
  printf("\t");
  equation::print();
}
