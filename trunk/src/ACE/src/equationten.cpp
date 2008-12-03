/************************************************************************
  			equationten.cpp 
**************************************************************************/

#include "equationten.h"
equationTEN::equationTEN(){;}
void equationTEN::print(){
  if (mLine > -1)
    printf("TEN%d\t",mLine);
  else
    printf("TEN\t");
 equation::print();
}
