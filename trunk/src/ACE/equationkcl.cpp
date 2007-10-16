/************************************************************************
  			equationkcl.cpp 
**************************************************************************/

#include "equationkcl.h"

equationKCL::equationKCL(int n){
  mNode = n;
}
void equationKCL::print(){
 printf("KCL%d\t",mNode);
 equation::print();
}

