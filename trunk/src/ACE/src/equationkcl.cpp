/************************************************************************
  			equationkcl.cpp 
**************************************************************************/

#include "equationkcl.h"

equationKCL::equationKCL(int n){
  mNode = n;
}
void equationKCL::print(){
 printf("KCL%d",mNode);
 if (mIsDyn)
   printf("*%d",mLine);
 printf("\t");
 equation::print();
}

