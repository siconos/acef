/************************************************************************
  			equation_nl.cpp 
**************************************************************************/

#include "equation_nl.h"
equation_NL::equation_NL(){;}
void equation_NL::print(){
  if (mLine > -1)
    printf("NL%d\t",mLine);
  else
    printf("NL\t");
  printf("\n");
}
