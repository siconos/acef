/************************************************************************
  			component.cpp 
**************************************************************************/


#include "unknown.h"
#include "equation.h"


component::component(){
  mU=0;
  mI=0;
  mEquation=0;
  mType=ACE_TYPE_NO;
  mNodePos=-1;
  mNodeNeg=-1;
  mName=0;

}


void component::stamp () {
  
}
/**
 * 
 */
void component::addUnknowns () {
  
}
/**
 * 
 */
void component::addEquations () {
  
}
void component::print () {
  char name[ACE_CHAR_LENGTH];
  ACE_TYPE_TO_CHAR(mType,name);
  printf("component type : %s \n",name);
  if (mName)
    printf("\t Name : %s\n",mName);
  printf("\tnodePos, nodeNeg : %d %d\n",mNodePos,mNodeNeg);
  
  
}


