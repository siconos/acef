#include"ace.h"
#include "algo.h"
int main(){
  if (!initParserLibrary()){
    ACE_INTERNAL_ERROR("initParserLibrary");
  }
  readFile("./LR.cir");
  printCircuit();
  algo *a=new algo();
  a->perform();
  a->printComponents();
  stopParserLibrary();
  return 0;
 
}
