#include"ace.h"
#include "algo.h"
int main(int argc, char **argv){
  if (!initParserLibrary()){
    ACE_INTERNAL_ERROR("initParserLibrary");
  }
  readFile(argv[1]);
  printCircuit();
  algo *a=new algo();
  a->perform();
  stopParserLibrary();
  return 0;
 
}
