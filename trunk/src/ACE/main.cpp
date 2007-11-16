#include"ace.h"
#include "algo.h"
int main(int argc, char **argv){
  if (argc<2){
    printf("usage : toto file.cir\n");
    return 0;
  }
  algo *a=new algo(argv[1]);
  a->perform();
  return 0;
 
}
