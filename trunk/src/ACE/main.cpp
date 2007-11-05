#include"ace.h"
#include "algo.h"
int main(int argc, char **argv){
  algo *a=new algo(argv[1]);
  a->perform();
  return 0;
 
}
