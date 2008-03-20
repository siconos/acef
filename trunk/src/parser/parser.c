#include "parser.h"
#include <dlfcn.h>
#include <stdio.h>

typedef void (*fct_print)();
typedef int (*fct_init)();
typedef int (*fct_nextComponent)(void*);
typedef int (*fct_readFile)(char *);
typedef int (*fct_initComponentList)(char *);
typedef int (*fct_getNbElementsOfType)(char *);
typedef int (*fct_getSourceValue)(char *,void *,double *);
typedef int (*fct_computeSourcesValues)(double);
typedef int (*fct_getTransValues)(double *,double *,double *);
typedef int (*fct_getICvalue)(int *,int *,double *);
typedef int (*fct_getPrintElem)(void **);
typedef int (*fct_initSimulation)(int,double);

static fct_print ptr_print=0;
static fct_readFile ptr_readFile=0;
static fct_nextComponent ptr_nextComponent=0;
static fct_initComponentList ptr_initComponentList=0;
static fct_getNbElementsOfType ptr_getNbElementsOfType=0;
static fct_getSourceValue ptr_getSourceValue=0;
static fct_computeSourcesValues ptr_computeSourcesValues=0;
static fct_getTransValues ptr_getTransValues=0;
static fct_init ptr_initICvalue=0;
static fct_getICvalue ptr_getICvalue=0;
static fct_init ptr_initPrintElem=0;
static fct_getPrintElem ptr_getPrintElem=0;
static fct_initSimulation ptr_initSimulation=0;

static void* module=0;


void printCircuit(){
  (*ptr_print)();
}

int initParserLibrary(){
  const char *error;
  module = dlopen("perform.so", RTLD_LAZY);
  if (!module) {
   fprintf(stderr, "Couldn't open perform.so: %s\n",dlerror());
   module = dlopen("/home/bipop/bonnefon/ACE/lib/perform.so", RTLD_LAZY);
   if (!module) {
     fprintf(stderr, "Couldn't open .../perform.so: %s\n",dlerror());
     return (0);
   }
  }
  ptr_print = dlsym(module, "MEprint");
  if ((error = dlerror())) {
   fprintf(stderr, "Couldn't find MEprint: %s\n", error);
   return (0);
  }
  ptr_initComponentList = dlsym(module, "initComponentList");
  if ((error = dlerror())) {
   fprintf(stderr, "Couldn't find initComponentList: %s\n", error);
   return (0);
  }
  ptr_nextComponent = dlsym(module, "nextComponent");
  if ((error = dlerror())) {
   fprintf(stderr, "Couldn't find nextComponent: %s\n", error);
   return (0);
  }
  ptr_getNbElementsOfType = dlsym(module, "getNbElementsOfType");
  if ((error = dlerror())) {
   fprintf(stderr, "Couldn't find getNbElementsOfType: %s\n", error);
   return (0);
  }
  ptr_readFile = dlsym(module, "readFile");
  if ((error = dlerror())) {
   fprintf(stderr, "Couldn't find readFile: %s\n", error);
   return (0);
  }
  
  ptr_getSourceValue = dlsym(module, "getSourceValue");
  if ((error = dlerror())) {
   fprintf(stderr, "Couldn't find getSourceValue: %s\n", error);
   return (0);
  }
  ptr_computeSourcesValues = dlsym(module, "computeSourcesValues");
  if ((error = dlerror())) {
   fprintf(stderr, "Couldn't find computeSourcesValues: %s\n", error);
   return (0);
  }
  
  ptr_getTransValues = dlsym(module, "getTransValues");
  if ((error = dlerror())) {
   fprintf(stderr, "Couldn't find getTransValues: %s\n", error);
   return (0);
  }

  ptr_initICvalue = dlsym(module, "initICvalue");
  if ((error = dlerror())) {
   fprintf(stderr, "Couldn't find initICvalue: %s\n", error);
   return (0);
  }

  ptr_getICvalue = dlsym(module, "getICvalue");
  if ((error = dlerror())) {
   fprintf(stderr, "Couldn't find getICvalue: %s\n", error);
   return (0);
  }

  ptr_initPrintElem = dlsym(module, "initPrintElem");
  if ((error = dlerror())) {
   fprintf(stderr, "Couldn't find initPrintElem: %s\n", error);
   return (0);
  }
  
  ptr_getPrintElem = dlsym(module, "getPrintElem");
  if ((error = dlerror())) {
   fprintf(stderr, "Couldn't find getPrintElem: %s\n", error);
   return (0);
  }
  ptr_initSimulation = dlsym(module, "initSimulation");
  if ((error = dlerror())) {
   fprintf(stderr, "Couldn't find initSimulation: %s\n", error);
   return (0);
  }

  return 1;
}

void stopParserLibrary(){
    dlclose(module);
}
int readFile(char* file){
  return (*ptr_readFile)(file);
}
int getNbElementsOfType(char *type){
  return (*ptr_getNbElementsOfType)(type);
}

int initComponentList(char *type){
  return (*ptr_initComponentList)(type);
}
int nextComponent(void * data){
  return (*ptr_nextComponent)(data);
}

int getSourceValue(char *type,void* id,double* value){
  return (*ptr_getSourceValue)(type,id,value);
}
int computeSourcesValues(double time){
  return (*ptr_computeSourcesValues)(time);
}
int getTransValues(double * step, double * stop, double * start){
  return (*ptr_getTransValues)(step,stop,start);
}
int initICvalue(){
  return (*ptr_initICvalue)();

}
int getICvalue(int * numNode,int * icGiven, double * icValue) {
  return (*ptr_getICvalue)(numNode,icGiven,icValue);
}
int initPrintElem(){
  return (*ptr_initPrintElem)();

}
int getPrintElem(void ** p){
  return (*ptr_getPrintElem)(p);
}
int initSimulation(int type,double val){
  return (*ptr_initSimulation)(type,val);
}
