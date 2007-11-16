#include "parser.h"
#include <dlfcn.h>
#include <stdio.h>

typedef void (*fct_print)();
typedef int (*fct_nextComponent)(void*);
typedef int (*fct_readFile)(char *);
typedef int (*fct_initComponentList)(char *);
typedef int (*fct_getNbElementsOfType)(char *);
typedef int (*fct_getSourceValue)(char *,void *,double *);
typedef int (*fct_computeSourcesValues)(double);

static fct_print ptr_print=0;
static fct_readFile ptr_readFile=0;
static fct_nextComponent ptr_nextComponent=0;
static fct_initComponentList ptr_initComponentList=0;
static fct_getNbElementsOfType ptr_getNbElementsOfType=0;
static fct_getSourceValue ptr_getSourceValue=0;
static fct_computeSourcesValues ptr_computeSourcesValues=0;

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
