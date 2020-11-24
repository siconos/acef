#include "Spiceparser.h"
#include "perform.h"
#include <stdio.h>


static void* module=0;


void ParserPrintCircuit(){
  MEprint();
}

int ParserInitLibrary(){
  return 1;
}

void ParserStopLibrary(){
}
int ParserReadFile(char* file){
  return readFile(file);
}
int ParserGetNbElementsOfType(char *type){
  return getNbElementsOfType(type);
}

int ParserInitComponentList(char *type){
  return initComponentList(type);
}
int ParserNextComponent(void * data){
  return nextComponent(data);
}

int ParserGetSourceValue(char *type,void* id,double* value){
  return getSourceValue(type,id,value);
}
int ParserGetSourceValue2(int inttype,void* id,double* value){
  return getSourceValue2(inttype,id,value);
}

//int getISRCValue(void* id,double* value){
//  *value = ((ISRCinstance *) id)->currentValue;
//  return 1;
//}
//int getVSRCValue(void* id,double* value){
//  *value = ((VSRCinstance *) id)->currentValue;
//  return 1;
//}

int ParserComputeSourcesValues(double time){
  return computeSourcesValues(time);
}
int ParserGetTransValues(double * step, double * stop, double * start){
  return getTransValues(step,stop,start);
}
int ParserInitICvalue(){
  return initICvalue();

}
int ParserGetICvalue(int * numNode,int * icGiven, double * icValue) {
  return getICvalue(numNode,icGiven,icValue);
}
int ParserInitPrintElem(){
  return initPrintElem();

}
int ParserGetPrintElem(void ** p){
  return getPrintElem(p);
}
int ParserInitSimulation(int type,double val){
  return initSimulation(type,val);
}
void ParserSetIputFromId(char* type,unsigned int id,double v){
  if (id)
    setIputFromId(type, id, v);
  
}
unsigned int ParserGetId(char* type,char* name){
  return getId( type, name);
}
