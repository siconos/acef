#ifndef __SICONOS_PERFORM_H
#define __SICONOS_PERFORM_H



int readFile(char *file);
void initParserLibrary();
void MEperform();
void MEprint();
int initComponentList(char *type);
int nextComponent(void * data);
int getNbElementsOfType(char* type);
//input : component's Id
//output : value.
//if succes, return 1 else return 0 
int getSourceValue(char *type,void* id,double* value);
int computeSourcesValues(double time);
int getTransValues(double * step, double * stop, double * start);
int getICvalue(int * numNode,int * icGiven, double * icValue);
int initICvalue();

void setPrintStr(char * str);
int myCKTNodeId(char * name);
void freePrintElem();
int initPrintElem();
int getPrintElem(void ** p);

int initSimulation(int type,double val);

void setIputFromId(char* type,unsigned int id,double v);
unsigned int getId(char* type,char* name);


#endif
