#ifndef __SPICE_PARSER_API_H
#define __SPICE_PARSER_API_H



#ifdef __cplusplus
extern "C" {
#endif
#define PARSER_TSTEP 1
#define PARSER_TSTOP 2


  typedef struct {
  int node1;
  int node2;
  void * next;
  char * name;
} dataPrint;
  
int ParserReadFile(char *file);
int ParserInitLibrary();
void ParserStopLibrary();
int ParserInitComponentList(char *type);
int ParserNextComponent(void * data);
int ParserGetNbElementsOfType(char *type);
void ParserPrintCircuit();
  //input : component's Id
//output : value.
//if succes, return 1 else return 0 
int ParserGetSourceValue(char *type,void* id,double* value);
int ParserComputeSourcesValues(double time);

  /*Get .trans values*/
  int ParserGetTransValues(double * step, double * stop, double * start);
  /*Get initial values .IC*/
  int ParserInitICvalue();  
  int ParserGetICvalue(int * numNode,int * icGiven, double * icValue);

  /*.print tran management*/
  int ParserInitPrintElem();
  /* p must be casted in dataPrint ** */
  int ParserGetPrintElem(void ** p);
  
  int ParserInitSimulation(int type,double val);

  void ParserSetIputFromId(char* type,unsigned int id,double v);
  unsigned int ParserGetId(char* type,char* name);


typedef struct {
  char *name;
  int nodePos;
  int nodeNeg;
  double value;
  
} dataRES;

typedef struct {
  char *name;
  int nodePos;
  int nodeNeg;
  double value;
  
} dataCAP;

typedef struct {
  char *name;
  int nodePos;
  int nodeNeg;
  double value;
  
} dataIND;

typedef struct {
  char *name;
  int nodePos;
  int nodeNeg;
  double area;
} dataDIO;
typedef struct {
  char *name;
  int mode;
  /*int nodePos;
    int nodeNeg;*/
  int drain;
  int gate;
  int source;
  double w;
  double k;
  double vt;
} dataMOS1;
typedef struct {
  char *name;
  int mode;
  int collector;
  int base;
  int emitor;
} dataBJT;
/* structures to describe Voltage Controlled Current Sources */
typedef struct {
  char *name;
  int nodePos;
  int nodeNeg;
  double coef;
  int nodeDriverPos;
  int nodeDriverNeg;
  double value;
} dataVCCS;
/* structures to describe Voltage Controlled Voltage Sources */
typedef struct {
  char *name;
  int nodePos;
  int nodeNeg;
  double coef;
  int nodeDriverPos;
  int nodeDriverNeg;
  double value;
} dataVCVS;

typedef struct {
  char *name;
  int nodePos;
  int nodeNeg;
  int type;
} dataARB;

/*structures to describe independent voltage sources*/
typedef struct {
  char *name;
  int nodePos;
  int nodeNeg;
  double value;
  void* id;
} dataVSRC;


/* information needed per source instance */
typedef struct {
  char *name;
  int nodePos;
  int nodeNeg;
  double value;
  double *coef;
  void* id;
} dataISRC;

  /* information needed per comparator */
  typedef struct {
    char *name;
    int nodePos;
    int nodeNeg;
    int nodeOut;
    double vmin;
    double vmax;
    double vepsilon;
    double voffset;
  } dataCOMP;

#ifdef __cplusplus
}
#endif

#endif
