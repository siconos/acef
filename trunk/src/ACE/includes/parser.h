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
  
int readFile(char *file);
int initParserLibrary();
void stopParserLibrary();
int initComponentList(char *type);
int nextComponent(void * data);
int getNbElementsOfType(char *type);
void printCircuit();
  //input : component's Id
//output : value.
//if succes, return 1 else return 0 
int getSourceValue(char *type,void* id,double* value);
int computeSourcesValues(double time);

  /*Get .trans values*/
  int getTransValues(double * step, double * stop, double * start);
  /*Get initial values .IC*/
  int initICvalue();  
  int getICvalue(int * numNode,int * icGiven, double * icValue);

  /*.print tran management*/
  int initPrintElem();
  /* p must be casted in dataPrint ** */
  int getPrintElem(void ** p);
  
  int initSimulation(int type,double val);

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
typedef struct {
  char *name;
  int nodePos;
  int nodeNeg;
  int nodeS;
  double epsilon;
  double vplus;
  double vmoins;
  int type;
} dataCOMP;
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
#ifdef __cplusplus
}
#endif
