

#define PARSER_TSTEP 1
#define PARSER_TSTOP 2


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
