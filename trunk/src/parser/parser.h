#ifdef __cplusplus
extern "C" {
#endif

int readFile(char *file);
int initParserLibrary();
void stopParserLibrary();
int initComponentList(char *type);
int nextComponent(void * data);
int getNbElementsOfType(char *type);
void printCircuit();

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

/* structures to describe Voltage Controlled Current Sources */
typedef struct {
  char *name;
  int nodePos;
  int nodeNeg;
  double *coef;
  double value;
} dataVCCS;


/*structures to describe independent voltage sources*/
typedef struct {
  char *name;
  int nodePos;
  int nodeNeg;
  double value;
} dataVSRC;


/* information needed per source instance */
typedef struct {
  char *name;
  int nodePos;
  int nodeNeg;
  double value;
  double *coef;
} dataISRC;
#ifdef __cplusplus
}
#endif
