#ifndef ACE_H
#define ACE_H

#include "extern.h"
#include <vector>

class unknown;
class equation;
class component;
#include "SiconosAlgebra.h"

#define ACE_INF 0.00000000000000001
#define ACE_NULL_COEF_MAT 0.00000000000000001
#define ACE_NULL 0.000001
#define ACE_DOUBLE double
#define ACE_BIG_FLOAT 1e9



#define ACE_TYPE_NO 0
#define ACE_TYPE_RES 1
#define ACE_TYPE_CAP 2
#define ACE_TYPE_IND 3
#define ACE_TYPE_VSRC 4
#define ACE_TYPE_ISRC 5
#define ACE_TYPE_DIO 6
#define ACE_TYPE_COMP 7
#define ACE_TYPE_MOS 8
#define ACE_TYPE_VCVS 9
#define ACE_TYPE_ARB 10
#define ACE_TYPE_VCCS 11



#define ACE_TYPE_U 108
#define ACE_TYPE_V 109
#define ACE_TYPE_I 110

#define ACE_CHAR_LENGTH 124

//SOLVER TYPE
#define ACE_SOLVER_ENUM 0
#define ACE_SOLVER_SIMPLEX 1
#define ACE_SOLVER_PATH 2
#include "acetime.h"


#define ACE_TIMER_MAIN 0
#define ACE_TIMER_EQUATION 1
#define ACE_TIMER_SIMULATION 2
#define ACE_TIMER_SOLVER 3
#define ACE_TIMER_LS_STEP 4
#define ACE_TIMER_SOLVE_ENUM 5
#define ACE_TIMER_SOLVE_SIMPLEX 6
#define ACE_TIMER_SOLVE_PATH 7
#define ACE_TIMER_SOLVE_GUESS 8
#define ACE_TIMER_SOLVE_LU 9
#define ACE_TIMER_DIRECT 10
#define ACE_TIMER_LU_DIRECT 11
#define ACE_TIMER_SIMPLEX_FIRST 12
#define ACE_TIMER_SIMPLEX_GUESS 13
#define ACE_TIMER_SIMPLEX_TREE 14
#define ACE_TIMER_SIMPLEX_TRY_NODE 15
#define ACE_TIMER_COMPUTE_VAR 16
#define ACE_TIMER_PROD_MAT 17
#define ACE_TIMER_TEST 18
#define ACE_TIMER_TEST_1 19
#define ACE_TIMER_TEST_2 20
#define ACE_TIMER_TEST_3 21
#define ACE_TIMER_LAST 22

typedef std::vector<unknown *> unknowns;
typedef std::vector<equation *> equations;
typedef std::vector<component *> components;

typedef std::vector<component *>::iterator componentsIt;



extern char ACE_name[];
extern aceTime ACE_times[];
extern int ACE_SOLVER_TYPE;
extern UBLAS_TYPE ACE_MAT_TYPE;
extern int ACE_MUET_LEVEL; //0 verbose .... 10 muet
#define ACE_MUET 10


void ACE_INIT();
bool ACE_IS_NULL(ACE_DOUBLE d);
void ACE_MESSAGE(char * mess,int level=0);
void ACE_ERROR(char * mess);
void ACE_WARNING(char * mess);
void ACE_INTERNAL_ERROR(char *mess);
void ACE_INTERNAL_WARNING(char *mess);
void ACE_TYPE_TO_CHAR(int type,char* name);

//Raise an error if not b
//check b is true
void ACE_CHECK_IERROR(bool b,char* mess);
void ACE_CHECK_IWARNING(bool b,char* mess);
void ACE_CHECK_WARNING(bool b,char* mess);
void ACE_CHECK_ERROR(bool b,char* mess);


//TIME FUNCTION
void ACE_INIT_TIME();
void ACE_PRINT_TIME();
void ACE_STOP_SOLVER_TIME();
#endif //ACE_H
