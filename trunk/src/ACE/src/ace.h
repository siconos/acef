/*
 *
 */
/**
 *\file ace.h
  \brief This file contains declations of constants and tools used in the ACEF module.


*/

#ifndef ACE_H
#define ACE_H

/**
 *
 * define PRE_COMPUTE_ADAPTIVE to use the precompute adaptive time stepping.
 *
 *
 */

//#define PRE_COMPUTE_ADAPTIVE

#include "extern.h"
#include <vector>
#include <fstream>


class unknown;
class equation;
class component_NONE_LINEAR_NS;
class component;
class node;
#include "SiconosAlgebra.hpp"

#define ACE_INF 1e-17   
#define ACE_NULL_COEF_MAT 1e-17 
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
#define ACE_TYPE_BJT 12
#define ACE_TYPE_RELAY 13
#define ACE_TYPE_MOS_NL 14
#define ACE_TYPE_MOS_NL2 15




#define ACE_TYPE_U 108
#define ACE_TYPE_V 109
#define ACE_TYPE_I 110

#define ACE_CHAR_LENGTH 124

//X'=....
#define ACE_FORMULATION_SEMI_EXPLICT 0
//MX'= . with cap current unknowns
#define ACE_FORMULATION_WITHOUT_INVERT 1
//MNA
#define ACE_FORMULATION_MNA 2
//MNA_V
#define ACE_FORMULATION_MNA_V 3
//MNA_V
#define ACE_FORMULATION_STAMP_ONLY 4

//Semi-explicit index 1(?)
#define ACE_FORMULATION_SE1 5

//DAE, SEMI-EXPLICT, MNA
extern int ACE_FORMULATION;
//Adaptatve time stepping or SEMI-EXPLICT
extern int ACE_WITH_ADAPTATIVE_TIME_STEPPING;//0 or 1
extern ACE_DOUBLE ACE_ATOL_LOCAL;
extern ACE_DOUBLE ACE_RTOL_LOCAL;

/*options discretization*/
extern double ACE_THETA_X;
extern double ACE_THETA_ZS;
extern double ACE_THETA_P;

/*option for the MOS component*/
extern int ACE_USE_NL_MOS;
extern int ACE_USE_SMOOTH_MOS;
extern int  ACE_MOS_NB_HYP;
extern ACE_DOUBLE ACE_MOS_POWER_SUPPLY;

/*diode options*/
extern ACE_DOUBLE ACE_DIODE_THRESHOLD;

#ifdef PRE_COMPUTE_ADAPTIVE
#define ACE_NB_ADAPT_STEP 7
#else
#define ACE_NB_ADAPT_STEP 0
#endif

extern int ACE_CMP_ADAT[ACE_NB_ADAPT_STEP+1];


/*

 *ACE_CUR_STEP means  ==> mH* 2 pow ACE_CUR_STEP
 */
extern int ACE_CUR_STEP;

//SOLVER TYPE
#define ACE_SOLVER_ENUM 0
#define ACE_SOLVER_SIMPLEX 1
#define ACE_SOLVER_PATH 2
#define ACE_SOLVER_NUMERICS_DIRECT_ENUM 3
#define ACE_SOLVER_NUMERICS_DIRECT_SIMPLEX 4
#define ACE_SOLVER_NUMERICS_DIRECT_PATH 5
#define ACE_SOLVER_FB 6
#define ACE_SOLVER_NUMERICS_PATH_ENUM 7
#define ACE_SOLVER_NUMERICS_DIRECT_PATH_ENUM 8

#include "acetime.h"


#define ACE_TIMER_MAIN 0
#define ACE_TIMER_EQUATION 1
#define ACE_TIMER_SIMULATION 2
#define ACE_TIMER_LS_STEP 3
#define ACE_TIMER_SOLVE_NUMERICS 4
#define ACE_TIMER_COMPUTE_VAR 5
#define ACE_TIMER_PROD_MAT 6
#define ACE_TIMER_TEST 7
#define ACE_TIMER_TEST_1 8
#define ACE_TIMER_TEST_2 9
#define ACE_TIMER_TEST_3 10
#define ACE_TIMER_TEST_4 11
#define ACE_TIMER_TEST_5 12
#define ACE_TIMER_TEST_6 13
#define ACE_TIMER_TEST_7 14
#define ACE_TIMER_TEST_8 15
#define ACE_TIMER_TEST_9 16
#define ACE_TIMER_TEST_10 17
#define ACE_TIMER_TEST_11 18
#define ACE_TIMER_LAST 19

typedef std::vector<unknown *> unknowns;
typedef std::vector<equation *> equations;
typedef std::vector<node *> nodes;
typedef std::vector<component_NONE_LINEAR_NS *> component_NLs;
typedef std::vector<component *> components;

typedef std::vector<component *>::iterator componentsIt;



extern char ACE_name[];
extern aceTime ACE_times[];
extern int ACE_SOLVER_TYPE;
extern UBLAS_TYPE ACE_MAT_TYPE;
extern int ACE_MUET_LEVEL; //0 verbose .... 10 muet
#define ACE_MUET 10

/**
 * 
 * \brief Must be call at the begining.
 */

void ACE_INIT();
void ACE_STOP();
bool ACE_IS_NULL(ACE_DOUBLE d);
void ACE_MESSAGE(char * mess,int level=0);
void ACE_ERROR(char * mess);
void ACE_WARNING(char * mess);
void ACE_INTERNAL_ERROR(char *mess);
void ACE_INTERNAL_WARNING(char *mess);

/**
 * 
 * \fn void ACE_TYPE_TO_CHAR(int type,char* name)
 * \brief convert an integer to char
 * \param[in] type is an integer.
 * \param[out] the name of the type.
 *
 */
void ACE_TYPE_TO_CHAR(int type,char* name);

ofstream & ACE_GET_LOG_STREAM();
ofstream & ACE_GET_LOG1_STREAM();

/**
 * 
 * \brief Check b is true, raise an error if not b
 */
void ACE_CHECK_IERROR(bool b,char* mess);
void ACE_CHECK_IWARNING(bool b,char* mess);
void ACE_CHECK_WARNING(bool b,char* mess);
void ACE_CHECK_ERROR(bool b,char* mess);


//TIME FUNCTION
void ACE_INIT_TIME();
void ACE_PRINT_TIME();
void ACE_STOP_SOLVER_TIME();

/*
 *Functions to set input
 *
 */
unsigned int ACE_GET_ID(char* type,char* name);
void ACE_SET_INPUT(char* type,unsigned int id, double v);
#endif //ACE_H
