#ifndef ACE_H
#define ACE_H

#include "extern.h"
#include <vector>

class unknown;
class equation;
class component;

#define ACE_INF 0.00000000000000001
#define ACE_NULL 0.000001
#define ACE_DOUBLE double



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

typedef std::vector<unknown *> unknowns;
typedef std::vector<equation *> equations;
typedef std::vector<component *> components;

typedef std::vector<component *>::iterator componentsIt;

extern char ACE_name[];

bool ACE_IS_NULL(ACE_DOUBLE d);
void ACE_MESSAGE(char * mess);
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
#endif //ACE_H
