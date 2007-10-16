#ifndef ACE_H
#define ACE_H

#include "extern.h"
#include <vector>
class unknown;
class equation;
class component;

#define ACE_INF 0.000000001
#define ACE_DOUBLE double

#define ACE_TYPE_NO 0
#define ACE_TYPE_RES 1
#define ACE_TYPE_CAP 2
#define ACE_TYPE_IND 3
#define ACE_TYPE_VSRC 4
#define ACE_TYPE_ISRC 5
#define ACE_TYPE_DIO 6
#define ACE_TYPE_COMP 7

#define ACE_TYPE_U 108
#define ACE_TYPE_V 109
#define ACE_TYPE_I 110

#define ACE_CHAR_LENGTH 24

typedef std::vector<unknown *> unknowns;
typedef std::vector<equation *> equations;
typedef std::vector<component *> components;



bool ACE_IS_NULL(ACE_DOUBLE d);
void ACE_ERROR(char * mess);
void ACE_WARNING(char * mess);
void ACE_INTERNAL_ERROR(char *mess);
void ACE_TYPE_TO_CHAR(int type,char* name);

#endif //ACE_H
