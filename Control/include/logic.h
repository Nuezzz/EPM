#ifndef _LOGIC_H
#define _LOGIC_H
//
// Definition of logical type
typedef int bool;
enum{false = 0,true = 1}; 
// MACRO
// define logic xor
#define xor(a,b) ((!(a)&&(b))||((a)&&!(b)))
//
//
#endif
