// MODULE FOR SAFE MEMORY ALLOCATION
// The allocation functions of this module
// will legally
//
#ifndef _MEMORY_H
#define _MEMORY_H
//
#include <stdio.h>
#include <unistd.h>
//
#define SafeFree(p) if(p)free(p)
#define FileExists(f) (access(f, F_OK) != -1)
//
// safe malloc
void *SafeMalloc(const unsigned long int n,const unsigned int size);
void *SafeRealloc(void *ptr,const unsigned long int n,const unsigned int size);
void *SafeCalloc(const unsigned long int n,const unsigned int size);
FILE *SafeFOpen(const char *filename,const char *mode);
void CreateFolder(const char *folder);
//
#endif // __MEMORY_H
