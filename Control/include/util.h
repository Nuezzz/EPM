#ifndef _UTIL_H
#define _UTIL_H
//
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
//
#define INPUT_BUFSIZE 1024
//
#define StringAlloc(s,nc) s = SafeCalloc(nc+1,sizeof(char))
#define StringClone(src,dst) src = SafeCalloc(strlen(dst)+1,sizeof(char)); strcpy(src,dst)
#define StringCat(s,a,b) s = SafeCalloc(strlen(a)+strlen(b)+1,sizeof(char)); strcpy(s,a); strcat(s,b)
#define StringFree(s) SafeFree(s)
//
// Assumes that at leas 1 char is in the string
#define StringUp(s,i0) i0=0;do{s[i0] = toupper(s[i0]);}while(s[++i0])
//
#define NextLine(fileptr,buffer,num_token,line) do{fgets(buffer,INPUT_BUFSIZE,fileptr);line++;num_token=StrTok(buffer,'#');}while(!num_token)
//
//
void ProgressBar(const char *msg,double perc);
char *FullPath(const char *folder,const char *filename);
void StrCat(char **dst,const char *str0);
// String uppercase
void StrToUpper(char *str);
// String tokenization
int StrTok(char *str,char comment);
char *GetToken(char *str,int i);
char *NextToken(char *str);
//
void StrExtTrim(char **string);

double Dot( double *a, double *b);
//
//
#endif
