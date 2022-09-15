#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "memory.h"
#include "util.h"
// File defining general utility functions
//
double Dot( double *a, double *b){
    double r=0;
    r = a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
    return r;
}
// Function to print a message and a progress bar on the screen
// It takes as input a double between 0 and 1
// message cannot be longer than 54 chars
void ProgressBar(const char *msg,double perc)
{
    char str[82];
    char per[5];
    int  nhashes;
    int  terminated;
    int i0,i1;
    // compute number of hashes
    nhashes = (perc*100.0+0.1)/5.0;
    sprintf(per,"%3d%%",((int)(perc*100.0+0.1)));
    str[0] = '\r';
    terminated = 0;
    i1 = 0;
    for(i0 = 1; i0 < 55; i0++)
    {
        terminated = !msg[i1];
        if(!terminated)
        {
            str[i0] = msg[i1++];
        }
        else
        {
            str[i0] = ' ';
        }
    }
    str[i0++] = '[';
    for(;i0 < 76; i0++)
    {
        if(nhashes)
        {
            str[i0] = '=';
            nhashes--;
        }
        else
        {
            str[i0] = ' ';
        }
    }
    str[i0++] = ']';
    str[i0++] = per[0];
    str[i0++] = per[1];
    str[i0++] = per[2];
    str[i0++] = per[3];
    str[i0] = 0;
    fputs(str,stdout);
    fflush(stdout);
}
//
//
char *FullPath(const char *folder,const char *filename)
{
    char *str;
    const int folder_len = strlen(folder);
    if(!folder_len)
    {
        str = SafeCalloc(strlen(filename)+1,sizeof(char));
        strcpy(str,filename);
        return str;
    }
    if(folder[folder_len-1] == '/')
    {
        StringCat(str,folder,filename);
    }
    else
    {
        str = SafeCalloc(folder_len+1+strlen(filename)+1,sizeof(char));
        strcpy(str,folder);
        strcat(str,"/");
        strcat(str,filename);
    }
    return str;
}
//
//
//
////////////////////////////////////////////////////////////////////////////////
// Concatenate two string and store the content in the first one
// It performs memory allocation and deallocation
////////////////////////////////////////////////////////////////////////////////
void StrCat(char **dst,const char *str0)
{
    if(!(*dst))
    {
        *dst = SafeCalloc(strlen(str0)+1,sizeof(char));
        strcpy(*dst,str0);
    }
    else
    {
        *dst = SafeRealloc(*dst,strlen(str0)+strlen((*dst))+1,sizeof(char));
        strcat(*dst,str0);
    }
}
//
//
////////////////////////////////////////////////////////////////////////////////
// Function that tokenize a string. It remove multiple spaces,
// tabs, and newline characters and separates token by a NULL character
// Return the number of tokens in the string
// Exclude any part of the string that has comments
// Comments starts with character 'comment'
////////////////////////////////////////////////////////////////////////////////
int StrTok(char *str,char comment)
{
    int i0 = 0;
    int i1;
    int num = 0;
    int set = 0;
    // Eliminate trailing chars
    i1 = 0;
    while(str[0] == ' ' || str[0] == '\t' || str[0] == '\n')
    {
        while(str[i1] && str[i1] != comment)
        {
            str[i1] = str[i1+1];
            i1++;
        }
        i1 = 0;
    }
    while(str[i0] && str[i0] != comment)
    {
        if(str[i0] != ' ' && str[i0] != '\t' && str[i0] != '\n')
        {
            if(!set)
            {
                set = 1;
                num++;
            }
        }
        else
        {
            i1 = i0;
            while(str[i1+1] == ' ' || str[i1+1] == '\t' || str[i1+1] == '\n')
            {
                while(str[i1] && str[i1] != comment)
                {
                    str[i1] = str[i1+1];
                    i1++;
                }
                i1 = i0;
            }
            set = 0;
            str[i0] = 0;
        }
        i0++;
    }
    str[i0] = 0;
    return num;
}
//
//
////////////////////////////////////////////////////////////////////////////////
// Get the next token in the string 
// the function doesn't check on the effective string length
////////////////////////////////////////////////////////////////////////////////
char *NextToken(char *str)
{
    int i0=0;
    while(str[i0++]) ;
    return &str[i0];
}
//
////////////////////////////////////////////////////////////////////////////////
// Extract the i-th token from a string
// there is no check on the string length 
////////////////////////////////////////////////////////////////////////////////
//
char *GetToken(char *str,int i)
{
    int i0=0;
    while(i--)
    {
        while(str[i0++]);
    }
    return &str[i0];
}
//
//
//
//
////////////////////////////////////////////////////////////////////////////////
// Convert a string to uppercase
////////////////////////////////////////////////////////////////////////////////
void StrToUpper(char *str)
{
    int i0=0;
    while(str[i0]){str[i0] = toupper(str[i0]); i0++;}
}
//
//
//
//
////////////////////////////////////////////////////////////////////////////////
// Remove extension including the point
////////////////////////////////////////////////////////////////////////////////
void StrExtTrim(char **string)
{
    char *str = *string;
    int i0=strlen(str)-1;
    while(i0 && str[i0]!='.') i0--;
    if(i0 == 0 && str[i0]!='.') return;
    *string = SafeRealloc(*string,i0+1,sizeof(char));
    (*string)[i0]=0;
}
// static inline int equal(const long n0[2],const long n1[2]);
// static inline int greater(const long n0[2],const long n1[2]);
// static inline int lesser(const long n0[2],const long n1[2]);
// static inline int leq(const Segment *t0,const Segment *t1);
//
static inline long select_pivot(double *a,const long l,const long r)
{
    const long m = (l+r)/2;
    const double lp = a[l];
    const double mp = a[m];
    const double rp = a[r];
    if((lp <= mp && mp <= rp) 
    || (rp <= mp && mp <= lp)) return m;
    if((mp <= rp && rp <= lp)
    || (lp <= rp && rp <= mp)) return r;
    if((rp <= lp && lp <= mp)
    || (mp <= lp && lp <= rp)) return l;
    // required to avoid warnings
    // but all cases have been covered by
    // the three previous if statements
    return m;
}
static inline long parallel_partition(double *a,double *less,double *more,long l,long r)
{
    double pivot;
    long i, jmin, jmax;
    long j;
    double t;
    // choose pivot and 
    // move it in the right buffer
    i = select_pivot(a,l,r);
    pivot = a[i];
    t = a[i];
    a[i] = a[r];
    a[r] = t;
    jmin = -1;
    jmax = -1;
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) private(j)
    #endif
    for(i = l; i < r; i++)
    {
        if(a[i] <= pivot)
        {
            #ifdef _OPENMP
            #pragma omp atomic capture
            #endif
            j = ++jmin;
            less[j] = a[i];
        }
        else
        {
            #ifdef _OPENMP
            #pragma omp atomic capture
            #endif
            j = ++jmax;
            more[j] = a[i];
        }
    }
    // the r-th element is the pivot
    // Its position is crucial so we
    // need to swap it outside the loop
    i = r;
    j = ++jmin;
    less[j] = a[i];
    // concatenate the arrays
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for(i = 0; i <= jmin; i++)
    {
        a[l+i] = less[i];
    }
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for(i = 0; i <= jmax; i++)
    {
        a[l+i+jmin+1] = more[i];
    }
    return jmin+l;
}
// Insertion sort
static inline void InsertionSort(double *a,const long l,const long r)
{
    long i,j;
    double t;
    for(i = l; i <= r; i++)
    {
        for(j = i-1; j >= l && (a[j] > a[j+1]); j--)
        {
            t = a[j+1];
            a[j+1] = a[j];
            a[j] = t;
        }
    }
}
//
//
void DoubleArrayQuickSort(double *a,const long l,const long r)
{
    long *min;
    long *max;
    long idx;
    long imin,imax,imid;
    double *less;
    double *more;
    //
    if(r-l+1 < 1) return ;
    min = SafeMalloc(r-l+1,sizeof(long));
    max = SafeMalloc(r-l+1,sizeof(long));
    less = SafeMalloc(r-l+1,sizeof(double));
    more = SafeMalloc(r-l+1,sizeof(double));
    idx = 0;
    min[idx] = l;
    max[idx] = r;
    idx++;
    while(idx)
    {
        // pop indices
        idx--;
        imin = min[idx];
        imax = max[idx];
        if(imax - imin + 1 <= 1e5)
        {
            InsertionSort(a,imin,imax);
        }
        else
        {
            imid = parallel_partition(a,less,more,imin,imax);
            // push indices
            min[idx] = imin;
            max[idx] = imid-1;
            idx++;
            min[idx] = imid+1;
            max[idx] = imax;
            idx++;
        }
    }
    free(min);
    free(max);
    free(less);
    free(more);
}
//
static inline int leq(const double t0,const double t1)
{
    if( t0 <= t1)
      return 1;
    else
      return 0;
}
