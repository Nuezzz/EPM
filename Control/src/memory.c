#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/stat.h>
#include <errno.h>
#include "error.h"
//
//
// PARAMETERS TO ALLOCATION SAFE FUNCTIONS
// - n: is the number of elements
// - size: is the size of each element
//
// safe malloc
void *SafeMalloc(const unsigned long int n,const unsigned int size)
{
    void *ptr = 0;
    const void *null_ptr = 0;
    ptr = malloc(size*n);
    assert(ptr != null_ptr);
    return ptr;
}
// safe calloc
void *SafeCalloc(const unsigned long int n,const unsigned int size)
{
    void *ptr = 0;
    const void *null_ptr = 0;
    ptr = calloc(n,size);
    assert(ptr != null_ptr);
    return ptr;
}
// safe realloc
// For safety reasons use this function to reduce the effective size of
// a vector, not to expand it.
void *SafeRealloc(void *ptr,const unsigned long int n,const unsigned int size)
{
    const void *null_ptr = 0;
    ptr = realloc(ptr,size*n);
    assert(ptr != null_ptr);
    return ptr;
}
// safe fopen
FILE *SafeFOpen(const char *filename,const char *mode)
{
    char err_str[256];
    FILE *file_ptr;
    const FILE *null = 0;
    file_ptr = fopen(filename,mode);
    if(file_ptr == null)
    {
        sprintf(err_str, "Error opening file %s.\n",filename);
        printf("Error opening file %s.\n",filename);
        fflush(stdout);
        ERROR_THROW(INPUT_ERROR, err_str);
        ERROR_CATCH;
    }
    //assert(file_ptr != null);
    return file_ptr;
}
//
//
void CreateFolder(const char *folder)
{
    if(mkdir(folder, S_IRWXU) == -1)
    {
        switch(errno)
        {
            case EEXIST:
                printf("Directory %s already exists.\n",folder);
                fflush(stdout);
                break;
            case EACCES:
                printf("Error: Cannot create directory %s. Permission denied.\n",folder);
                fflush(stdout);
                abort();
            case ENOSPC:
                printf("Error: Cannot create directory %s. Not enough space on disk.\n",folder);
                fflush(stdout);
                abort();
            default:
                printf("Error: Cannot create directory %s. Obscure forces prevented the creation of the directory.\n",folder);
                fflush(stdout);
                abort();
        }
    }
}
