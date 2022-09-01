#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#ifdef _OPENMP
#include <omp.h>
#endif // _OPENMP
//
// default error file name
// - if the user doesn't supply one
//   the default value will be used
//
#define DEFAULT_ERR_MESSAGE "Error detected."
//
// Static variables to handle errors
// only functions from this module can exit the
// program through an exit() function
//
//
static FILE *error_stream=0;
static char error_file_name[64];
static int  error_code; // last error code
//
//
//
// Initialize static variables
// to be called at the beginning of main()
void ErrorStreamOpen(const char *fname)
{
    const FILE *null = 0;
    #ifdef _OPENMP
    #pragma omp critical (errorset)
    #endif
    {
        if(fname)
        {
            error_stream = fopen(fname,"w");
            strcpy(error_file_name,fname);
            //printf("Error stream redirected to %s.\n",fname);
            //fflush(stdout);
        }
        else
        {
            error_stream = stderr;
            strcpy(error_file_name,"standard error stream");
            //printf("Using standard error stream.\n");
            //fflush(stdout);
        }
        assert(error_stream != null);
        error_code = 0;
    }
}
//
//
void ErrorThrow(const int code,const char *message,const char *file,const char *func,const int line)
{
    #ifdef _OPENMP
    #pragma omp critical (errorset)
    #endif
    {
        fprintf(error_stream,"Error %d in %s, %s() at %d: %s\n",code,file,func,line,(message[0])? message : DEFAULT_ERR_MESSAGE);
        fflush(error_stream);
        error_code = code;
    }
}
//
//
//
// Catching an error means taking action
//  - The action is report the error code and abort execution
// flush the error stream and exit with current error_code
void ErrorCatch()
{
    #ifdef _OPENMP
    #pragma omp critical (errorset)
    #endif
    {
        printf("\nError flag raised.\nProgram execution stopped.\nSee %s for details.\n",error_file_name);
        fflush(stdout);
        fprintf(error_stream,"Exit code: %d\n",error_code);
        fflush(error_stream);
        exit(error_code);
    }
}
// check if an error has been thrown
int ErrorCheck()
{
    int code;
    #ifdef _OPENMP
    #pragma omp atomic read
    #endif
    code = error_code;
    return code;
}
//
void ErrorStreamClose()
{
    #ifdef _OPENMP
    #pragma omp critical (errorset)
    #endif
    {
        if(error_stream)
        {
            if(!error_code)
                fprintf(error_stream,"No errors detected.\n");
            fflush(error_stream);
            if(error_stream != stderr)
                fclose(error_stream);
        }
    }
}
