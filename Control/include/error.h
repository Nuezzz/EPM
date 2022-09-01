//
// ERROR HANDLING MODULE
// Errors are efficiently handled
// and saved onto an output file for
// further inspection
//
#ifndef _ERROR_H
#define _ERROR_H
//
//
//
// Header file for handling errors
//
#define ERROR_THROW(code,message) ErrorThrow(code,message,__FILE__,__FUNCTION__,__LINE__) 
#define ERROR_CATCH ErrorCatch()
#define ERROR_CHECK ErrorCheck()
//
// define the code for memory allocation error
#define ALLOC_ERROR     1
#define FOPEN_ERROR     2
#define MESH_ERROR      3
#define READER_ERROR    4
#define INPUT_ERROR     5
#define OUTPUT_ERROR    6
#define SIM_ERROR       7
#define UMFPACK_ERROR   8
//
//
//
// to call at the beginning of main
void ErrorStreamOpen(const char *fname);
// return 0 if no error was thrown
void ErrorThrow(const int code,const char *message,const char *file,const char *func,const int line);
// check
int ErrorCheck();
// function to catch error
void ErrorCatch();
// to call before exiting the main()
void ErrorStreamClose();
//
#endif // __ERROR_H__
