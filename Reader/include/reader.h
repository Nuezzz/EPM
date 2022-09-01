#ifndef _READER_H
#define _READER_H
//
//
#include "list.h"
#include <stdlib.h>
#include <stdio.h>
//
// in this data structure
// a file object is represented by a list
//  - each entry of a list is another list representing a line
//  - a line is made by string separated by spaces
//  - multiple spaces are ignored
//  - the character # makes a comment
//  - empty lines are not saved
//  - if you wish to include the ' ' char in an entry just enclose
//    the record in '...'
typedef struct reader
{
    unsigned int    length;     // number of lines
    struct list     record;     // lines
} Reader;
//
Reader          *ReaderReadFile     (const char *fname);
char            *ReaderGetEntry     (Reader *fr,unsigned int row,unsigned int column);
unsigned int    ReaderGetFileLength (Reader *fr);
unsigned int    ReaderGetLineLength (Reader *fr,unsigned int line);
void            ReaderFree          (Reader *fr);
void            ReaderPrintFile     (FILE *stream,Reader *fr);
// Utility functions
int ReaderFindKeyword(Reader *fr,const char *key);
//
#endif
