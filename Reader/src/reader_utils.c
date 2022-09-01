// This module defines top level functions for the parser
// like localization of specific strings accross a read file
// or parsing of integers and floating point values 
#include "reader.h"
#include "error.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
//
//
//
//
// Look into the read file for the line
// having first word equal to the keyword
// returns -1 if there is no such line beginning
// with the given string
int ReaderFindKeyword(Reader *fr,const char *key)
{
    int i0;
    char *s;
    for(i0 = 0; i0 < ReaderGetFileLength(fr); i0++)
    {
        s = ReaderGetEntry(fr,i0,0);
        if(!strcmp(s,key))
        {
            return i0;
        }
    }
    return -1;
}
//

