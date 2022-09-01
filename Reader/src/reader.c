// file reader interface
// Implement a list based input file reader
//
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "list.h"
#include "memory.h"
#include "reader.h"
//
#define BUF_SIZE 1025
//
// message buffer used to 
// personalize error messages
static char msg[128];
//
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//
// Free the content of the reader
void ReaderFree(Reader *fr)
{
    List *line;
    while(fr->record.size)
    {
        line = ListPopBack(&(fr->record));
        ListFree(line,free);
        free(line);
    }
    fr->length = 0;
}
//
static char TrimSpace(FILE *fp)
{
    char ch;
    if(feof(fp)) return 0;
    while(1)
    {
        ch = getc(fp);
        if(ch!=' '&&ch!='\t') return ch;
    }
}
static char SaveString(char *buf,unsigned int *idx,FILE *fp,unsigned int linenum)
{
    char ch;
    if(feof(fp))
    {
        sprintf(msg,"Missing \' char to terminate string at line %u",linenum);
        ERROR_THROW(READER_ERROR,msg);
        ERROR_CATCH;
    }
    while(1)
    {
        ch = getc(fp);
        if(ch == '\''|| ch == '\"')
        {
            // found the string closing character
            buf[(*idx)++] = 0;
            break;
        }
        if(feof(fp)||ch=='\n')
        {
            sprintf(msg,"Missing \' char to terminate string at line %u",linenum);
            ERROR_THROW(READER_ERROR,msg);
            ERROR_CATCH;
        }
        buf[(*idx)++] = ch;
    }
    return getc(fp);
}
char SaveRecord(char *buf,unsigned int *idx,FILE *fp)
{
    char ch;
    while(1)
    {
        if(feof(fp)) return 0;
        ch = getc(fp);
        if(ch == '\t' || ch == ' ' || ch == '\n' || ch == '\'' || ch == '\"' || ch == '!') return ch;
        buf[(*idx)++] = ch;
    }
}
void JumpComment(FILE *fp)
{
    char ch;
    if(feof(fp)) return;
    while(1)
    {
        ch = getc(fp);
        if(feof(fp)||ch=='\n') return;
    }
}
//
//
// return the created list address
List *NewLine(Reader *reader,List *current_line)
{
    List *record;
    if(current_line)
    {
        // if nothing has been saved on last line return it
        if(current_line->size==0) return current_line;
    }
    record = SafeCalloc(1,sizeof(List));
    ListPushBack(&(reader->record),record);
    reader->length++;
    return record;
}
void NewRecord(Reader *reader,List *current_line,char *buf)
{
    char *str;
    // pushes a string only if has a length different than 0
    if(strlen(buf)==0) return;
    str = SafeCalloc(strlen(buf)+1,sizeof(char));
    strcpy(str,buf);// save string
    memset(buf,0,BUF_SIZE);//reset buffer
    ListPushBack(current_line,str);//push string in record list
}
//
// read and parse an input file
Reader *ReaderReadFile(const char *fname)
{
    Reader *reader;
    List *current_line;
    unsigned int line_number;
    FILE *fp;
    char c,*buf;
    unsigned int idx,space_flag,comment_flag,string_flag,begin_flag,record_flag;
    // allocate data and open text file
    fp = SafeFOpen(fname,"r");
    reader = SafeCalloc(1,sizeof(Reader));
    buf = SafeCalloc(BUF_SIZE,sizeof(char));
    // start reading
    // reset status flags
    idx = 0;
    string_flag = 0;
    comment_flag = 0;
    space_flag = 0;
    begin_flag = 1;
    record_flag = 0;
    line_number = 0;
    // pointer initialization
    current_line = 0;
    while(1) // read up to EOF
    {
        if(feof(fp)) break;
        if(begin_flag)
        {
            c = TrimSpace(fp); // trim initial space
            if(!current_line || current_line->size)// check exceptional cases
            {
                // if current line doesn't exists or is not empty
                current_line = NewLine(reader,current_line); // initialize new line
                line_number++;
            }
            idx = 0; // reset buffer index
            begin_flag = 0; // exit newbie state
            if(feof(fp)) break;
        }
        if(string_flag)
        {
            c=SaveString(buf,&idx,fp,line_number);
            buf[idx++]=0; // terminate string
            NewRecord(reader,current_line,buf); // save record
            idx = 0;
            string_flag = 0;
            if(feof(fp)) break;
        }
        if(comment_flag)
        {
            JumpComment(fp);
            c = '\n';
            comment_flag = 0;
            if(feof(fp)) break;
        }
        if(record_flag)
        {
            buf[idx++] = c;
            c = SaveRecord(buf,&idx,fp);
            buf[idx++] = 0;
            NewRecord(reader,current_line,buf);
            idx = 0;
            record_flag = 0;
            if(feof(fp)) break;
        }
        if(space_flag)
        {
            c = TrimSpace(fp); // trim initial space
            space_flag = 0;
            if(feof(fp)) break;
        }
        // set the flags to control the state flow
        switch(c)
        {
        case '\n':
            begin_flag = 1;
            break;
        case ' ':
        case '\t':
            space_flag = 1;
            break;
        case '#':
            comment_flag = 1;
            break;
        case '\'':
            string_flag = 1;
            break;
        case '\"':
            string_flag = 1;
            break;
        default :
            record_flag = 1;
            break;
        }
    }
    fclose(fp);
    free(buf);
    if(!current_line->size)
    {
        // check if last line is empty
        free(ListPopBack(&(reader->record)));
        reader->length--;
    }
    return reader;
}

void ReaderPrintInfo(FILE *stream, Reader *fr)
{
    List *line;
    unsigned int i0;
    fprintf(stream,"Number of lines read: %u\n",fr->length);
    for(i0 = 0; i0 < fr->length; i0++)
    {
        line = ListGetData(&(fr->record),i0);
        fprintf(stream,"Line %u: %u records saved\n",i0,line->size);
    }
}

void ReaderPrintFile(FILE *stream,Reader *fr)
{
    List *line;
    char *str;
    unsigned int i0,i1;
    for(i0 = 0; i0 < fr->length; i0++)
    {
        line = ListGetData(&(fr->record),i0);
        for(i1=0; i1<line->size; i1++)
        {
            str = ListGetData(line,i1);
            fprintf(stream,"'%s' ",str);
        }
        fprintf(stream,"\n");
    }
}


char *ReaderGetEntry(Reader *fr,unsigned int line,unsigned int index)
{
    List *l;
    char *record;
    l = ListGetData(&(fr->record),line);
    if(!l)
    {
        sprintf(msg,"Line %u does not exist.",line);
        ERROR_THROW(READER_ERROR,msg);
        ERROR_CATCH;
    }
    record = ListGetData(l,index);
    if(!record)
    {
        sprintf(msg,"Entry %u at line %u does not exist.",index,line);
        ERROR_THROW(READER_ERROR,msg);
        ERROR_CATCH;
    }
    return record;
}

unsigned int ReaderGetFileLength(Reader *fr)
{
    return fr->length;
}

unsigned int ReaderGetLineLength(Reader *fr,unsigned int line)
{
    List *l;
    l = ListGetData(&(fr->record),line);
    if(!l)
    {
        sprintf(msg,"Line %u does not exist.",line);
        ERROR_THROW(READER_ERROR,msg);
        ERROR_CATCH;
    }
    return l->size;
}
