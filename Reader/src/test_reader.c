#include <stdio.h>
#include <stdlib.h>
#include "error.h"
#include "memory.h"
#include "openmp.h"
#include "timer.h"
#include "version.h"
#include "reader.h"


int main(int argc,char **argv)
{
    Timer *perfbm;
    Reader *fr;
    //
    ErrorStreamOpen(0);
    perfbm = SafeCalloc(1,sizeof(Timer));
    if(ErrorCheck())
    {
        ErrorCatch();
    }
    TimerStart(perfbm);
    printf("MC3D v.%s\n",STR(VERSION));
    printf("Test run for the Reader module.\n");
    printf("Input file to parse: %s\n",argv[1]);
    fr = ReaderReadFile(argv[1]);
    ReaderPrintFile(stdout,fr);
    ReaderFree(fr);
    TimerStop(perfbm);
    TimerReport(perfbm);
    ErrorStreamClose();
    free(fr);
    free(perfbm);
    return 0;
}
