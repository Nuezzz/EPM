#ifndef _TIMER_H
#define _TIMER_H
//
#include <time.h>
//
// if openmp is not active wall clock
// and cpu clocks are equal
typedef struct timer
{
    clock_t cstart, cend;   // cpu clocks
    double  wtime1, wtime2; // wall clocks
    double  ctime;          // time interval
    double  tot_wtime;      // total wall time
    double  tot_ctime;      // total cpu time

} Timer;
//
//
void TimerStart(Timer *t);
void TimerStop(Timer *t);
void TimerReport(Timer *t,FILE *fp);
void TimerETA(Timer *t, double current_perc);
void TimerPrintProgress(Timer *t,double current_progress);
#endif
