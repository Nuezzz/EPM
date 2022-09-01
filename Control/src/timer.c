#include <stdio.h>
#include <time.h>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "timer.h"
//
// start the timer
void TimerStart(Timer *t)
{
    // prepare for the timing of the calculation
    #ifdef _OPENMP
    t->wtime1 = omp_get_wtime();
    #endif
    t->cstart = clock();
}
//    
// stop the timer
void TimerStop(Timer *t)
{
    #ifdef _OPENMP
    t->wtime2 = omp_get_wtime();
    t->tot_wtime += ((double)(t->wtime2 - t->wtime1));
    #endif
    t->cend = clock();
    t->ctime = ((double)(t->cend - t->cstart))/CLOCKS_PER_SEC;
    t->tot_ctime += t->ctime;
}
//
//
// report on the time lapse
void TimerReport(Timer *t,FILE *fp)
{
    //printf("MC simulation has finished! \n");
	int cputime_sec, cputime_min, cputime_hr;
	int walltime_sec, walltime_min, walltime_hr;

	cputime_sec = (int) t->ctime;
	cputime_hr = (int) (cputime_sec/3600);
	cputime_sec = cputime_sec - (cputime_hr*3600);
	cputime_min = (int) (cputime_sec/60);
	cputime_sec = cputime_sec - (cputime_min*60);

	walltime_sec = (int) t->wtime2-t->wtime1;
	walltime_hr = (int) (walltime_sec/3600);
	walltime_sec = walltime_sec - (walltime_hr*3600);
	walltime_min = (int) (walltime_sec/60);
	walltime_sec = walltime_sec - (walltime_min*60);

    printf("CPU time spent is: %d h, %d min, %d sec.\n",cputime_hr,cputime_min,cputime_sec);
    if (fp != NULL)
    {
        fprintf(fp,"CPU time spent is: %d h, %d min, %d sec.\n",cputime_hr,cputime_min,cputime_sec);
    }
    #ifdef _OPENMP
    printf("Wall time spent is: %d h, %d min, %d sec.\n",walltime_hr,walltime_min,walltime_sec);
    printf("Performance factor: %.3f\n",t->ctime/(t->wtime2-t->wtime1));
    printf("Number of threads deployed: %d\n\n",omp_get_max_threads());
    if (fp != NULL)
    {
        fprintf(fp, "Wall time spent is: %d h, %d min.\n",(int)floor(t->wtime2-t->wtime1)/3600,(int)ceil(fmod(t->wtime2-t->wtime1,3600.0)/60.0));
        fprintf(fp, "Performance factor: %.3f\n",t->ctime/(t->wtime2-t->wtime1));
        fprintf(fp, "Number of threads deployed: %d\n\n",omp_get_max_threads());
    }
    #endif
    fflush(stdout);
    if (fp != NULL) fflush(fp);
}
//
//
// Estimate ETA
#ifdef _OPENMP
void TimerETA(Timer *t,double current_percentage)
{
    double cwtime,eta;
    double remaining_perc;
    // get current time
    cwtime = (omp_get_wtime() - t->wtime1);
    eta = cwtime/(current_percentage);
    remaining_perc = 100.0 - current_percentage;
    eta *= remaining_perc;
    printf("ETA(h:m): %03d:%02d",(int)floor(eta)/3600,(int)ceil(fmod(eta,3600.0)/60.0));
    fflush(stdout);
}
#else
void TimerETA(Timer *t,double current_percentage)
{
    double cwtime,eta;
    double remaining_perc;
    // get current time
    cwtime = ((double)(clock() - t->cstart))/CLOCKS_PER_SEC;
    eta = cwtime/(current_percentage);
    remaining_perc = 100.0 - current_percentage;
    eta *= remaining_perc;
    printf("ETA(h:m): %d:%d",(int)floor(eta)/3600,(int)ceil(fmod(eta,3600.0)/60.0));
    fflush(stdout);
}
#endif
//
//
void TimerPrintProgress(Timer *t,double current_progress)
{
    if(current_progress == 0)
    {
        printf("\r[  0.0%%] - ETA(h:m): Inf");
    }
    else
    {
        printf("\r[%3.1lf%%] - ",current_progress);    
        TimerETA(t,current_progress);
    }
    fflush(stdout);
}
