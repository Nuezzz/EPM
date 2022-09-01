#include "openmp.h"
#ifdef _OPENMP
#include <omp.h>
#endif // _OPENMP
// define global variable
// that set the number of
// parallel threads
static unsigned int omp_thread_num = 1;
// single thread by default
// this function has to be called before any other
// parallelization routine
void OMPSetThreadNum(unsigned int n)
{
    #ifdef _OPENMP
    if(n <= 0 || n > omp_get_thread_limit())
    {
        omp_thread_num = 1;
    }
    else
        omp_thread_num = n;
    omp_set_dynamic(0);     // Explicitly disable dynamic teams
    omp_set_num_threads(omp_thread_num); // Use omp_thread_num threads for all consecutive parallel regions
    #endif // _OPENMP
}
//
//
// set dynamic threading
void OMPSetDynamicThreading()
{
    #ifdef _OPENMP
    omp_set_dynamic(1);     // Explicitly enable dynamic teams
    #endif // _OPENMP
}
//
unsigned int OMPGetThreadNum()
{
    return omp_thread_num;
}
