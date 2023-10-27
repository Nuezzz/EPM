#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "openmp.h"
#include "atom.h"
#include "error.h"
#include "memory.h"
#include "util.h"
#include "ham.h"
#include <complex.h>

#define  KMAX  40
//#define  STEPS 20

/**
 * @brief Found the G_vec sets corresponding to specific K_vec within the E_cut ,
 * build the Hamitonian, solve it's eigen value and print out the band
 * @param band 
 * @param L 
 * @param Estack 
 * @param Hstack 
 * @param k 
 * @param E_cut 
 * @param n_threads 
 */
static inline void KPatch(FILE *band, Lattice *L, double **k, double E_cut, int n_threads )
  {
    
    Eigen* Estack[n_threads];
    double complex* Hstack[n_threads];
    int n_band = L->a_set->n_atoms*2+8;
    int i,thread;
    int m[n_threads];

    TimerStart(&L->Form_time);
    #ifdef _OPENMP
    #pragma omp parallel private(thread)
    #endif
    {
        #ifdef _OPENMP
            thread = omp_get_thread_num();
        #else
            thread = 0;
        #endif
            Estack[thread]=GVecInit( L,  KMAX,k[thread], E_cut);
    }
    printf("Finish Building the G_vec mesh with \n");
    fflush(stdout);


    for(i=0;i<n_threads;i++)
    {
        Hstack[i]=HTot_loc(L,Estack[i],k[i]);
        printf("Finish Building the %d th Hamitonian with %d rank\n", i+1, Estack[i]->NG);
        fflush(stdout); 
    }
    TimerStop(&L->Form_time);
       
    

        // #ifdef _OPENMP
        // #pragma omp parallel private(thread)
        // #endif
        // {
        //     #ifdef _OPENMP
        //     thread = omp_get_thread_num();
        //     #else
        //     thread = 0;
        //     #endif
        //     m[thread]=CalcBand( Estack[thread], Hstack[thread], n_band);
        // }
        
        //print_rmatrix( "Selected eigenvalues", 1, m, L->E, 1 );
        TimerStart(&L->Solve_time);
        for(i=0;i<n_threads;i++)
        {
            m[i]=CalcBand( Estack[i], Hstack[i], n_band);
        }
        TimerStop(&L->Solve_time);

        for(i=0;i<n_threads;i++)
        {
            PrintEigen(band, Estack[i]->E, k[i],m[i]);
            BandFinish(Estack[i]);

        }


  }

static inline double *KBuild(Lattice *L, int N)
{
    double unit_a   = 2.00*PI/L->A0;
    double k_0[3];
    double k_1[3];
    double dk[3];
    double *k_vector =SafeCalloc(3*N,sizeof(double));

    k_0[0]=0;                   k_0[1]=0;                   k_0[2]=0;
    k_1[0]=unit_a/2.0;          k_1[1]=unit_a/2.0;          k_1[2]=unit_a/2.0;
    dk [0]=(k_1[0]-k_0[0])/N;   dk [1]=(k_1[1]-k_0[1])/N;   dk [2]=(k_1[2]-k_0[2])/N;
    for(int i=0; i<N; i++)
    {
        k_vector[i*3+0]= k_0[0] + i*dk[0];
        k_vector[i*3+1]= k_0[1] + i*dk[1];
        k_vector[i*3+2]= k_0[2] + i*dk[2];
    }
    return k_vector;
    
}


//create G_vec array for each K_vec in the first patch, each patch should 
//contain n_treads of K_vec
int main(int argc, char **argv)
{
    FILE *band;
    Lattice *L;

    Timer   tot_time;

    int     i,j;
    char   *foldername;
    char   *simname = argv[1];
    int    n_threads = atoi(argv[2]);
    int    STEPS=atoi(argv[3]);
    
    StringClone(foldername,simname);
    StrCat(&foldername, "_band_structure");
    

    // double complex *H_tot;
    
    double  E_cut = 80.0;
    double  *k[n_threads];
    double  *k_path;

       


	printf("Empirical Sudopotential Calculation\n ");
    fflush(stdout);
    ErrorStreamOpen("error.log");



    TimerStart(&tot_time);
    if(argc == 4)
    {
        CreateFolder(foldername);
        band = OpenBandFile(foldername);
        OMPSetThreadNum(n_threads);


        L = LatticeInitial(simname);

        L->Form_time.tot_wtime=0;
        L->Form_time.tot_ctime=0;
        L->Solve_time.tot_wtime=0;
        L->Solve_time.tot_ctime=0;

        k_path=KBuild(L, STEPS*n_threads);
        //PPtest(L,foldername);

        for(i=0;i<STEPS;i++)
        {
            for(j=0;j<n_threads;j++)
            {
                k[j]= k_path+3*(i*n_threads+j);
            }
            KPatch(band, L,k, E_cut, n_threads );
            printf("The %d patch finished, %f %% \n ",i+1, (float)(i+1)/(float)(STEPS)*100);
            fflush(stdout);
        }

        printf("The time spent to biuld the Hamitonian: %f\n",L->Form_time.tot_wtime);
        printf("The time spent to solve the eigen problem: %f\n",L->Solve_time.tot_wtime);
    
        fclose(band);
    }

    else
    {
        printf("No tittle was selected\n");
        return 1;
    }
    
    TimerStop(&tot_time);
    printf("Total time spend for %d atoms, and %d k points:\n",L->a_set->n_atoms,STEPS*n_threads);
    TimerReport(&tot_time,NULL);


    ErrorStreamClose();
    printf("Band calculate finished\n");
    fflush(stdout);
    return 0;
}
