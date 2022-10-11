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

#define  KMAX  20
#define  STEPS 10

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
    Eigen   **Estack;
    double complex **Hstack;
    int n_band = L->a_set->n_atoms*2+8;
    int i,thread;
    int m[n_threads];

    Estack=SafeCalloc(n_threads, sizeof(Eigen*));
    Hstack=SafeCalloc(n_threads, sizeof(double complex*));

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
    printf("Finish Building the G_vec mesh\n");
    fflush(stdout);

    for(i=0;i<n_threads;i++)
    {
        Hstack[i]=HTot(L,Estack[i],k[i]);
        free(Estack[i]->G_vec);
        printf("Finish Building the Hamitonian for %d th \n", i+1);
        fflush(stdout); 
    }
    
       
    

        #ifdef _OPENMP
        #pragma omp parallel private(thread)
        #endif
        {
            #ifdef _OPENMP
            thread = omp_get_thread_num();
            #else
            thread = 0;
            #endif
            m[thread]=CalcBand( Estack[thread], Hstack[thread], n_band);
        }
        
        //print_rmatrix( "Selected eigenvalues", 1, m, L->E, 1 );

        for(i=0;i<n_threads;i++)
        {
            PrintEigen(band, Estack[i]->E, k[i],m[i]);
            free(Estack[i]);
            free(Hstack[i]);
        }
        free(Estack);
        free(Hstack);

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


    int     i,j;
    char   *foldername;
    char   *simname = argv[1];
    int    n_threads = atoi(argv[2]);
    
    StringClone(foldername,simname);
    StrCat(&foldername, "_band_structure");
    

    // double complex *H_tot;
    
    double  E_cut = 80.0;
    double  *k[n_threads];
    double  *k_path;

       


	printf("Empirical Sudopotential Calculation\n ");
    fflush(stdout);
    ErrorStreamOpen("error.log");

    if(argc == 3)
    {
        CreateFolder(foldername);
        band = OpenBandFile(foldername);
        OMPSetThreadNum(n_threads);


        L = LatticeInitial(simname);
        k_path=KBuild(L, STEPS*n_threads);
        PPtest(L,foldername);

        for(i=0;i<STEPS;i++)
        {
            for(j=0;j<n_threads;j++)
            {
                k[j]= k_path+3*(i*STEPS+j);
            }
            KPatch(band, L,k, E_cut, n_threads );
            printf("The %d patch finished, %f %% \n ",i+1, (float)(i+1)/(float)(STEPS)*100);
            fflush(stdout);
        }
    
        fclose(band);
    }

    else
    {
        printf("No tittle was selected\n");
        return 1;
    }
    
    ErrorStreamClose();
    printf("Band calculate finished\n");
    fflush(stdout);
    return 0;
}
