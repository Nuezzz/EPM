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
static inline void KPatch_loc(FILE *band, Lattice *L, double *k, double E_cut)
{    
    Eigen* Estack;
    double complex* Hstack;
    int n_band = L->a_set->n_atoms*5;
    int m=0;

    //TimerStart(&L->Form_time);

    Estack = GVecInit( L,  KMAX,k, E_cut);
    Hstack = H_tot_loc(L,Estack,k);

    //print_matrix( "Selected eigenvectors (stored columnwise)", 20, 10, Hstack, Estack->NG);
    m=CalcBand( Estack, Hstack, n_band, Estack->NG );
    PrintEigen( band, Estack->E, k,m );
    BandFinish( Estack );
  }

/**
 * @brief Found the G_vec sets corresponding to specific K_vec within the E_cut ,
 * build the Hamitonian including spin-obital, solve it's eigen value and print out the band
 * @param band 
 * @param L 
 * @param Estack 
 * @param Hstack 
 * @param k 
 * @param E_cut 
 * @param n_threads 
 */
static inline void KPatch_so(FILE *band, Lattice *L, double *k, double E_cut)
{
    
    Eigen* Estack;
    double complex* Hstack;
    int n_band = L->a_set->n_atoms*4+L->a_set->n_atoms*5;
    int m;

    //TimerStart(&L->Form_time);

    

    Estack = GVecInit( L,  KMAX,k, E_cut);
    Hstack = H_tot_so(L,Estack,k);

    //print_matrix( "Selected eigenvectors (stored columnwise)", 20, 10, Hstack, Estack->NG*2);
    m =CalcBand( Estack, Hstack, n_band, Estack->NG*2);

    PrintEigen(band, Estack->E, k,m);
    BandFinish(Estack);



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
    //print out the how many K vector is created
    printf("K vector created: %d\n",N);
    fflush(stdout);
    return k_vector;
    
}


//create G_vec array for each K_vec in the first patch, each patch should 
//contain n_treads of K_vec
//write a main function with three arguments, simname, simtype, STEPS
//simtype can be SO or LOC
//STEPS is the number of K_vec in the first patch
//the program will create a folder named simname_band_structure
//and write the band structure in i


int main(int argc, char **argv)
{
    FILE *band;
    Lattice *L;

    Timer   tot_time;

    int     i;
    char   *foldername;
    char   *simname = argv[1];
    char   *simtype = argv[2];
    int    STEPS=atoi(argv[3]);
    int    n_threads = 1;
    
    StringClone(foldername,simname);
    StrCat(&foldername, "_band_structure");
    

    // double complex *H_tot_loc;
    
    double  E_cut = 80.0;
    double  *k;
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
        if(!strcmp(simtype,"SO"))
        {
            L->pot_type=1;
            for(i=0;i<STEPS;i++)
            {
                k= k_path+3*i;
            
                KPatch_so(band, L,k, E_cut);
                printf("The %dth K vector finished, %f %% \n ",i+1, (float)(i+1)/(float)(STEPS)*100);
                fflush(stdout);
            }
        }

        else if(!strcmp(simtype,"LOC"))
        {
            L->pot_type=0;
            for(i=0;i<STEPS;i++)
            {
                k= k_path+3*i;
            
                KPatch_loc(band, L,k, E_cut);
                printf("The %dth K vector finished, %f %% \n ",i+1, (float)(i+1)/(float)(STEPS)*100);
                fflush(stdout);
            }
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
