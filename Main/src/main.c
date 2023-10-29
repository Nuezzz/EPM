#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#ifdef _OPENMP
#include <omp.h>
#endif
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
 * build the Hamitonian, solve it's eigen value and print out the Band_FIle
 * @param L  
 * @param k 
 * @param E_cut 
 * @param m 
 */
static inline Eigen* KPatch_loc( Lattice *L, double *k, double E_cut,  int* m)
{    
    Eigen* Estack;
    double complex* Hstack;
    int n_band = L->a_set->n_atoms*5;

    //TimerStart(&L->Form_time);

    Estack = GVecInit( L,  k, E_cut);
    Hstack = H_tot_loc(L,Estack,k);

    //print_matrix( "Selected eigenvectors (stored columnwise)", 20, 10, Hstack, Estack->NG);
    *m=CalcBand( Estack, Hstack, n_band, Estack->NG );
    return Estack;
  }

/**
 * @brief Found the G_vec sets corresponding to specific K_vec within the E_cut ,
 * build the Hamitonian including spin-obital, solve it's eigen value and print out the band 
 * @param L  
 * @param k 
 * @param E_cut 
 * @param m 
 */
static inline Eigen* KPatch_so(Lattice *L, double *k, double E_cut,int *m)
{
    Eigen* Estack;
    double complex* Hstack;
    int n_band = L->a_set->n_atoms*4+L->a_set->n_atoms*5;

    //TimerStart(&L->Form_time);

    
    Estack = GVecInit( L,  k, E_cut);
    Hstack = H_tot_so(L,Estack,k);

    //print_matrix( "Selected eigenvectors (stored columnwise)", 20, 10, Hstack, Estack->NG*2);
    *m =CalcBand( Estack, Hstack, n_band, Estack->NG*2);

    return Estack;
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
    FILE *Band_File, *NG_File;
    Lattice *L;

    Timer   tot_time;

    int     i;
    char   *foldername;
    char   *simname = argv[1];
    char   *simtype = argv[2];
    int    STEPS=atoi(argv[3]);
    int    n_threads = atoi(argv[4]);
    
    StringClone(foldername,simname);
    StrCat(&foldername, "_band_structure");

    // double complex *H_tot_loc;
    
    double  E_cut = 100.0;
    double  *k;
    double  *k_path;

    //***********test line*************************
    int N_G[STEPS];

    //**********************************************
	printf("Empirical Sudopotential Calculation\n ");
    fflush(stdout);
    ErrorStreamOpen("error.log");


    TimerStart(&tot_time);
    if(argc == 5)
    {
        CreateFolder(foldername);
        Band_File = OpenBandFile(foldername);
        NG_File   = OpenNGFile(foldername);
        OMPSetThreadNum(n_threads);

        L = LatticeInitial(simname);

        L->Form_time.tot_wtime=0;
        L->Form_time.tot_ctime=0;
        L->Solve_time.tot_wtime=0;
        L->Solve_time.tot_ctime=0;

        k_path=KBuild(L, STEPS);
        //PPtest(L,foldername);
        if(!strcmp(simtype,"SO"))
        {
            L->pot_type=1;
            #ifdef _OPENMP
            #pragma omp parallel for schedule(static) private(i,k)
            #endif
            for(i=0;i<STEPS;i++)
            {
                Eigen* Estack;
                int m=0;
                k= k_path+3*i;
            
                Estack = KPatch_so( L,k, E_cut,&m);

                #ifdef _OPENMP
                #pragma omp critical
                #endif
                PrintEigen(Band_File, Estack->E, k,m);

                N_G[i]=Estack->NG;                     
                BandFinish(Estack);

                //printf("The %dth K vector finished, %f %% \n ",i+1, (float)(i+1)/(float)(STEPS)*100);
                //fflush(stdout);
            }
            PrintNG(NG_File, N_G, k_path, STEPS);
        }

        else if(!strcmp(simtype,"LOC"))
        {
            L->pot_type=0;
            #ifdef _OPENMP
            #pragma omp parallel for schedule(static) private(i,k)
            #endif
            for(i=0;i<STEPS;i++)
            {
                k= k_path+3*i;

                int m=0;
                Eigen* Estack;
                Estack = KPatch_loc( L,k, E_cut,&m);
            
                #ifdef _OPENMP
                #pragma omp critical
                #endif
                PrintEigen(Band_File, Estack->E, k,m);

                BandFinish(Estack);
                //printf("The %dth K vector finished, %f %% \n ",i+1, (float)(i+1)/(float)(STEPS)*100);
                //fflush(stdout);
            }
        }

        fclose(Band_File);
        fclose(NG_File);
    }

    else
    {
        printf("No tittle was selected\n");
        return 1;
    }

    TimerStop(&tot_time);
    printf("Total time spend for %d atoms, and %d k points:\n",L->a_set->n_atoms,STEPS);
    TimerReport(&tot_time,NULL);

    ErrorStreamClose();
    printf("Band calculate finished\n");
    fflush(stdout);
    return 0;
}
