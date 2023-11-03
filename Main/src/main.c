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
#include "epm.h"
#include <complex.h>





// /**
//  * @brief Found the G_vec sets corresponding to specific K_vec within the E_cut ,
//  * build the Hamitonian, solve it's eigen value and print out the Band_FIle
//  * @param L  
//  * @param k 
//  * @param E_cut 
//  * @param m 
//  */
// static inline Eigen* KPatch_so(Lattice *L, double *k, double E_cut,int *m)
// {
//     Eigen* Estack;
//     double complex* Hstack;
//     int n_band = L->a_set->n_atoms*4+L->a_set->n_atoms*5;

//     //TimerStart(&L->Form_time);

    
//     Estack = GVecInit;
//     Hstack = H_tot_so(L,Estack,k);

//     //print_matrix( "Selected eigenvectors (stored columnwise)", 20, 10, Hstack, Estack->NG*2);
//     *m =CalcBand( Estack, Hstack, n_band, Estack->NG*2);

//     return Estack;
// }


//create G_vec array for each K_vec in the first patch, each patch should 
//contain n_treads of K_vec
//write a main function with three arguments, simname, simtype, STEPS
//simtype can be SO or LOC
//STEPS is the number of K_vec in the first patch
//the program will create a folder named simname_band_structure
//and write the band structure in i

int main(int argc, char **argv)
{
    int     i;
    char   *inputfile = argv[1];
    

    Timer   tot_time;
    EPM *s = SafeCalloc(1,sizeof(EPM));
    BAND   *b;
    ErrorStreamOpen("error.log");
    // if the input file not exist
    if(!FileExists(inputfile))
    {
        free(inputfile);
        printf("No input file found.\n ");
        fflush(stdout);
        return 1;
    }

    TimerStart(&tot_time);
    EPMReadInput(s,inputfile);
    EPMSimInit(s);

    switch (s->simwave)
    {
    case 0://don't save wave function
        
        for(i=0;i<s->n_kpath;i++)
        {   
            b=s->band+i;
            Calc_EPM_Patch(b->k_list,b->n_kpoint,s->Emax,s->n_band,s->pottype,s->L,b->Band_File);
            EPMSimFree(b,0);
        }
        break;

    case 1://save wave function
        for(i=0;i<s->n_kpath;i++)
        {   
            b=s->band+i;
            Calc_EPM_Patch_Wave (b->k_list,b->n_kpoint,s->Emax,s->pottype,s->L,b->Band_File,b->Wave_File);
            EPMSimFree(b,1);
        }
        break;
    
    default:
        break;
    }
    TimerStop(&tot_time);
    printf("Total time spend for %d atoms, and %d k points:\n",s->L->a_set->n_atoms,s->n_kpath*s->n_k_p_path);
    TimerReport(&tot_time,NULL);

    ErrorStreamClose();
    printf("Band calculate finished\n");
    fflush(stdout);
    return 0;
}
