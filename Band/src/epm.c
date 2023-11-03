#include <stdio.h>  
#include <complex.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#endif
#include "openmp.h"
#include "atom.h"
#include "error.h"
#include "ham.h"
#include "epm.h"
#include "memory.h"
#include "util.h"
#include "reader.h"

/**
 * @brief Build the K vector mesh in reciprocal space
 *  
 * @param L pointer to lattice
 * @param N number of K vector
 * @param k_start start point of K vector
 * @param k_end end point of K vector
 * 
 * @return k_vector pointer to the K vector mesh
 **/
static inline double *K_Build(Lattice *L, int N, double *k_start, double *k_end)
{
    double unit_a   = 2.00*PI/L->A0;
    double k_0[3];
    double k_1[3];
    double dk[3];
    double *k_vector =SafeCalloc(3*N,sizeof(double));

    k_0[0]=unit_a*k_start[0];
    k_0[1]=unit_a*k_start[1];
    k_0[2]=unit_a*k_start[2];

    k_1[0]=unit_a*k_end[0];
    k_1[1]=unit_a*k_end[1];
    k_1[2]=unit_a*k_end[2]; 

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

/**
 * @brief Found the G_vec sets corresponding to specific K_vec within the E_cut ,
 * build the Hamitonian, solve it's eigen value and print out the Band_FIle
 * @param L  
 * @param k 
 * @param E_cut 
 * @param m 
 **/
static inline Eigen* Calc_K_loc( Lattice *L, double *k, double E_cut,  int n_band)
{    
    Eigen* eigen;
    double complex* Ham;
    int m=0;
    unsigned int ptype=0;//potential type LOC

    //TimerStart(&L->Form_time);
    eigen = GVecInit( L,  k, E_cut, ptype);
    Ham = H_tot_loc(L,eigen,k);

    //print_matrix( "Selected eigenvectors (stored columnwise)", 20, 10, Hstack, Estack->NG);
    m=CalcBand( eigen, Ham, n_band, eigen->NG );
    if(m==n_band)
    {
        return eigen;
    }
    else
    {
    char errstr[256];
    printf("Number of band %d calculated is not equal to the number of band %d requested\n", m,n_band);
    sprintf(errstr,"Number of band %d calculated is not equal to the number of band %d requested\n", m,n_band);
    ERROR_THROW(LAPACK_ERROR, errstr);
    ERROR_CATCH;
    }

}

/**
 * @brief Found the G_vec sets corresponding to specific K_vec within the E_cut ,
 * build the Hamitonian, solve it's eigen value and print out the Band_FIle
 * @param L  
 * @param k 
 * @param E_cut 
 * @param m 
 **/
static inline Eigen* Calc_K_so( Lattice *L, double *k, double E_cut,  int n_band)
{    
    Eigen* eigen;
    double complex* Ham;
    int m=0;
    unsigned int ptype=1;//potential type SO

    //TimerStart(&L->Form_time);
    eigen = GVecInit( L,  k, E_cut, ptype);
    Ham = H_tot_so(L,eigen,k);

    //print_matrix( "Selected eigenvectors (stored columnwise)", 20, 10, Hstack, Estack->NG);

    m=CalcBand( eigen, Ham, n_band, eigen->NG *2);
    if(m==n_band)
    {
        return eigen;
    }
    else
    {
    char errstr[256];
    printf("Number of band %d calculated is not equal to the number of band %d requested\n", m,n_band);
    sprintf(errstr,"Number of band %d calculated is not equal to the number of band %d requested\n", m,n_band);
    ERROR_THROW(LAPACK_ERROR, errstr);
    ERROR_CATCH;
    }

}

/*******
 * Initialize the EPM simulation
 * Create the output directory and output file
 * Read the Lattice information
 * ******/
void EPMSimInit(EPM *s)
{
    //
    BAND *subband;
    unsigned int    n_path;
    char   *foldername;
    StringClone(foldername,s->simname);
    StrCat(&foldername, "_result");
    printf("Start Empirical Psuedopotential Calculation\n");


    CreateFolder(foldername);

    s->L = LatticeInitial(s->simname);
    s->n_band= s->pottype*s->L->a_set->n_atoms*(NVALANCE+ NCONDUCTION);
    OMPSetThreadNum(s->n_threads);

    //for BAND type of simulation
    if(!strcmp(s->simtype,"BAND"))
    {
        double* k_start = 0;
        double* k_end   = 0;

        n_path  =   s->n_kpath;
        s->band =   SafeCalloc(n_path, sizeof(BAND));

        for(int i=0;i<n_path;i++)
        {
            k_start = s->k_symetry+3*i;
            k_end   = s->k_symetry+3*i+3;
            subband           = s->band+i;
            subband->n_kpoint = s->n_k_p_path;
            subband->n_band   = s->n_band;
            subband->k_list   = K_Build(s->L, subband->n_kpoint,k_start,k_end);
            subband->Emax     = s->Emax;
            subband->L        = s->L;
            subband->pottype  = s->pottype;
            subband->simwave  = s->simwave;
            subband->Band_File= OpenBandFile(foldername, i);
            if(s->simwave)
                subband->Wave_File= OpenWaveFile(foldername, i);
        } 
        
    }
    else if(!strcmp(s->simtype,"KLIST"))
    {
        //#code
        subband =   SafeCalloc(n_path, sizeof(BAND));

        subband->n_kpoint = s->n_kpoint;
        subband->n_band   = s->n_band;
        subband->k_list   = s->k_list;
        subband->Emax     = s->Emax;
        subband->L        = s->L;
        subband->pottype  = s->pottype;
        subband->simwave  = s->simwave;
        subband->Band_File= OpenBandFile(foldername, 0);
        if(s->simwave)
            subband->Wave_File= OpenWaveFile(foldername, 0);

        s->band = subband;

    }
}

/**
 * @brief Free the memory of the EPM simulation
 * 
 * @param s 
 * @param simwave indicate whether the wave function is saved
 
*/
void EPMSimFree(BAND *s, unsigned int   simwave)
{
fclose(s->Band_File);
if(simwave) fclose(s->Wave_File);
free(s->k_list);

printf("Band calculate finished Successfully\n");
fflush(stdout);
}


void Calc_EPM_Patch(double  *k_path, int n_k, double E_cut , int n_bands, unsigned int pottype, Lattice *L, FILE *Band_File)
{
    int i;
    double  *k;

    switch (pottype)
    {
    case 1://LOC
        #ifdef _OPENMP
        #pragma omp parallel for schedule(static) private(i,k)
        #endif
        for(i=0;i<n_k;i++)
        {
            Eigen *Eg;
            k= k_path+3*i;
            
            Eg = Calc_K_loc( L, k, E_cut,n_bands);

            #ifdef _OPENMP
            #pragma omp critical
            #endif
            PrintEigen(Band_File, Eg->E, k,n_bands);                     
            BandFinish(Eg);
            }
        
        break;
    case 2://SO
        #ifdef _OPENMP
        #pragma omp parallel for schedule(static) private(i,k)
        #endif

        for(i=0;i<n_k;i++)
        {
            Eigen *Eg;
            k= k_path+3*i;
            
            Eg = Calc_K_so( L, k, E_cut,n_bands);

            #ifdef _OPENMP
            #pragma omp critical
            #endif
            PrintEigen(Band_File, Eg->E, k,n_bands);                     
            BandFinish(Eg);
        }
        break;
        
    default:
        printf("Wrong simulation type\n");
        fflush(stdout);
        exit(0);
        break;
    }
    
}

void Calc_EPM_Patch_Wave(double  *k_path, int n_k, double E_cut , unsigned int pottype, Lattice *L, FILE *Band_File, FILE *Wave_File)
{
    //code
}