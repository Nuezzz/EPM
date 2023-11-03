
#ifndef _BAND
#define _BAND
#include "atom.h"

#define SIMNAME     "SIMULATION"
#define SIMWAVE     "SIMULATION_WAVE"
#define NTHREADS    "N_THREADS"
#define SIMTYPE     "SIMULATION_TYPE"
#define POTTYPE     "POTENTIAL_TYPE"
#define KPATH       "K_PATH"
#define ECUT        "E_CUT"

#define NVALANCE    2
#define NCONDUCTION 4
typedef struct band
{
    char   *simtype;//simulation type

    unsigned int   pottype;//potential type         1/SO or 0/LOC
    unsigned int   simwave;//do not save wave function if simwave = 0, save wave function if simwave = 1

    Lattice *L; //lattice information 
   
    unsigned int     n_band;//number of band
    unsigned int     n_kpoint;//number of k points

    double Emax; //cut off energy (in ev)
    double *k_list;//list point of k_point in one k_path

    FILE *Band_File;
    FILE *Wave_File;
} BAND;

typedef struct epm
{
    char    *simname;   //simulation name 
    char    *simtype;   //simulation type      BAND or KLIST
    char    *kpath_name;//file name for stored k_path

    unsigned int   pottype;//potential type         2/SO or 1/LOC
    unsigned int   simwave;//do not save wave function if simwave = 0, save wave function if simwave = 1

    Lattice *L; //lattice information
    /*************** IF SIMTYPE is BAND ****************/
    unsigned int     n_kpath;//number of kpath
    unsigned int     n_k_p_path;//number of k points per k path
    double  *k_symetry;//array of high symetry k point
    BAND *band;//array of subband


    /*************** IF SIMTYPE is KLIST ****************/
    unsigned int     n_kpoint;//number of k points
    double  *k_list;//array of high symetry k point    

    unsigned int     n_band;//number of band
    unsigned int     n_threads;//number of threads

    double Emax; //cut off energy (in ev)
} EPM;

void EPMReadInput(EPM *s,const char *filename);
void EPMSimInit(EPM *s);
void EPMSimFree(BAND *s, unsigned int   simwave);
void Calc_EPM_Patch(double  *k_path, int n_k, double E_cut , int n_bands, unsigned int pottype, Lattice *L, FILE *Band_File);
void Calc_EPM_Patch_Wave(double  *k_path, int n_k, double E_cut , unsigned int pottype, Lattice *L, FILE *Band_File, FILE *Wave_File);

#endif