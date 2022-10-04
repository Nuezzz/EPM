#ifndef _LATTICE
#define _LATTICE

#include <complex.h>
typedef struct atom
{
    double tau[3]; //direct atomic position
    unsigned int spe; // specie of the atom
    unsigned int *n_neighbor; //array of number of neighbors from various type of atoms
    unsigned int *neighbor_spe;//array of neighbor atom species
    unsigned int n_spe; // number of neighbor atom species
}Atom;

typedef struct atomstack
{
    /*atoms array*/
    Atom *atom_storage;
    Atom **atom_array;
    unsigned int n_atoms;

    /*neighbor information array*/
    unsigned int *neigh_storage;
    // unsigned int **neigh_arary;
}AtomStack;

// typedef struct param
// {
// 	double p_a[4];      // Ryd
// 	double p_b[4];      // 1/a.u.
// 	double p_c[4];      // a.u.^2
//     double p_v;         // a.u.^3
// } Param;

// typedef struct paramstack
// {
//     Param *par_storage; 
//     Param **par_array; // array of pointers to the bond atom type
//     /* data */
// }ParamStack;


typedef struct lattice
{
    /*Info of atom species*/
	unsigned int n_spe; // number of atom spycies
    unsigned int *natom_spe;
    char    **mat_name;
    char    *mat_storage;

	double a[3][3]; //lattice vector (A)
    double b[3][3]; //reciprocal lattice vector (1/A)

    double A0; // lattice constant
    double vol;// the volume of unit cell  (A^3)

    unsigned int n_kstep; // number of strps per k-path
    double Emax; //cut off energy (in Rydberg)
    double **G_vec;//array of G vectors

    double complex *Phi;//eigen State
    double          *E;// eigen energy

	AtomStack   *a_set;// Atom sets
    // ParamStack  *pset_storage;
    // ParamStack  **pset_array; // array of pointers to the base atom type
} Lattice;

Lattice *LatticeInitial ( char *title, int K_max );
void    FindNeighbor    (double A_dis, Lattice *l);
void    PrintGvec       (Lattice *s, char* filename, int N);
int     BuildG          (Lattice *s, double E_cut,  int Kmax, double *k_vec);

#endif 