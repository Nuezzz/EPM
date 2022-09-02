typedef struct atom
{
    double tau[3]; //direct atomic position
    unsigned int spe; // specie of the atom
    unsigned int *n_neighbor; //array of number of neighbors from various type of atoms

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

typedef struct param
{
	double q[4];
	double p_s;
	double p_m;
} Param;

typedef struct paramstack
{
    Param *par_stotage; 
    Param **par_array; // array of pointers to the bond atom type
    /* data */
}ParamStack;


typedef struct lattice
{
    char   title[60];
    /*Info of atom species*/
	unsigned int n_spe; // number of atom spycies
    unsigned int *natom_spe;
    char    **mat_name;
    char    *mat_storage;

	double a[3][3]; //lattice vector

    unsigned int n_kstep; // number of strps per k-path
    double Emax; //cut off energy (in Rydberg)
	AtomStack *a_set;// Atom sets
    ParamStack  *pset_storage;
    ParamStack  **pset_array; // array of pointers to the base atom type
} Lattice;