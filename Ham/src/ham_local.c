#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ham.h"
#include "memory.h"
#include "util.h"
#include "constants.h"

/**
 * @brief print out the band structure of certain K point
 *
 * @param E
 * @param simname
 * @param k
 * @param N
 */

static inline FILE *OpenPPFile(char *simname)
{
    FILE *fp;
    char filename[128];
    sprintf(filename, "%s/PP.dat", simname);
    fp = SafeFOpen(filename, "w");
    return fp;
}

/**
 * @brief Compute the atomic form factor Vf(q) for respect atomic species
 *
 * @param i the species of designated atom
 * @param j the species of neighbor atom
 * @param s the lattice pointer
 * @param q_norm modulus of q in (a.u.)^-1
 * @return form factor in (eV)
 */
static inline double FormFactor(unsigned int i, unsigned int j, Lattice *s, double q_norm)
{
    /// double Omega = 2.0e-30*s->vol/pow(BOHR,3.00)/s->a_set->n_atoms;	/* atomic volume, i.e. fcc cube/4	*/
    double Omega = 1.0e-30 * s->vol / pow(BOHR, 3.00);
    double Va; /* the atomic potential			*/
    char *bond;
    StringClone(bond, s->mat_name[i]);
    StrCat(&bond, s->mat_name[j]);
    if (!strcmp(bond, "GaAs"))
    {
        Va = (-1.24498 * exp(-1.52748 * pow((q_norm - 0), 2)) + 0.0366517 * exp(-0.959082 * pow((q_norm - 2.09782), 2)) + 0.0464357 * exp(-0.574047 * pow((q_norm - 2.01935), 2)) - 0.0133385 * exp(-11.2708 * pow((q_norm - 2.93581), 2))) * 131.4 / Omega;
        free(bond);
        return (Va * RYD); /* factor converts Rydberg --> eV */
    }

    if (!strcmp(bond, "AsGa"))
    {
        Va = (-1.0582 * exp(-0.959327 * pow((q_norm - 0), 2)) - 0.00217627 * exp(-6.53145 * pow((q_norm - 2.46808), 2)) - 0.0434312 * exp(-2.94679 * pow((q_norm - 0.851644), 2)) + 0.10569 * exp(-0.820922 * pow((q_norm - 1.22436), 2))) * 145.2 / Omega;
        free(bond);
        return (Va * RYD); /* factor converts Rydberg --> eV */
    }
    free(bond);
    return 0;
}

void PPtest(Lattice *s, char *simname)
{
    FILE *fp = OpenPPFile(simname);
    double Vf1, Vf0;
    // double max  = 6.0;
    // double step = max/100.0;
    double q[4];
    q[0] = sqrt(3);
    q[1] = sqrt(4);
    q[2] = sqrt(8);
    q[3] = sqrt(11);
    fprintf(fp, "%.9e \n ", 2.0e-30 * s->vol / pow(BOHR, 3.00) / s->a_set->n_atoms);
    for (int i = 0; i < 4; i++)
    {
        Vf1 = FormFactor(1, 0, s, q[i] * 2 * PI / s->A0 * 1.0e10 * BOHR);
        Vf0 = FormFactor(0, 1, s, q[i] * 2 * PI / s->A0 * 1.0e10 * BOHR);
        fprintf(fp, " %.9e, %.9e, %.9e \n", q[i], Vf0, Vf1);
    }
    fclose(fp);
}

// static inline double complex PotentialMix(Lattice *s,double *q)
// {
//     Atom **atoms= s->a_set->atom_array;
//     int n_atom=s->a_set->n_atoms;
//     int n_spe =s->n_spe;
//     int i,j;
//     double pot;
//     double complex V_loc;
//     double Vf[n_spe][n_spe];

//     double q_sqr = sqrt(Dot(q,q))*1.0e10*BOHR; // in unit of 1/Bohr
//     for( i=0; i<n_spe; i++)
//         for ( j = 0; j < n_spe; j++)
//         {
//             //Vf[i][j]= (i==j)? 0.0:FormFactor(i,j,s,q_sqr);
//             Vf[i][j]= 1;
//         }

//     for( i=0; i<n_atom; i++)
//     {
//         pot = 0;
//         for( j=0; j < atoms[i]->n_spe ; j++)
//         {
//             pot+=Vf[atoms[i]->spe][atoms[i]->neighbor_spe[j]]*atoms[i]->n_neighbor[j];
//         }
//         V_loc += cexp(-I*Dot(q,atoms[i]->tau))*pot/4.0;

//     }

//     return V_loc;
// }

/// @brief Use alloy linear combination to calculate the atomic potential at certain q point
/// @param s lattice vector for the system
/// @param q q_vector in unit of 1/Bohr within first Brillouin zone
/// @param V_loc the pointer to pass the output value of local potential
/// @param NG number of G vectors
double complex PotentialMix(Lattice *s, double q[3])
{

    double q_sqr = sqrt(Dot(q, q)) * 1.0e10 * BOHR; // in unit of 1/Bohr

    Atom **atoms = s->a_set->atom_array;
    Atom *atom;
    int n_atom = s->a_set->n_atoms;
    int n_spe = s->n_spe;
    // double pot[n_atom];
    double pot;
    double complex V_tmp = 0;
    double Vf[n_spe][n_spe];

    for (int i = 0; i < n_spe; i++)
        for (int j = 0; j < n_spe; j++)
        {
            Vf[i][j] = (i == j) ? 0.0 : FormFactor(i, j, s, q_sqr);
            // Vf[i][j]= 1;
        }

    for (int i = 0; i < n_atom; i++)
    {
        pot = 0;
        atom = atoms[i];
        for (int j = 0; j < atom->n_spe; j++)
        {
            // pot[i]+=Vf[atoms[i]->spe][atoms[i]->neighbor_spe[j]]*atoms[i]->n_neighbor[j];
            pot += Vf[atom->spe][atom->neighbor_spe[j]] * atom->n_neighbor[j];
        }
        // pot[i]/=4;
        pot /= 4;
        // V_loc[m*NG+n] += cexp(-I*Dot(q,atoms[i]->tau))*pot[i];
        V_tmp += cexp(-I * Dot(q, atom->tau)) * pot;
    }
    return V_tmp;
}

/// @brief
/// @param s
/// @param d
/// @return
double complex *HLocal(Lattice *s, Eigen *d)
{
    int NG = d->NG;
    double **G_vec = d->G_vec;
    double complex *V_loc = SafeCalloc(NG * NG, sizeof(double complex));

    int i, j;
    double q[3] = {0, 0, 0};

    for (i = 0; i < NG; i++)
    {
        for (j = 0; j <= i; j++)
        {
            q[0] = G_vec[i][0] - G_vec[j][0];
            q[1] = G_vec[i][1] - G_vec[j][1];
            q[2] = G_vec[i][2] - G_vec[j][2];
            V_loc[i*NG+j]=PotentialMix(s, q);
        }
    }
    return V_loc;
}


// {
//     int        n=4;
//     int        m;
//     double	   w[n];
//     double complex *Z = SafeCalloc(n*n, sizeof( double complex ));
//     double complex A[4*4] =
//     {
//         6.51 +0.00*  I,  0.00 +0.00*  I,  0.00+ 0.00*  I,  0.00+ 0.00* I,
//        -5.92 -9.53*  I, -1.73 +0.00*  I,  0.00+ 0.00*  I,  0.00+ 0.00* I,
//        -2.46 -2.91*  I,  6.50 -2.09*  I,  6.90+ 0.00*  I,  0.00+ 0.00* I,
//         8.84 -3.21*  I,  1.32 -8.81*  I, -0.59- 2.47*  I, -2.85+ 0.00* I
//     };
//     m = LapackEigenSolve(3, 4, A, w, Z);
//     printf("Number of eigen value calculated %d \n",m);
//     print_rmatrix( "Selected eigenvalues", 1, m, w, 1 );
//     print_matrix( "Selected eigenvectors (stored columnwise)", n, m, Z, n );

// }
