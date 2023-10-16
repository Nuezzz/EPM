#include <stdio.h>  
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ham.h"
#include "util.h"
#include "memory.h"
#include "constants.h"


/// @brief Constructed the hamitonian matrix for local potential
/// @param s pointer to the lattice
/// @param d 
/// @param Kvec 
/// @return 
double complex *H_tot_loc ( Lattice *s, Eigen *d, double *Kvec)
{
    int i;
    int NG = d->NG;
    double GpK[3];
    double **G_vec = d->G_vec;
    double complex *H= HLocal(s, d);

    // Declare variables for eigen solver

    // double	   w[NG];
    // double complex Z[NG*NG];

    for(i=0; i< NG; i++)
    {
        GpK[0]= G_vec[i][0]+Kvec[0];
        GpK[1]= G_vec[i][1]+Kvec[1];
        GpK[2]= G_vec[i][2]+Kvec[2];
        H[i*NG+i]+= H2M0Q0*Dot(GpK,GpK)*1.0e20;
    }

    return H;
    // print_matrix( "Selected eigenvectors (stored columnwise)", NG, NG, H, NG);
    // m = LapackEigenSolve(10, NG, H, w, Z);
    // print_rmatrix( "Selected eigenvalues", 1, m, w, 1 );

}

double complex *H_tot_so ( Lattice *s, Eigen *d, double *Kvec)
{
    int i;
    int NG = d->NG;
    double Ek_tmp;
    double GpK[3];
    double **G_vec = d->G_vec;
    double complex *H= HSO(s, d, Kvec);

    // Declare variables for eigen solver

    // double	   w[NG];
    // double complex Z[NG*NG];

    for(i=0; i< NG; i++)
    {
        GpK[0]= G_vec[i][0]+Kvec[0];
        GpK[1]= G_vec[i][1]+Kvec[1];
        GpK[2]= G_vec[i][2]+Kvec[2];
        Ek_tmp = H2M0Q0*Dot(GpK,GpK)*1.0e20;
        H[(i*4)*NG+2*i]+= Ek_tmp;
        H[(i*4+2)*NG+2*i+1]+= Ek_tmp;
    }
    //print_matrix( "Selected eigenvectors (stored columnwise)", 20, 10, H, 2*NG);
    return H;
    
}


/**
 * @brief Calculate the eigen energy and eigen vector of vertain k point
 * 
 * @param d pointer to the structure d that contains G_vector mesh and storage for eigen value, eigen vector
 * @param H hamitonial matrix
 * @param bands 
 * @param N_RANK
 * @return number of bands calulated
*/
int CalcBand(Eigen *d, double complex *H, int bands, int N_RANK)
{
    int     m=0;
    double	   *w= d->E;
    double complex *Z= d->Phi;

    //print_matrix( "Selected eigenvectors (stored columnwise)", 20, 10, H, N_RANK);
    m = LapackEigenSolve(bands, N_RANK, H, w, Z);

    free(H);
    return m;

}

/**
 * @brief print out the band structure of certain K point
 * 
 * @param fp file name for output
 * @param E  pointer to the eigen value
 * @param k  pointer to the k vector
 * @param N  number of bands been printed out
 */
void PrintEigen(FILE *fp, double *E, double *k, int N)
{

    int i;
    double mode;

    mode=sqrt(Dot(k,k));
    fprintf(fp, "%.15e %.15e %.15e %.15e", k[0],k[1],k[2],mode);
    for(i=0;i<N;i++)
    {
        fprintf(fp," %.15e",E[i]);
    }
    fprintf(fp,"\n");
}

FILE *OpenBandFile( char* simname)
{
    FILE *fp;
	char filename[128];
	sprintf(filename,"%s/Band_Stracture.dat",simname);
    fp = SafeFOpen(filename, "w");
    return fp;
}