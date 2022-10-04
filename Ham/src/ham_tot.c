#include <stdio.h>  
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ham.h"
#include "util.h"
#include "memory.h"
#include "constants.h"



double complex *HTot ( Lattice *s, int NG, double *Kvec)
{
    int i,j;
    double GpK[3];
    double **G_vec = s->G_vec;
    double complex *H =  SafeCalloc(NG*NG, sizeof( double complex ));
    double complex *V_loc= HLocal(s, NG);

    // Declare variables for eigen solver

    // double	   w[NG];
    // double complex Z[NG*NG];


    for(i=0; i< NG; i++)
    {
        for(j=0; j< NG; j++)
        {
            H[i*NG+j] = V_loc[i*NG+j];
        }
    }
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

/**
 * @brief Calculate the eigen energy and eigen vector of vertain k point
 * 
 * @param s the lattice pointer
 * @param NG Number of valid G vectors within cutoff energy
 * @param Kvec k point in the reciprocal space
 * @param bands Number of selected bands needs to be calculated
 * @return number of bands calulated
 */
int CalcBand( Lattice *s, double *Kvec, int NG, int bands)
{
    int     m;
    double	   *w= s->E;
    double complex *Z= s->Phi;
    double complex *H = HTot(s,NG,Kvec);
    //print_matrix( "Selected eigenvectors (stored columnwise)", NG, NG, H, NG);
    m = LapackEigenSolve(bands, NG, H, w, Z);

    return m;

}

/**
 * @brief print out the band structure of certain K point
 * 
 * @param E 
 * @param simname 
 * @param k 
 * @param N 
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