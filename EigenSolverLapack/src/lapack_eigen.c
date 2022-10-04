/*
	Direct linear solver using LU decomposition in pardiso.
	Sparse matrix operations are handeled using intel mkl routines
	Pardiso performs the sparse LU factorization and solve

	Written Nov 2021
 */
#include <stdio.h>  
#include <stdlib.h>
#include <string.h>
#include "error.h"

#include "eigenlapack.h"


void print_matrix( char* desc, int m, int n, double complex *a, int lda );
void print_rmatrix( char* desc, int m, int n, double *a, int lda );

int LapackEigenSolve(int NUM_BANDS, int N_RANK, double complex  *A, double *w, double complex *z)
{


	lapack_int n  		= N_RANK;
	lapack_int lda		= N_RANK;
	double	   vl,vu;

	lapack_int il		= 1;//lower limit of eigenvalue index
	lapack_int iu		= NUM_BANDS;//upper limit of eigenvalue index

	double 	   abstol   = -1.0;

	lapack_int m;
	lapack_int ldz 	  	= N_RANK;
	


	lapack_int ifail[ldz];
	lapack_int info;

	//print_matrix( "Original matrix is (stored columnwise)", n, n, A, ldz );
	//info = LAPACKE_zheevx( LAPACK_ROW_MAJOR, 'N', 'I', 'L', n, a, lda, vl, vu,  il, iu, abstol, &m, w, z, ldz, ifail);
	info = LAPACKE_zheevx(  LAPACK_ROW_MAJOR, 'V',  'I', 'L', n, A, lda, vl, vu, il, iu, abstol, &m, w, z, ldz, ifail );


	if (info != 0) {
		printf("\nERROR during solution: %lld\n", info);
		exit(3);
	}

	return m;

}

void print_matrix( char* desc, int m, int n, double complex *a, int lda ) {
        int i, j;
        printf( "\n %s\n", desc );
        for( i = 0; i < m; i++ ) {
                for( j = 0; j < n; j++ )
                        printf( " (%6.2f,%6.2f)", creal(a[i*lda+j]), cimag(a[i*lda+j]) );
                printf( "\n" );
        }
}

/* Auxiliary routine: printing a real matrix */
void print_rmatrix( char* desc, int m, int n, double* a, int lda ) {
        int i, j;
        printf( "\n %s\n", desc );
        for( i = 0; i < m; i++ ) {
                for( j = 0; j < n; j++ ) printf( " %6.2f", a[i*lda+j] );
                printf( "\n" );
        }
}

