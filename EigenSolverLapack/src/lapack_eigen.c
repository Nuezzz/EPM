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

static inline void LapackComplexConvert(double complex *A, lapack_complex_double *a, int rank);

void LapackEigenSolve(int NUM_BANDS, int N_RANK, double complex  *A)
{


	lapack_int n  		= N_RANK;
	lapack_complex_double a[N_RANK*N_RANK];
	lapack_int lda		= N_RANK;
	double	   vl,vu;

	lapack_int il		= 1;//lower limit of eigenvalue index
	lapack_int iu		= NUM_BANDS;//upper limit of eigenvalue index

	double 	   abstol   = -1.0;

	lapack_int m;
	lapack_int ldz 	  	= N_RANK;
	double	   w[ldz];
	lapack_complex_double z[ldz*NUM_BANDS];

	lapack_int ifail[ldz];
	lapack_int info;

	LapackComplexConvert(A, a ,N_RANK);

	info = LAPACKE_zheevx( LAPACK_ROW_MAJOR, 'N', 'I', 'L', n, a, lda, vl, vu,  il, iu, abstol, &m, w, z, ldz, ifail);
	//info = LAPACKE_zheevx(  LAPACK_ROW_MAJOR, "V",  "I", "L", n, a, lda, vl, vu, il, iu, abstol, &m, w, z, ldz, ifail );


	if (info != 0) {
		printf("\nERROR during solution: %lld\n", info);
		exit(3);
	}

	2for(int i=0; i<m; i++){
		printf("the %d eigenvalue is %lf\n",i,w[i]);
	}

}

/**
 * @brief convert the commplex type into lapack complex type
 * 
 * @param A complex pointer to rank*rank size array
 * @param rank the size of input and outpur array
 * @return lapack_complex_double type of aarray
 */
static inline void LapackComplexConvert(double complex *A, lapack_complex_double *a, int rank)
{
	int i,j;
	
	for (i=0; i<rank; i++)
	{
		for (j=0; j<rank; j++)
		{
			a[i*rank+j].real= creal(A[i*rank+j]);
			a[i*rank+j].imag= cimag(A[i*rank+j]);
		}
	}
	return a;
}
