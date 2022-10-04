/* Header file for pardiso
Copies the format of umfsolver.h
 */

#ifndef _PARDISOSOLVER_H
#define _PARDISOSOLVER_H

#include "mkl.h"
#include <complex.h>

int LapackEigenSolve(int NUM_BANDS, int N_RANK, double complex  *A, double *w, double complex *z);
void print_matrix( char* desc, int m, int n, double complex *a, int lda );
void print_rmatrix( char* desc, int m, int n, double *a, int lda );
#endif
