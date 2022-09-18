/* Header file for pardiso
Copies the format of umfsolver.h
 */

#ifndef _PARDISOSOLVER_H
#define _PARDISOSOLVER_H

#include "mkl.h"
#include "complex.h"

#define PARDISO_DEBUGGING_MSGLVL 0
void LapackEigenSolve(int NUM_BANDS, int N_RANK, double complex  *A);

#endif
