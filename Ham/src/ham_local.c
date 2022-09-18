#include <stdio.h>  
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include "ham.h"
#include "constants.h"
#include "atom.h"


void HLocal()
{
    double complex A[4*4] =
    {
        6.51 +0.00*  I,  0.00 +0.00*  I,  0.00+ 0.00*  I,  0.00+ 0.00* I,
       -5.92 -9.53*  I, -1.73 +0.00*  I,  0.00+ 0.00*  I,  0.00+ 0.00* I,
       -2.46 -2.91*  I,  6.50 -2.09*  I,  6.90+ 0.00*  I,  0.00+ 0.00* I,
        8.84 -3.21*  I,  1.32 -8.81*  I, -0.59- 2.47*  I, -2.85+ 0.00* I
    };
    LapackEigenSolve(2, 4, A);

}



