#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "atom.h"
#include "error.h"
#include "memory.h"
#include "util.h"
#include "ham.h"

//
//
int main(int argc, char **argv)
{
    Lattice *L;
    int     N;
    char   *filename;
    StringClone(filename, "GaAs_64.csv");
    int    Kmax =20;
    double E_cut = 50.0;
    double k[3];
    k[0]=0; k[1]=0; k[2]=0;

	printf("Empirical Sudopotential Calculation\n ");
    fflush(stdout);
    ErrorStreamOpen("error.log");

    if(argc == 2)
    {
        L = LatticeInitial(argv[1],2*Kmax+1);
        N = BuildG(L, E_cut, Kmax, k);
        HLocal();
    }

    else
    {
        printf("No tittle was selected\n");
    }
    PrintGvec(L, filename, N);
    ErrorStreamClose();
    return 0;
}
