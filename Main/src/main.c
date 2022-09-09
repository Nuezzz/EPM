#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "atom.h"
#include "error.h"
#include "memory.h"
#include "util.h"

//
int main(int argc, char **argv)
{
    Lattice *L;
	printf("Empirical Sudopotential Calculation\n ");
    fflush(stdout);
    ErrorStreamOpen("error.log");
    if(argc == 2)
    {
        L = LatticeInitial(argv[1]);
    }
    else
    {
        printf("No tittle was selected\n");
    }
    ErrorStreamClose();
    return 0;
}
