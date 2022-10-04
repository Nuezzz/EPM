#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "atom.h"
#include "error.h"
#include "memory.h"
#include "util.h"
#include "ham.h"
#include <complex.h>


//
//
int main(int argc, char **argv)
{
    FILE *band;
    Lattice *L;
    int     m,N;
    char   *foldername;

    
    StringClone(foldername, argv[1]);
    StrCat(&foldername, "_band_structure");
    

    // double complex *H_tot;
    int    Kmax =20;
    double E_cut = 80.0;
    double k[3];
    k[0]=0; k[1]=0; k[2]=0;



	printf("Empirical Sudopotential Calculation\n ");
    fflush(stdout);
    ErrorStreamOpen("error.log");

    if(argc == 2)
    {
       
        CreateFolder(foldername);
        band = OpenBandFile(foldername);
        L = LatticeInitial(argv[1],2*Kmax+1);
        N = BuildG(L, E_cut, Kmax, k);
        //PrintGvec(L, foldername, N);
        //double complex *Z = SafeCalloc(N*N, sizeof( double complex ));
        PPtest(L,foldername);
        int n_band = L->a_set->n_atoms*2+8;
        m=CalcBand( L, k, N, n_band);
        //print_rmatrix( "Selected eigenvalues", 1, m, L->E, 1 );
        PrintEigen(band, L->E, k,m);



        fclose(band);
    }

    else
    {
        printf("No tittle was selected\n");
    }
    
    ErrorStreamClose();
    return 0;
}
