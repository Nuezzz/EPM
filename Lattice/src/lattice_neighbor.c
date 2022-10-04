#include <stdio.h>  
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include "memory.h"
#include "util.h"
#include "reader.h"
#include "atom.h"

static inline void RecordNeighbor (Lattice *l, unsigned int **neigh)
{
    char    *bond;
    int     i,j,k;
    int     n_atom = l->a_set->n_atoms;
    Atom    **atoms= l->a_set->atom_array;
    for(i=0;i<n_atom; i++)
    {
        for(k=0, j=0; j<l->n_spe; j++)
        {
            if(neigh[i][j]!= 0)   k++;
            
        }
        atoms[i]->neighbor_spe  =  SafeCalloc(k, sizeof(unsigned int));
        atoms[i]->n_neighbor    =  SafeCalloc(k, sizeof(unsigned int));
        atoms[i]->n_spe         = k;
        for(k=0, j=0; j<l->n_spe; j++)
        {
            if(neigh[i][j]!= 0)
            {
                atoms[i]->neighbor_spe[k] = j;
                atoms[i]->n_neighbor[k]   = neigh[i][j];
                k++;
            }
        }
        
    }
}


void FindNeighbor (double A_dis, Lattice *l)
{
    int     i,j,k,IREP,NI,NII;
    int     n_atom = l->a_set->n_atoms;
    double  DX,DY,DZ,D,DIS;
    double  a_vec[3][3];
    double  tmp[3];
    unsigned int **neigh;
    Atom **atoms= l->a_set->atom_array;

    neigh = SafeCalloc(n_atom, sizeof(unsigned int*));
    for(i=0; i<n_atom; i++)
    {
        neigh[i] = SafeCalloc(l->n_spe, sizeof(unsigned int));
    }
    a_vec[0][0]= l->a[0][0]; a_vec[0][1]= l->a[0][1]; a_vec[0][2]= l->a[0][2];
    a_vec[1][0]= l->a[1][0]; a_vec[1][1]= l->a[1][1]; a_vec[1][2]= l->a[1][2];
    a_vec[2][0]= l->a[2][0]; a_vec[2][1]= l->a[2][1]; a_vec[2][2]= l->a[2][2];

    DIS= pow(A_dis, 2.0);


    // Search all the 8 nearest cells for system smaller than 1000 
    if(n_atom<1000)
        IREP=1;
    else
        IREP=0;

    for ( NI = 0; NI < n_atom; NI++)
    {
        for ( NII = NI; NII < n_atom; NII++)
        {
            for (i=-IREP;i<IREP+1;i++)
            {
                for (j=-IREP;j<IREP+1;j++)
                {
                    for (k=-IREP;k<IREP+1;k++)
                    {
						DX = remainder(atoms[NI]->tau[0]-atoms[NII]->tau[0]+10.5,1.0) -0.5+i;
						DY = remainder(atoms[NI]->tau[1]-atoms[NII]->tau[1]+10.5,1.0) -0.5+j;
						DZ = remainder(atoms[NI]->tau[2]-atoms[NII]->tau[2]+10.5,1.0) -0.5+k;

                        tmp[0] = DX*a_vec[0][0]+DY*a_vec[1][0]+DZ*a_vec[2][0];
                        tmp[1] = DX*a_vec[0][1]+DY*a_vec[1][1]+DZ*a_vec[2][1];
                        tmp[2] = DX*a_vec[0][2]+DY*a_vec[1][2]+DZ*a_vec[2][2];
                        D =  tmp[0]*tmp[0]+tmp[1]*tmp[1]+tmp[2]*tmp[2];      

                        if (NII!=NI && D < DIS )
                        {
                            neigh[NI][atoms[NII]->spe]++;
                            neigh[NII][atoms[NI]->spe]++;
                        }
            
                    }
                }
            }
        }
        
    }
    RecordNeighbor (l, neigh);

}
