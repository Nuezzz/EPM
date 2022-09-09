#include <stdio.h>  
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include "memory.h"
#include "util.h"
#include "reader.h"
#include "atom.h"

void FindNeighbor (double A_dis, Lattice *l)
{
    int     i,j,k,IREP,NI,NII;
    int     n_atom = 0;
    double  DX,DY,DZ,D,DIS;
    double  a_vec[3][3];
    double  tmp[3];
    Atom **atoms= l->a_set->atom_array;

    a_vec[0][0]= l->a[0][0]; a_vec[0][1]= l->a[0][1]; a_vec[0][2]= l->a[0][2];
    a_vec[1][0]= l->a[1][0]; a_vec[1][1]= l->a[1][1]; a_vec[1][2]= l->a[1][2];
    a_vec[2][0]= l->a[2][0]; a_vec[2][1]= l->a[2][1]; a_vec[2][2]= l->a[2][2];

    DIS= pow(A_dis, 2.0);

    for(i=0; i<l->n_spe; i++)
    {
        n_atom+=l->natom_spe[i];
    }
    // Search all the 8 nearest cells for system smaller than 1000 
    if(n_atom<1000)
        IREP=1;
    else
        IREP=0;

    for ( NI = 0; NI < n_atom; NI++)
    {
        for ( NII = 0; NII < n_atom; NII++)
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
                            atoms[NI]->n_neighbor[atoms[NII]->spe]++;
                            atoms[NII]->n_neighbor[atoms[NI]->spe]++;
                        }
            
                    }
                }
            }
        }
        
    }

}