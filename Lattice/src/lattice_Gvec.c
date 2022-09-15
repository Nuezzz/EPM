#include <stdio.h>  
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include "memory.h"
#include "constants.h"
#include "util.h"
#include "reader.h"
#include "atom.h"

void PrintGvec(Lattice *s, char* filename, int N)
{
    FILE *fp;
    int i;
    double mode;
    fp = SafeFOpen(filename, "w");
    for(i=0;i<N;i++)
    {
        mode=s->G_vec[i][0]*s->G_vec[i][0]+s->G_vec[i][1]*s->G_vec[i][1]+s->G_vec[i][2]*s->G_vec[i][2];
        fprintf(fp, "%.15e %.15e %.15e %.15e\n", s->G_vec[i][0],s->G_vec[i][1],s->G_vec[i][2],mode);
    }
    fclose(fp);

}
/**
 * @brief Build the G vector mesh in reciprcal space for 
 * specific K_vector, within the cutoff energy
 * 
 * @param s lattice pointer
 * @param E_cut cut off energy (eV)
 * @param Kmax  maximun of G vector index
 * @param k_vec pointer of k_vector use as the center of G vector
 */
int BuildG(Lattice *s, double E_cut,  int Kmax, double *k_vec)
{
	int i,j,k;

	int NG = 0; //count of number of valid G-vectors 
	int G_int[3]; 
	double G_tmp[3]; 

	double Gk_mod; // square of length of current G+k vector 
	double G_max; // square of max length of G vector
	double b[3][3];//reciprocal lattice vector bases

	G_max = E_cut * IH2M0Q0* 1.0e-20 ;//CONVERT A^-2 to m^-2
	for(i=0; i<3; i++)
	{
		b[i][0]=s->b[i][0];
		b[i][1]=s->b[i][1];
		b[i][2]=s->b[i][2];
	}
	for(i=0; i<2*Kmax+1; i++)
	{
		G_int[0]=i-Kmax;
		for(j=0; j<2*Kmax+1; j++)
		{
			G_int[1]=j-Kmax;
			for(k=0; k<2*Kmax+1; k++)
			{
				G_int[2]=k-Kmax;
				G_tmp[0] = G_int[0]*b[0][0]+G_int[1]*b[1][0]+G_int[2]*b[2][0];
				G_tmp[1] = G_int[0]*b[0][1]+G_int[1]*b[1][1]+G_int[2]*b[2][1];
				G_tmp[2] = G_int[0]*b[0][2]+G_int[1]*b[1][2]+G_int[2]*b[2][2];

				Gk_mod= (G_tmp[0]+k_vec[0])*(G_tmp[0]+k_vec[0])+(G_tmp[1]+k_vec[1])*(G_tmp[1]+k_vec[1])+(G_tmp[2]+k_vec[2])*(G_tmp[2]+k_vec[2]);
				
				if(Gk_mod < G_max)
				{
					s->G_vec[NG][0]=G_tmp[0];
					s->G_vec[NG][1]=G_tmp[1];
					s->G_vec[NG][2]=G_tmp[2];
					NG ++;
				}
			}
		}
	}
	return NG;
}