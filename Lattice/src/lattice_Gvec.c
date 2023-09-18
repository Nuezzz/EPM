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

void PrintGvec(Eigen *s, char* simname, int N)
{
    FILE *fp;
    int i;
    double mode;
	char filename[128];
	sprintf(filename,"%s/G_vectors.csv",simname);
    fp = SafeFOpen(filename, "w");
    for(i=0;i<N;i++)
    {
        mode=s->G_vec[i][0]*s->G_vec[i][0]+s->G_vec[i][1]*s->G_vec[i][1]+s->G_vec[i][2]*s->G_vec[i][2];
        fprintf(fp, "%.15e %.15e %.15e %.15e\n", s->G_vec[i][0],s->G_vec[i][1],s->G_vec[i][2],mode);
    }
    fclose(fp);
}



/**
 * @brief Initialize and calloc the arrays for storing eigen value and eigen ectors
 * 
 * @param s pointer to lattice 
 * @param N number of valid G vectors for specific k point
 */
static inline void BandInit(Eigen *s, int N)
{
	s->G_vec 	=  SafeCalloc(N, sizeof(double*));
	s->G_stack 	=  SafeCalloc(3*N, sizeof(double));
	for(int i=0;i<N;i++)
	{
		s->G_vec[i]=s->G_stack+3*i;
	}
  	s->E 		=  SafeCalloc(N, sizeof(double));
	s->Phi 		=  SafeCalloc(N*N, sizeof(double complex));
}

void BandFinish(Eigen *s)
{
	free(s->G_vec);
	free(s->G_stack);
	free(s->E);
	free(s->Phi);
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
static inline int BuildG(Lattice *s,Eigen *d, double E_cut,  int Kmax)
{
	int i,j,k;

	int NG = 0; //count of number of valid G-vectors 
	int G_int[3]; 
	double G_tmp[3]; 

	double Gk_mod; // square of length of current G+k vector 
	double G_max; // square of max length of G vector
	double b[3][3];//reciprocal lattice vector bases
	double *G_bulk   =  SafeCalloc(3*Kmax*Kmax*Kmax, sizeof(double) );
	
	/*---------------------------------*/
	// FILE *fp;
    // int mode;
	// char filename[128];
	// sprintf(filename,"G_vectors.csv");
    // fp = SafeFOpen(filename, "w");
	 /******************************/

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

				Gk_mod= (G_tmp[0]+d->k_vec[0])*(G_tmp[0]+d->k_vec[0])+(G_tmp[1]+d->k_vec[1])*(G_tmp[1]+d->k_vec[1])+(G_tmp[2]+d->k_vec[2])*(G_tmp[2]+d->k_vec[2]);
				
				if(Gk_mod < G_max)
				{
					//mode=G_int[0]*G_int[0]+G_int[1]*G_int[1]+G_int[2]*G_int[2];
    				//fprintf(fp, "%d, %d, %d, %d, %.15e, %.15e, %.15e, %.15e\n", mode, G_int[0], G_int[1], G_int[2], G_tmp[0]/(2*PI/5.65359),G_tmp[1]/(2*PI/5.65359),G_tmp[2]/(2*PI/5.65359), Gk_mod/(2*PI/5.65359)/(2*PI/5.65359));
					
					G_bulk[NG*3+0]=G_tmp[0];
					G_bulk[NG*3+1]=G_tmp[1];
					G_bulk[NG*3+2]=G_tmp[2];
					NG += 1;
				}
			}
		}
	}

	BandInit(d,NG);

	for ( i = 0; i < NG; i++)
	{
		d->G_vec[i][0]=G_bulk[i*3+0];
		d->G_vec[i][1]=G_bulk[i*3+1];
		d->G_vec[i][2]=G_bulk[i*3+2];
	}
	
	//fclose(fp);
	free(G_bulk);
	return NG;
}


/**
 * @brief Build the reciprocal lattice mesh within the cut off, initialize 
 * the eigen energy and eigen fucntion storage
 * 
 * @param L lattice pointer that has already been fulfilled
 * @param K_max maximun number of k vector nodes in one direction, typically larger than  11
 * @param NG pointer to output the number of G_Vec within the cut off
 * @param k_vec pointer to pass the central k_vec
 * @param E_cut cut off energy
 * @return Eigen* 
 */
Eigen *GVecInit( Lattice *L, int K_max,double *k_vec, double E_cut)
{
	Eigen *s=SafeCalloc(1, sizeof(Eigen));
	int N = 2*K_max*((int)pow(L->a_set->n_atoms,1.0/3.0)+1)+1;
	s->k_vec[0]=k_vec[0]; s->k_vec[1]=k_vec[1]; s->k_vec[2]=k_vec[2];
	s->NG= BuildG(L,s,E_cut,N);

	if(s->NG ==  0)
	{
		printf(" Fail to select the G_mesh \n");
		exit(3);
	}
	return s;
}