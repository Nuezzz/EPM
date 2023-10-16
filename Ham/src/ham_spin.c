#include <stdio.h>  
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ham.h"
#include "memory.h"
#include "util.h"
#include "constants.h"


static inline double *lambda(Lattice *s)
{
	int i;
	int n_spe = s->n_spe;
	
	double *lamb = SafeCalloc(n_spe, sizeof(double));
	for(i=0;i<n_spe;i++)
	{
    // TODO: Make this a switch-case over an enum?
    	if(!strcmp(s->mat_name[i],"Si"))  lamb[i]=(0.000106*RYD);
    	else if(!strcmp(s->mat_name[i],"Ge"))  lamb[i]=(0.00058*RYD);
    	else if(!strcmp(s->mat_name[i],"Ga"))	lamb[i]=(0.000402*RYD);
    	else if(!strcmp(s->mat_name[i],"As"))	lamb[i]=(0.000402*1.38*RYD);
    	else if(!strcmp(s->mat_name[i],"Cd"))	lamb[i]=(0.055*(0.0250-0.0090)*RYD);
    	else if(!strcmp(s->mat_name[i],"Te"))	lamb[i]=(0.055*(0.0250+0.0090)*RYD);
    	else if(!strcmp(s->mat_name[i],"CDTEcb"))	lamb[i]=(0.343*(0.002-0.0002)*RYD);
    	else if(!strcmp(s->mat_name[i],"TECDcb"))	lamb[i]=(0.343*(0.002+0.0002)*RYD);
    	else
    	{
			printf("Undefined Lambda for material %s",s->mat_name[i]);
    	}
	}
	return lamb;

}


/// @brief print out the spin-obit coupling potential
/// 
/// @param s  
/// @param d 
/// @param Kvec 
/// @return 
double complex *HSO ( Lattice *s, Eigen *d, double *k_vec)
{
    double q[3], cross[3];
    double k_pg[3],kp_pg[3];
	double k_mod;
    double complex vso_tmp;
	double complex vloc_tmp;
	
    int i,j,k,n,m;

	Atom **atoms= s->a_set->atom_array;
	int n_atom=s->a_set->n_atoms;
    int     NG =d->NG;
	double Tpiba_sq = (2*PI/s->A0)*(2*PI/s->A0);
    double  **G_vec = d->G_vec;
    double complex *V_sp  = SafeCalloc(NG*NG*4, sizeof(double complex));
	double *v_s = lambda(s);

	for(i=0;i<NG;i++)
	{
		for(j=0;j<=i;j++)
		{
			vso_tmp = 0;
			// compute normalized k+G
			for(k=0;k<3;k++) k_pg[k]= k_vec[k]+G_vec[i][k];
			// compute normalized k+G'
			for(k=0;k<3;k++) kp_pg[k]=k_vec[k]+G_vec[j][k];

			for(k=0;k<3;k++) q[k]=G_vec[i][k]-G_vec[j][k];

			vloc_tmp = PotentialMix(s,q);
			//Calculate the value of (G+k) X (G'+k)
			cross[0] = (k_pg[1]*kp_pg[2]-k_pg[2]*kp_pg[1])/(Tpiba_sq);
			cross[1] = (k_pg[2]*kp_pg[0]-k_pg[0]*kp_pg[2])/(Tpiba_sq);
			cross[2] = (k_pg[0]*kp_pg[1]-k_pg[1]*kp_pg[0])/(Tpiba_sq);

			for(k=0;k<n_atom;k++)
			{
				m=atoms[k]->spe;
				vso_tmp = vso_tmp + v_s[m]*cexp(-I*Dot(q,atoms[k]->tau));
			}
			vso_tmp /= n_atom;

			//multiply the spin-orbit coupling potential with the phase factor
			double complex v_test1 = 0.00001*(-I*cross[2]);
			double complex v_test2 = 0.00001*(-cross[1]-I*cross[0]);
			//print out the complex number v_test
			//printf("%f:  %f + %fi\n",cross[2], creal(v_test), cimag(v_test));

			// V_sp[(i*2)*2*NG+2*j] 	= vloc_tmp+v_test1;
			// V_sp[ (2*i+1)*2*NG+2*j+1] 	= vloc_tmp+conj(v_test1);
			// V_sp[(i*2)*2*NG+2*j+1 ]		= v_test2;   //Block2: s= 1/2,s'=-1/2
			// V_sp[ (2*i+1)*2*NG+2*j] 	= conj(v_test2);    //Block3: s=-1/2,s'= 1/2
		
			V_sp[(i*2)*2*NG+2*j] 		= vloc_tmp+vso_tmp*(-I*cross[2]);         //Block1: s= 1/2,s'= 1/2
			V_sp[(i*2)*2*NG+2*j+1 ]		= vso_tmp*(-cross[1]-I*cross[0]);   //Block2: s= 1/2,s'=-1/2
			V_sp[ (2*i+1)*2*NG+2*j] 	= vso_tmp*(cross[1]-I*cross[0]);    //Block3: s=-1/2,s'= 1/2
			V_sp[ (2*i+1)*2*NG+2*j+1] 	= vloc_tmp+vso_tmp*(I*cross[2]);          //Block4: s=-1/2,s'=-1/2
		}
	}
	return V_sp;
}
