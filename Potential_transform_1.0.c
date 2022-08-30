/*****************************************************************************/
//gcc -O3 -lm -lrt -lrfftw -o transform Potential_transform.c && time ./transform



#include <stdio.h>  
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <fftw3.h>



#define NX 200
#define NY 200
#define NZ 200

#define WIDTH	5
#define HEIGHT	200

void finput		( int width, int height, double *in);   	 /* read the data*/ 
void foutput	( fftw_complex* out);		 /*write the data*/  

int main()  
{  


    int i,j;                            
    fftw_real in[NX][NY][2*(NZ/2+1) ];
    fftw_complex *out;
	fftw_plan p;
	p = fftw_plan_dft_r2c_3d(NX, NY, NZ,
	                        FFTW_REAL_TO_COMPLEX, 
	                        FFTW_ESTIMATE | FFTW_IN_PLACE);
	                       
/****************Alias the output pointer to the dead of input array****************/
	out = (fftw_complex) &in[0][0][0];
    
    
    finput(WIDTH,HEIGHT,in);
    fftw_execute(p, &in[0][0][0], &in[0][0][0]);
    fftw_destroy_plan(p);
    foutput(out);
 
    return 0;  
}  
/**************read data**************/
void finput( int width, int height, double *dst)
{
	int			i,j;
	double		tmp;
	double		*indx;
	FILE		*fp;
	
	fp = fopen( "real_potential.dat", "r" );
	
	for(i = 0; i < height; i++){
		for (j = 0 ; j < width; j++){
			indx = dst + i*width + j;
			fscanf(fp,"%f", &tmp);
			*indx = tmp; 
			}
	}
	
	if ( fclose( fp ))  printf ("error when close real_potential.dat");
	
}

 
void output( fftw_complex* out)  
{  
	int        i,j,k;
	FILE       *fp;
	fp = fopen( "reciprocal_potential.dat", "w" );
	fprintf(fp,"X		Y		Z		REAL		IMAGE\n");
	for(i=0;	i<NX;	i++){
		for(j=0;	j<NY;	j++){
			for(k=0;	k<(NZ/2+1);	k++){
				fprintf (fp, "%d		%d		%d		%lf		%lf \n",
						 i,j,k, creal(out[i][j][k]),cimag(out[i][j][k]));
				}
			}
		}
	
	if ( fclose( fp ))  printf ("error when close real_potential.dat");
	 
}  
  
  

