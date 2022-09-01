/*****************************************************************************/
//gcc -O3 -lm -lrt -lfftw3 -o transform Potential_transform.c && time ./transform


#include <stdio.h>  
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <time.h>

#define QPI 0.785398163
#define TPI 6.283185307


#define SOURCE "LOCPOT-InAs-RSopt"
#define CAT_DATA "cation_potential_ZB_1.dat"
#define AN_DATA "anion_potential_ZB_1.dat"
#define MAP	"V(G).dat"

#define WIDTH	5

/* Create abstract data type for vector */
typedef struct {
	double q[4];
	double pan;
	double pcat;
} pot_q;


void build_la	(FILE *fp);
//void finput		(FILE* fp, int width, int* height, double* in);    /* read the data */ 
void finput (FILE *fp, double* dst);

void foutput	( pot_q* data, int width, int height );	 /* write the data */  

void ffill		( fftw_complex* data );			 /* fill the symetry */
void fcheck		( double* in, double* check, int a, int b, int c);
void fconvert	( fftw_complex* out, pot_q* data);

void build_re	();
void fsort	( pot_q* data);

double	a[3][3];
double	b[3][3];

int		NX,NY,NZ;


int main()  
{  
    int 		HEIGHT;                            
    double*		in;
    pot_q*		prt_p; 
	FILE		*fp;
    fftw_complex*	out;
    fftw_plan		p1;
//    fftw_plan		p2;



	fp = fopen( SOURCE, "r") ;
	if(fp==NULL){
		printf("Fail to open source file");
		exit(1);
	}
	else	printf("Source File Opened!\n");

	build_la(fp); 	// read title and lattice information
	build_re();

	HEIGHT	=	NX*NY*NZ/WIDTH + 1;		// calculate # of columns in file
/******************Allocate input and output array********************/

    in		= 	(double *) fftw_malloc (sizeof(double)*NX*NY*NZ);
    prt_p	= 	(pot_q *) malloc (sizeof(pot_q)*NX*NY*NZ);
    out		= 	(fftw_complex *) fftw_malloc (sizeof(fftw_complex)*NX*NY*NZ);
    out		= 	(fftw_complex *) fftw_malloc (sizeof(fftw_complex)*NX*NY*(NZ/2+1));

    p1 = fftw_plan_dft_r2c_3d(NX, NY, NZ, in, out,
	                     FFTW_MEASURE);         

//    finput	(fp, WIDTH, HEIGHT,in);
    finput	(fp, in);
    fftw_execute(p1);   

    ffill	(out);
    fconvert	(out,prt_p);
    foutput	(prt_p, WIDTH, HEIGHT);
/********************** Check if FFT works correctly ****************************/
/*    p2 = fftw_plan_dft_c2r_3d(NX,NY,NZ,out,check,
				 FFTW_ESTIMATE);
    fftw_execute(p2);    
    fcheck(in,check,3,3,3); */


/************************free all space*****************************/
    fftw_destroy_plan(p1);
//    fftw_destroy_plan(p2);    

    fftw_free(in);
    fftw_free(out);
    fftw_free(prt_p);
    if ( fclose( fp ))  printf ("error when close real_potential.dat");
    return 0;  
}  



/**************read data**************/
void finput (FILE *fp, double* dst)
//void finput (FILE* fp, int width, int height, double* dst)
{
	int			i,j;
	double*		indx = dst;

	for (i=0; i < NX*NY*NZ ; i++){
		fscanf(fp,"%lf",indx);
		indx += 1;
		}
	
	/*for(i = 0; i < height; i++){
		for (j = 0 ; j < width; j++){
            if ((i*width+j)<(NX*NY*NZ)){
				indx = dst + i*width + j;
			//	fscanf(fp,"%lf", &tmp);
			//	*indx = tmp;
				fscanf(fp,"%lf",indx);
				} // end if
			}// end i
	}// end j*/
	

	
}

 

/****************read the structure's direct lattice vector*********************/
void build_la	(FILE *fp){
	int			i,j;
	double		lc;
	char 		title[60];
	
//	fscanf(fp, "%[^\n]", title);
	if(fgets(title,60,fp)==NULL){
		puts(title);
	};
	printf("%s\n",title);

	fscanf(fp,"%lf",&lc);
   
   	for(i = 0; i < 3; i++){
		for (j = 0 ; j < 3; j++){
				fscanf(fp,"%lf",&(a[i][j]));
				a[i][j] *= lc;
			}// end j
	}//end i
    

    for(i = 0; i < 3; i++){
		printf("a_%d:\n",i+1);
		printf("( %lf, %lf, %lf)\n", a[i][0],a[i][1],a[i][2]);
		}

	fgets(title,60,fp);	
	fgets(title,60,fp);	
	fgets(title,60,fp);	
	fgets(title,60,fp);	
	fgets(title,60,fp);	
	fgets(title,60,fp);	

	fscanf(fp,"%d   %d   %d", &NX, &NY, &NZ);
	printf("The real space potential has grid of %d X %d X %d\n",NX,NY,NZ);
}
	
/******************************************************/
void foutput ( pot_q* data, int width, int height)
{  
	int        	i,j;
	double		tmp_cat,tmp_an;
	pot_q*		tmp = data;


	
	FILE      	*fp1;
	FILE      	*fp2;
	FILE      	*fp3;

	fp1 = fopen( AN_DATA, "w" );
	fp2 = fopen( CAT_DATA, "w" );
	fp3 = fopen( MAP, "a" );
	
	/*fprintf (fp1, "Anion Potential \n %12.8lf\n", 1.000000) ;                       
	fprintf (fp2, "Cation Potential\n %12.8lf\n", 1.000000) ;
	fprintf (fp3, "Mag_q,	Anion_p,	Cation_p\n");

	for(i = 0; i < 3; i++){
		for (j = 0 ; j < 3; j++){
			fprintf (fp1, "%12.8lf	", b[i][j]) ;                       
			fprintf (fp2, "%12.8lf	", b[i][j]) ;
			}
		fprintf(fp1,"\n");
		fprintf(fp2,"\n");
	}

	for(i=0;	i<height;	i++){
		for(j=0;	j<width;	j++){
			fprintf (fp1, "%12.8lf	", (tmp+i*width+j)->pan) ;
			fprintf (fp2, "%12.8lf	", (tmp+i*width+j)->pcat);
			//tmp=tmp+1;
			}
		fprintf(fp1,"\n");
		fprintf(fp2,"\n");
		}*/

	fsort(data);

	for(i=0; (tmp+i)->q[3]<10.000 ; ){
		tmp_cat = 0;
		tmp_an  = 0;
		for(j=0;(tmp+i+j)->q[3]-(tmp+i)->q[3]<0.0001;j++){
			tmp_cat += (tmp+i+j)->pcat;
			tmp_an  += (tmp+i+j)->pan ;
		}// end j

		fprintf (fp3, "%12.8lf,	", (tmp+i)->q[3]) ;
		fprintf (fp3, "%12.8lf,	", tmp_an /(double)j);
		fprintf (fp3, "%12.8lf,	", tmp_cat/(double)j);
		i+=j;

		/*fprintf (fp3, "%12.8lf,	", (tmp+i)->pan ) ;
		fprintf (fp3, "%12.8lf	", (tmp+i)->pcat) ;
		i++;*/

		fprintf (fp3,"\n") ;


	}// end i
		
	if ( fclose( fp1 ))  printf ("error when close anion_potential.dat" );
	if ( fclose( fp2 ))  printf ("error when close cation_potential.dat");
	if ( fclose( fp3 ))  printf ("error when close V(G).dat");
	 
}  
  
  
void ffill ( fftw_complex* data)
{
	int	i,j,k;
	for(i=0;	i<NX;	i++){
		for(j=0;	j<NY;	j++){
			for(k=0;	k<NZ/2;	k++){
			data[i*NY*NZ+j*NZ+NZ-k-1] =conj(data[i*NY*NZ+j*NZ+k]);
			}
		}
	}
	
}


/*********************** convert V_sym and V_anti into V_cat and V_an**************************/
void fconvert	( fftw_complex* out, pot_q* data){
     int		i,j,k;
	 double		arc;
	 pot_q*		tmp = data;
     fftw_complex*	ps	= out;
		

	for(i=0;i<NX;i++){
		for(j=0;j<NY;j++){
			for(k=0;k<NZ;k++){

				if((i+j+k)%2){
					arc = QPI*(i+j+k);
				
					tmp->q[0] = (double)i*b[0][0]+(double)j*b[1][0]+(double)k*b[2][0];
								
					tmp->q[1] = (double)i*b[0][1]+(double)j*b[1][1]+(double)k*b[2][1];
				
					tmp->q[2] = (double)i*b[0][2]+(double)j*b[1][2]+(double)k*b[2][2];
				
					tmp->q[3] = sqrt(tmp->q[0]*tmp->q[0]+tmp->q[1]*tmp->q[1]+tmp->q[2]*tmp->q[2]);


/*			v_an [i*NZ*NY+j*NZ+k]=  cimag(ps[i*NZ*NY+j*NZ+k])
						/2.0/NX/NY/NZ;
			v_cat[i*NZ*NY+j*NZ+k]=  creal(ps[i*NZ*NY+j*NZ+k])
						/2.0/NX/NY/NZ;*/

					tmp->pan	=	( creal(ps[i*NZ*NY+j*NZ+k])/cos(arc)
									+ cimag(ps[i*NZ*NY+j*NZ+k])/sin(arc) )
									/2.0/NX/NY/NZ;
					tmp->pcat	=	( creal(ps[i*NZ*NY+j*NZ+k])/cos(arc)
									- cimag(ps[i*NZ*NY+j*NZ+k])/sin(arc) )
									/2.0/NX/NY/NZ;
					}// end if
				else{
					tmp->q[3]=9999999999;	// asign an absurd value when it's unsolvable				
					}// end else
				tmp += 1;

			}// end k
		}// end j
		printf("V-r(%d,%d,%d): %lf\n",i,j,k,creal(ps[i*NZ*NY+j*NZ+k]));
		printf("V-i(%d,%d,%d): %lf\n",i,j,k,cimag(ps[i*NZ*NY+j*NZ+k]));
	}// end i    
}

void fcheck	( double* in, double* check, int a, int b, int c){
	int i,j,k;
	
	printf("X	Y	Z	Origin	Transformed\n");
    
	for(i=0; i<a; i++){
		for(j=0; j<b; j++){
			for(k=0; k<c; k++){
			printf("%d	%d	%d	%12.8lf	%12.8lf\n", i, j, k,
				in[i*NY*NZ+j*NZ+k],check[i*NY*NZ+j*NZ+k]/(NX*NY*NZ));
			}
		}
	}  

}


	
/********************Find reciprocal lattice vector**********************/
void build_re	(){
	int i,j;
	double vol;
	
	b[0][0] = a[1][1]*a[2][2] - a[1][2]*a[2][1];
	b[0][1] = a[1][2]*a[2][0] - a[1][0]*a[2][2];
	b[0][2] = a[1][0]*a[2][1] - a[1][1]*a[2][0];
	
	b[1][0] = a[2][1]*a[0][2] - a[2][2]*a[0][1];
	b[1][1] = a[2][2]*a[0][0] - a[2][0]*a[0][2];
	b[1][2] = a[2][0]*a[0][1] - a[2][1]*a[0][0];
	
	b[2][0] = a[0][1]*a[1][2] - a[0][1]*a[1][1];
	b[2][1] = a[0][2]*a[1][0] - a[0][0]*a[1][2];
	b[2][2] = a[0][0]*a[1][1] - a[0][1]*a[1][0];

	vol	=	abs(a[1][1]*b[1][1]+a[1][2]*b[1][2]+a[1][0]*b[1][0]);

	for (i=0; i<3; i++){
		for (j=0; j<3; j++){
			b[i][j]=b[i][j]*TPI/vol;  
			} // end j
		}// end i
	for(i = 0; i < 3; i++){
		printf("b_%d:\n",i+1);
		printf("( %lf, %lf, %lf)\n", b[i][0],b[i][1],b[i][2]);
		}

	}


/************************ Sort potential by |q| **************************/

void fsort      ( pot_q* data){
        int i,j;
        double tmp_an, tmp_cat;
        double tmp_q[4];

        for (i=0; i < NX*NY*NZ ; i++){
                for (j=i; j < NX*NY*NZ ; j++){
                        if((data+i)->q[3] > (data+j)->q[3]){
                                tmp_an  = (data+i)->pan;
                                tmp_cat = (data+i)->pcat;
                                tmp_q[0]= (data+i)->q[0];
                                tmp_q[1]= (data+i)->q[1];
                                tmp_q[2]= (data+i)->q[2];
                                tmp_q[3]= (data+i)->q[3];


                                (data+i)->pan  = (data+j)->pan;
                                (data+i)->pcat = (data+j)->pcat;
                                (data+i)->q[0] = (data+j)->q[0];
                                (data+i)->q[1] = (data+j)->q[1];
                                (data+i)->q[2] = (data+j)->q[2];
                                (data+i)->q[3] = (data+j)->q[3];

                                (data+j)->pan  = tmp_an;
                                (data+j)->pcat = tmp_cat;
                                (data+j)->q[0] = tmp_q[0];
                                (data+j)->q[1] = tmp_q[1];
                                (data+j)->q[2] = tmp_q[2];
                                (data+j)->q[3] = tmp_q[3];
                        } // end if
                } // end j
        } // end i                      

}
