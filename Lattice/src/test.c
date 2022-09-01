//gcc -O3 -lm -lrt -lfftw3 -o STRUCTURE_FACTOR.EXE test.c &&  ./STRUCTURE_FACTOR.EXE

#include <stdio.h>  
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>

#define TPI 6.283185307179586477

#define SOURCE1 "LOCPOT-InSb-12-supercell_0.1243tau"
#define SOURCE2 "LOCPOT-SbIn-12-supercell_0.1243tau"

#define CAT_DATA "V+.dat"
#define AN_DATA "V-.dat"
#define MAP	"V(G).dat"

#define WIDTH	5
/* Create abstract data type for vector */




/* Function Declarations */
int		BUILD_LA					(FILE *fp);
void	BUILD_RPOT					(FILE *fp1, FILE *fp2, fftw_complex* dst_p, fftw_complex* dst_m);
void	BUILD_RE					();
void	READ_TAU					(FILE *fp, double* tau, int N);

void 	CONVERT						( fftw_complex* out_s, fftw_complex* out_m, pot_q* data,  double* tau, int N);
void 	FOUT 						( pot_q* data, int width);

/******************  Subfunctions **********************/
void	STRUCTURE_FACTOR			(double* G, double* tau, int N, double* sum, double* min);
//void 	SORT	   					( pot_q* data);
void 	quickSort					(pot_q arr[], int low, int high);
double	LINEAR_SOLVE				( fftw_complex V, double* fa);
double	DOT							(double* a, double* b, int n);

/* Global variables */
int		NX,NY,NZ;
double	a[3][3];
double	b[3][3];
double	vol;




int main()  
{  
    int 			i,N;  
	double			*tau;

	fftw_complex	*sum, *min;
    fftw_complex	*out_s, *out_m;
    fftw_plan		p1, p2;
    
    FILE			*fp1, *fp2;                       
    pot_q*			prt_p; 
      
      
/* Open two source files*/
    fp1 = fopen( SOURCE1, "r") ;
    fp2 = fopen( SOURCE2, "r") ;
	if(fp1==NULL || fp2==NULL){
		printf("Fail to open source file");
		exit(1);
	}
	else	printf("Source File Opened!\n");
	

/* Read structure parameters and move two File pointers to the energy section*/
	N = BUILD_LA(fp1);
	N = BUILD_LA(fp2);
	
	for(i = 0; i < 3; i++){
		printf("a_%d:\n",i+1);
		printf("( %lf, %lf, %lf)\n", a[i][0],a[i][1],a[i][2]);
		}
	printf("\n**************************************\n\n");
	
	tau = (double*) malloc ( 3*N*sizeof(double));
	
    READ_TAU(fp1,tau,N);
    READ_TAU(fp2,tau,N);
    for(i = 0; i < N; i++){
		printf("tau_%d:	",i+1);
		printf("( %lf, %lf, %lf)\n", tau[3*i],tau[3*i+1],tau[3*i+2]);
		}
	printf("\n**************************************\n\n");

	
    BUILD_RE();

    		
	printf("*******************************************\n\n");
	printf("The real space potential has grid of %d X %d X %d\n",NX,NY,NZ);	
	

/* Allocate input and output arrays */

    prt_p	= 	(pot_q *) malloc (sizeof(pot_q)*NX*NY*NZ);
    sum		= 	(fftw_complex *) fftw_malloc (sizeof(fftw_complex)*NX*NY*NZ);   
    min		= 	(fftw_complex *) fftw_malloc (sizeof(fftw_complex)*NX*NY*NZ);
    out_s	= 	(fftw_complex *) fftw_malloc (sizeof(fftw_complex)*NX*NY*NZ);
    out_m	= 	(fftw_complex *) fftw_malloc (sizeof(fftw_complex)*NX*NY*NZ);
 
/* FFT potential map into reciprocal space */    
    p1 = fftw_plan_dft_3d(NX, NY, NZ, sum, out_s,FFTW_FORWARD,FFTW_MEASURE); 
    p2 = fftw_plan_dft_3d(NX, NY, NZ, min, out_m,FFTW_FORWARD,FFTW_MEASURE); 
	                     
	BUILD_RPOT( fp1, fp2, sum, min) ; 
	
	fftw_execute(p1);  
	fftw_execute(p2);  
	
/* Convert the potential and output the results */   
	CONVERT(out_s, out_m, prt_p, tau, N); 
	 
	FOUT( prt_p, WIDTH);
	
/* Free all the spece and pouinters */	
	fftw_free(sum);
	fftw_free(min);
    fftw_free(out_s);
    fftw_free(out_m);
    fftw_free(prt_p);
    fftw_free(tau);
    
    if ( fclose( fp1 ))  printf ("error when close %s", SOURCE1);
    if ( fclose( fp2 ))  printf ("error when close %s", SOURCE2);
    return 0; 
 }
 
 
 
 
 
 
/********************************************************************
 * Funtion															*
 * 		calculate structure factor related quantities				*
 * Input															*
 * 		G:		Reciprocal vector with rank 3						*
 * 		tau: 	Position vector of all atoms						*
 * 		N:		Number of atoms inside the unit cell				*
 * Output:															*
 * 		Alpha_p	=	Sum of cos(G*tau[n]) for N atoms				*
 *		Beta_p	=	Sum of sin(G*tau[n]) for N atoms				*
 * 		Alpha_m	=	Sum of (-1)^n *cos(G*tau[n]) for N atoms		*
 * 		Beta_m	=	Sum of (-1)^n *sin(G*tau[n]) for N atoms		*
 ********************************************************************/
void STRUCTURE_FACTOR(double* G, double* tau, int N, double* sum, double* min){
	int 	i;
	double*	tmp;
	double	complex str_1, str_2;
	double	complex element;
	
	tmp			=	tau;
	str_1		=	0;
	str_2		=	0;
	for(i = 0; i < N/2; i++){
		element	=	cexp(-I*DOT(G,tmp,3));
		str_1	+=	element; // accumulator for exp(-I*G*tau)
		str_2	+=	element; // accumulator for (-1)^n * exp(-I*G*tau)
		tmp +=	3;

		element	=	cexp(-I*DOT(G,tmp,3));
		str_1	+=	element;
		str_2	-=	element;
		tmp	+=	3;
		}
		
	sum[0]  = creal( str_1);
	sum[1]	= cimag( str_1);
	
	min[0]	= creal( str_2);
	min[1]	= cimag( str_2);
	}





/********************************************************************
 * Funtion															*
 * 		Read the direct lattice structure							*
 * Input															*
 * 		fp:		File pointer to the input file						*
 * 																	*
 *  Output															*
 * 		a[3]:	Direct lattice vector passing by global variable	*
 * 		N:		Number of atoms inside the unit cell				*
 ********************************************************************/
int BUILD_LA	(FILE *fp ){
	int			i,j,N;
	int			m=0;		// m is # of atom A
	int			n=0;		// n is # of atom B
	double		lc;			// lc is the lattice constant
	char 		title[60],A[2],B[2];
	

	if(fgets(title,60,fp)==NULL){
		puts(title);
	};
	printf("%s\n",title);
	
/******************	Read the lattice vector a[3] *********************/
	fscanf(fp,"%lf",&lc);
   	for(i = 0; i < 3; i++){
		for (j = 0 ; j < 3; j++){
				fscanf(fp,"%lf",&(a[i][j]));
				a[i][j] *= lc;
			}// end j
	}//end i
    


	fscanf(fp,"%s",A);
	fscanf(fp,"%s",B);
//	printf("Two types of atom are: %s and %s\n",A,B);	
	
/******************	Read the # of atoms ******************************/	
	fscanf(fp,"%d",&m);	// # of atoms A
	fscanf(fp,"%d",&n);	// # of atoms B
	N=m+n;

	return N;
}
	



/********************************************************************
 *  Funtion															*
 * 		Read the atoms' position vectors							*
 *  Input															*
 * 		fp:		File pointer to the input file						*
 * 		tau:	Pointer for atoms position array					*
 * 		n:		# of atoms inside the cell							*
 * 																	*
 *  Output															*
 * 		tau:	Pointer for atoms position array					*
 ********************************************************************/	

void	READ_TAU (FILE *fp, double* tau, int N){
	int			i,j;
	double*		tmp;
	double		frac[3];		// fraction coordinates for one atom
	char 		title[10];		// type of coordination
	/******************	Read the atoms position tau **********************/

	fscanf (fp,"%s",title);				// get type of coordination of atoms position
	printf("Type of coordination is %s\n\n",title);
	for(i = 0; i < N; i++){
		for (j = 0 ; j < 3; j++){
			tmp=tau+3*i+j;
			fscanf(fp,"%lf",tmp);
		}
	}
	
	if(strcmp(title, "Direct") == 0){	// exicute if when "title == Direct"
		for(i = 0; i < N; i++){
			for (j = 0 ; j < 3; j++){
				frac[j] = tau[3*i+j];
			}
			for (j = 0 ; j < 3; j++){
				tau[3*i+j] = frac[0]*a[0][j]+frac[1]*a[1][j]+frac[2]*a[2][j] ;
			}
		}
	}
	
	printf("position vector tau of %d atoms has been get\n",N);	

		
	fgets(title,10,fp);	
	

	fscanf(fp,"%d   %d   %d", &NX, &NY, &NZ);
	}	
	
	
	
	
	
/********************************************************************
 * Function															*
 * 		Read the real space potential map							*
 * Input															*
 * 		fp1:		File pointer to the input file					*
 * 		fp2:		File pointer to the atom exchanged input file	*
 * 		dst_p:	Pointer to symetric potential arrays				*
 * 		dst_m:	Pointer to anti-sym potential arrays				*
 * Outpot															*
 *		The conbination of two real space potential map				*
 ********************************************************************/
void BUILD_RPOT (FILE *fp1, FILE *fp2, fftw_complex* dst_p, fftw_complex* dst_m){
	int			i;
	double		v1,v2;
	fftw_complex*		pp = dst_p;
	fftw_complex*		pm = dst_m;

	for (i=0; i < NX*NY*NZ ; i++){
		fscanf(fp1,"%lf",&v1);	// Read the potential at i_th grid in crystal
		fscanf(fp2,"%lf",&v2);	// Read the potential at i_th grid of exchanged atom
		*pp = v1+v2;
		*pm = v1-v2;
		pp += 1;
		pm += 1;
		}
	printf("V+(r) and V-(r) has been builded\n");
	printf("*******************************************\n\n");
	
}

	
/********************************************************************
 *  Funtion															*
 * 		Build the reciprocal lattice structure						*
 *  Input															*
 * 																	*
 *  Output															*
 * 		b[3]:	Reciprocal lattice vector passing by global variable*
 ********************************************************************/
void BUILD_RE	(){
	int i,j;

	
	b[0][0] = a[1][1]*a[2][2] - a[1][2]*a[2][1];
	b[0][1] = a[1][2]*a[2][0] - a[1][0]*a[2][2];
	b[0][2] = a[1][0]*a[2][1] - a[1][1]*a[2][0];
	
	b[1][0] = a[2][1]*a[0][2] - a[2][2]*a[0][1];
	b[1][1] = a[2][2]*a[0][0] - a[2][0]*a[0][2];
	b[1][2] = a[2][0]*a[0][1] - a[2][1]*a[0][0];
	
	b[2][0] = a[0][1]*a[1][2] - a[0][1]*a[1][1];
	b[2][1] = a[0][2]*a[1][0] - a[0][0]*a[1][2];
	b[2][2] = a[0][0]*a[1][1] - a[0][1]*a[1][0];

	vol	=	fabs(a[1][1]*b[1][1]+a[1][2]*b[1][2]+a[1][0]*b[1][0]);

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




/********************************************************************
 *  Funtion															*
 * 		calculate the dot production of two vectors					*
 *  Input															*
 * 		a:		pointer of the first vector							*
 * 		b:		pointer of the second vector						*
 * 		n:		rank of the vectors									*
 *  Output															*
 * 		The inner product of two vectors							*
 ********************************************************************/
double DOT( double* a, double* b, int n){
	double	r=0;
	int		i;
	for(i=0; i<n; i++){
		r += (*(a+i))* (*(b+i));
		}
	return r;
	}



/********************************************************************
 *  Funtion															*
 * 		calculate the dot production of two vectors					*
 *  Input															*
 * 		a:		pointer of the first vector							*
 * 		b:		pointer of the second vector						*
 * 		n:		rank of the vectors									*
 *  Output															*
 * 		The inner product of two vectors							*
 ********************************************************************/

void CONVERT	( fftw_complex* out_s, fftw_complex* out_m, pot_q* data,  double* tau, int N){
     int			i,j,k,ii;
	 double			G[3], factor_s[2], factor_m[2];
	 pot_q*			tmp = data;
     fftw_complex*	ps1	= out_s;
	 fftw_complex*	ps2	= out_m;
			
	for(i=0;i<NX;i++){
		for(j=0;j<NY;j++){
			for(k=0;k<NZ;k++){
				for(ii = 0; ii < 3; ii++){
					G[ii]=i*b[0][ii]+j*b[1][ii]+k*b[2][ii];
					}
				
				tmp->q[0] = G[0];	
				tmp->q[1] = G[1];
				tmp->q[2] = G[2];
				tmp->q[3] = sqrt(G[0]*G[0]+G[1]*G[1]+G[2]*G[2]);

				STRUCTURE_FACTOR(G, tau, N, factor_s, factor_m);

				tmp->p_s	=	LINEAR_SOLVE(ps1[i*NZ*NY+j*NZ+k], factor_s);
				tmp->p_m	=	LINEAR_SOLVE(ps2[i*NZ*NY+j*NZ+k], factor_m);

				tmp += 1;

			}// end k
		}// end j
	}// end i 
	
	
	/*for(i=0;i<10;i++){
		for(j=0;j<2;j++){
			for(k=0;k<2;k++){	   
				printf("(%d,%d,%d)	%lf	:	",i,j,k,data[i*NY*NZ+j*NZ+k].q[3]);
				printf(" %lf,	%lf\n",data[i*NY*NZ+j*NZ+k].p_s,data[i*NY*NZ+j*NZ+k].p_m);
				}
			}
		}*/
	printf("V+(q) and V-(q) has been converted\n");
	printf("*******************************************\n\n");
}




/********************************************************************
 *  Funtion															*
 * 		calculate V+(G) and V-(G) with the result of corresponding	*
 * 		FFT.														*
 * 		fa[2]= {Alpha_p, Beta_p} for V+				 				*
 * 		fa[2]= {Alpha_m, Beta_m} for V-				 				*
 *  Input															*
 * 		V:		+/- Local potential for choosing G vector			*
 * 		fa:		pointer for linear equation set matrix elements		*
 *  Output															*
 * 		Corresponding V+(G) or V-(G)								*
 ********************************************************************/

double LINEAR_SOLVE	( fftw_complex V, double* fa){
	 double		res;
	 double		det;

	det = fa[0]*fa[0]+fa[1]*fa[1];
	res	= vol*(fa[0]*creal(V)+fa[0]*cimag(V))/det;
	return res;
	}




/********************************************************************
 *  Funtion															*
 * 		Sort potential by |q|										*
 *  Input															*
 * 		data:	pointer for structure array 			 			*
 *  Output															*
 * 		Descending order of data structure							*
 ********************************************************************/
/*void SORT      ( pot_q *data){
        int i,j;
//        double tmp_s, tmp_m;
//        double tmp_q[4];
		pot_q tmp;
		printf("Begin sorting\n");
		
        for (i=0; i < NX*NY*NZ ; i++){
                for (j=i; j < NX*NY*NZ ; j++){
                        if(data[i].q[3] > data[j].q[3]){
 //                               tmp_s 	= (data+i)->p_s;
 //                               tmp_m	= (data+i)->p_m;
 //                               tmp_q[0]= (data+i)->q[0];
 //                               tmp_q[1]= (data+i)->q[1];
 //                               tmp_q[2]= (data+i)->q[2];
 //                               tmp_q[3]= (data+i)->q[3];
								tmp=data[i];
								
//                                (data+i)->p_s	= (data+j)->p_s;
//                                (data+i)->p_m	= (data+j)->p_m;
//                                (data+i)->q[0]	= (data+j)->q[0];
//                                (data+i)->q[1]	= (data+j)->q[1];
//                                (data+i)->q[2]	= (data+j)->q[2];
//                                (data+i)->q[3]	= (data+j)->q[3];
								data[i] = data[j];

//                                (data+j)->p_s	= tmp_s;
//                                (data+j)->p_m	= tmp_m;
//                                (data+j)->q[0]	= tmp_q[0];
//                                (data+j)->q[1]	= tmp_q[1];
//                                (data+j)->q[2]	= tmp_q[2];
//                                (data+j)->q[3]	= tmp_q[3];
								data[j] = tmp;
                        } // end if
                } // end j
        } // end i                      
	printf("Input file has been sorted\n");
	printf("*******************************************\n\n");
}

*/



/********************************************************************
 *  Funtion															*
 * 		Print V+(G)	and V-(G) potential map							*
 * 		Print V+(|G|) and V-(|G|)
 *  Input															*
 * 		data:	pointer for structure array 			 			*
 * 		width:	number of points for single line		 			*
 * 		height:	number of lines in the data section		 			*
 *  Output															*
 * 		Descending order of data structure							*
 ********************************************************************/
void FOUT ( pot_q* data, int width){  
	int        	i,j;
	double		tmp_s, tmp_m;
	pot_q*		tmp		= 	data;
//	int			height	=	NX*NY*NZ/WIDTH + 1;		// calculate # of columns in file

//	printf("strat printing results \n");
	
	FILE      	*fp1;
	FILE      	*fp2;
	FILE      	*fp3;

	fp1 = fopen( AN_DATA, "w" );
	fp2 = fopen( CAT_DATA, "w" );
	fp3 = fopen( MAP, "a" );
	
	if(fp1==NULL || fp2==NULL || fp3==NULL){
		printf("Fail to open source file");
		exit(1);
	}
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
			fprintf (fp1, "%12.8lf	", (tmp+i*width+j)->p_s) ;
			fprintf (fp2, "%12.8lf	", (tmp+i*width+j)->p_m);
			//tmp=tmp+1;
			}
		fprintf(fp1,"\n");
		fprintf(fp2,"\n");
		}*/


	quickSort(tmp,0,NX*NY*NZ-1);

	for(i=0; tmp[i].q[3]<20.000 ; ){
		tmp_s = 0;
		tmp_m  = 0;
		for(j=0;(tmp+i+j)->q[3]-(tmp+i)->q[3]<0.0001;j++){
			tmp_s += (tmp+i+j)->p_s;
			tmp_m  += (tmp+i+j)->p_m ;
		}// end j

		fprintf (fp3, "%12.8lf,	", (tmp+i)->q[3]) ;
		fprintf (fp3, "%12.8lf,	", tmp_s /(double)j);
		fprintf (fp3, "%12.8lf\n", tmp_m/(double)j);
		i+=j;

		/*fprintf (fp3, "%12.8lf,	", (tmp+i)->pan ) ;
		fprintf (fp3, "%12.8lf	", (tmp+i)->pcat) ;
		i++;

		fprintf (fp3,"\n") ;*/


	}// end i
		
	if ( fclose( fp1 ))  printf ("error when close anion_potential.dat" );
	if ( fclose( fp2 ))  printf ("error when close cation_potential.dat");
	if ( fclose( fp3 ))  printf ("error when close V(G).dat");
	 
}  
 

// A utility function to swap two elements 
void swap(pot_q* a, pot_q* b) 
{ 
    pot_q t = *a; 
    *a = *b; 
    *b = t; 
} 
  
/* This function takes last element as pivot, places 
   the pivot element at its correct position in sorted 
    array, and places all smaller (smaller than pivot) 
   to left of pivot and all greater elements to right 
   of pivot */
int partition (pot_q arr[], int low, int high) 
{ 
    pot_q pivot = arr[high];    // pivot 
    int i = (low - 1);  // Index of smaller element 
  
    for (int j = low; j <= high- 1; j++) 
    { 
        // If current element is smaller than the pivot 
        if (arr[j].q[3] < pivot.q[3]) 
        { 
            i++;    // increment index of smaller element 
            swap(&arr[i], &arr[j]); 
        } 
    } 
    swap(&arr[i + 1], &arr[high]); 
    return (i + 1); 
} 
  
/* The main function that implements QuickSort 
 arr[] --> Array to be sorted, 
  low  --> Starting index, 
  high  --> Ending index */
  
void quickSort(pot_q arr[], int low, int high) 
{ 
    if (low < high) 
    { 
        /* pi is partitioning index, arr[p] is now 
           at right place */
        int pi = partition(arr, low, high); 
  
        // Separately sort elements before 
        // partition and after partition 
        quickSort(arr, low, pi - 1); 
        quickSort(arr, pi + 1, high); 
    } 
} 



 
