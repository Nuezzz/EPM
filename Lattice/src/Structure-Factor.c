//gcc -O3 -lm -lrt -o STRUCTURE_FACTOR.EXE test.c &&  ./STRUCTURE_FACTOR.EXE

#include <stdio.h>  
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#define TPI 6.283185307179586477
#define SOURCE "Position"

double	complex STRUCTURE_FACTOR	(double* G, double* tau, int N);
int		BUILD_LA					(FILE *fp);
void	BUILD_RE					();
void	READ_TAU					(FILE *fp, double* tau, int N);
double	DOT							( double* a, double* b, int n);

double	a[3][3];
double	b[3][3];

int main()  
{  
	double*			tau;
	double			G[3];
	double	complex S;
	int				X,Y,Z;	// integer coordination of G vecotr
    int 			i,j,n;     
    FILE			*fp;
    
    printf("Entering the intiger coodination of G:\n");
    scanf("%d, %d, %d",&X,&Y,&Z);
    
    
    fp = fopen( SOURCE, "r") ;
	if(fp==NULL){
		printf("Fail to open source file");
		exit(1);
	}
//	else	printf("Source File Opened!\n");

    n=BUILD_LA(fp);
    tau = (double*) malloc ( 3*n*sizeof(double));
    READ_TAU(fp,tau,n);
    
//    for(i = 0; i < n; i++){
//		printf("tau_%d:	",i+1);
//		printf("( %lf, %lf, %lf)\n", tau[3*i],tau[3*i+1],tau[3*i+2]);
//		}
//	printf("\n**************************************\n\n");
	
    BUILD_RE();
   	for(i = 0; i < 3; i++){
		for (j = 0 ; j < 3; j++){
			G[i]=X*b[0][i]+Y*b[1][i]+Z*b[2][i];
			}
		}

    		
	printf("\nG:	( %lf, %lf, %lf)\n\n*******************************************\n\n", G[0],G[1],G[2]);
		
    S=STRUCTURE_FACTOR (G, tau, n);
    printf("The structure factor of G (%d, %d, %d) is:\n", X,Y,Z);
    printf("(%12.8lf,%12.8lf)\n",creal(S),cimag(S));

 }
 
/******************************************************************
Funtion
		calculate structure factor
Input
		G:		Reciprocal vector with rank 3
		tau: 	Position vector of all atoms
		N:		Number of atoms inside the unit cell
Output:
		Sum of exp(-iG*tau[n]) for N atoms
 *****************************************************************/
double complex STRUCTURE_FACTOR(double* G, double* tau, int N){
	int 	i;
	double*	tmp;
	double	complex str_fctr = 0;
	for(i = 0; i < N; i++){
		tmp=tau+3*i;
		str_fctr += cexp(-I*DOT(G,tmp,3));
		printf("G*tau_%d:	%lf\n",i,DOT(G,tmp,3)/3.14159265359);///////////
		}

	return str_fctr;
	}


/*****************************************************************
 Funtion
		Read the direct lattice structure
 Input
		fp:		File pointer to the input file
		
 Output
		a[3]:	Direct lattice vector passing by global variable
		N:		Number of atoms inside the unit cell
******************************************************************/
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
    

    for(i = 0; i < 3; i++){
		printf("a_%d:\n",i+1);
		printf("( %lf, %lf, %lf)\n", a[i][0],a[i][1],a[i][2]);
		}
	printf("\n**************************************\n\n");
	fscanf(fp,"%s",A);
	fscanf(fp,"%s",B);
	printf("Two types of atom are: %s and %s\n",A,B);	
	
/******************	Read the # of atoms ******************************/	
	fscanf(fp,"%d",&m);	// # of atoms A
	fscanf(fp,"%d",&n);	// # of atoms B
	N=m+n;

	return N;
}
	


/*****************************************************************
 Funtion
		Read the atoms' position vectors
 Input
		fp:		File pointer to the input file
		tau:	Pointer for atoms position array
		n:		# of atoms inside the cell
		
 Output
		tau:	Pointer for atoms position array
******************************************************************/	

void	READ_TAU (FILE *fp, double* tau, int N){
	int			i,j;
	double*		tmp;
	double		frac[3];		// fraction coordinates for one atom
	char 		title[10];		// type of coordination
	/******************	Read the atoms position tau **********************/

	fscanf (fp,"%s",title);				// get type of coordination of atoms position
//	printf("Type of coordination is %s\n\n",title);
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
	
//	printf("position vector tau of %d atoms has been get\n",N);	

		
	fgets(title,10,fp);	
	

//	fscanf(fp,"%d   %d   %d", &NX, &NY, &NZ);
//	printf("The real space potential has grid of %d X %d X %d\n",NX,NY,NZ);
	}	
/*****************************************************************
 Funtion
		Build the reciprocal lattice structure
 Input

		
 Output
		b[3]:	Reciprocal lattice vector passing by global variable

******************************************************************/
void BUILD_RE	(){
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


/*****************************************************************
 Funtion
		calculate the dot production of two vectors
 Input
		a:		pointer of the first vector
		b:		pointer of the second vector
		n:		rank of the vectors
		
 Output
		The inner product of two vectors

******************************************************************/
double DOT( double* a, double* b, int n){
	double	r=0;
	int		i;
	for(i=0; i<n; i++){
		r += (*(a+i))* (*(b+i));
		}
	return r;
	}
