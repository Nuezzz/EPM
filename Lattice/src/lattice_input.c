#include <stdio.h>  
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "memory.h"
#include "util.h"
#include "reader.h"
#include "../include/atom.h"

void ReadPoscar	(Lattice *s, const char *title)
{
	int			i,j;
	double		lc;
	char 		*filename;

	filename = Fullpath
	
	if(!FileExists(filename))
    {
        free(filename);
        printf("No .region file found.\n Current mesh is regionless.\n");
        fflush(stdout);
        return;
    }

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