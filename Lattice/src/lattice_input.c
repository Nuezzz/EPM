#include <stdio.h>  
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include "memory.h"
#include "util.h"
#include "reader.h"
#include "../include/atom.h"

static inline void ReadLatticeVector( Lattice *s, FILE *fp, double lat_con)
{
	int	i,j;
	for(i = 0; i < 3; i++){
	for (j = 0 ; j < 3; j++){
			fscanf(fp,"%lf",&(s->a[i][j]));
			s->a[i][j] *= lat_con;
		}// end j
	}//end i
}

static inline void ReadMatInfo( Lattice *s, FILE *fp )
{
	int         i=0;
	int			mat_count = 0;
	int			storage_count = 0;
	char 		buff[INPUT_BUFSIZE];
	char		*pch;
	// count number of material name
	fgets(buff,INPUT_BUFSIZE,fp);	

	pch = strtok (buff," ,.-");
  	while (pch != NULL)
  	{
		mat_count ++;
		storage_count += strlen(pch);
    	pch = strtok (NULL, " ,.-");
  	}
	free(pch);
	s->n_spe = mat_count;
	s->mat_name = SafeCalloc(mat_count, sizeof( char*));
	s->mat_storage = SafeCalloc(storage_count, sizeof( char));
	// read again and record the names
	pch = strtok (buff," ,.-");
  	while (pch != NULL)
  	{
		s->mat_name[i++]=pch;
    	pch = strtok (NULL, " ,.-");
  	}

}

static inline void ReadAtomNum(Lattice *s, FILE *fp)
{
	s->natom_spe = SafeCalloc(s->n_spe, sizeof(int));
	for(int i=0; i < s->n_spe; i++)
	{
		fscanf(fp,"%d",&(s->natom_spe[i]));
	}
}


static inline void ReadAtomPos(Lattice *s, FILE *fp)
{
	Atom *atom;
	char *type;
	int  i,j,k;
	int	 N = s->a_set->n_atoms;
	double tmp[3];
	double lv[3][3];

	for(i=0; i<3; i++)
	{
		lv[i][0]=s->a[i][0];
		lv[i][1]=s->a[i][1];
		lv[i][2]=s->a[i][2];
	}


	/* get the type of atom position*/
	fgets(type,60,fp);

	// get type of coordination of atoms position
	if (strcmp(type, "Direct"))
	{
		for(k=0; k<s->n_spe; k++)
		{
			N = s->natom_spe[k];
				for(i = 0; i < N; i++)
				{
					atom= s->a_set->atom_array[i];

					for (j = 0 ; j < 3; j++)
					{
						fscanf(fp,"%lf",(tmp+j));
					}

					atom->spe = k;
					atom->tau[0] = lv[0][0]*tmp[0]+lv[1][0]*tmp[1]+lv[2][0]*tmp[2];
					atom->tau[1] = lv[0][1]*tmp[0]+lv[1][1]*tmp[1]+lv[2][1]*tmp[2];
					atom->tau[2] = lv[0][2]*tmp[0]+lv[1][2]*tmp[1]+lv[2][2]*tmp[2];
				}
		}
	}
	else
	{
		for(k=0; k<s->n_spe; k++)
		{
			N = s->natom_spe[k];
				for(i = 0; i < N; i++)
				{
					atom= s->a_set->atom_array[i];

					for (j = 0 ; j < 3; j++){
						fscanf(fp,"%lf",(tmp+j));
					}

					atom->spe = k;
					atom->tau[0]=tmp[0];
					atom->tau[1]=tmp[1];
					atom->tau[2]=tmp[2];
				}
		}

	}
	
}


AtomStack *AsetInitial( unsigned int Nspe, unsigned int *Natomspe)
{
	/*Sum up the total atom number*/
	AtomStack *Aset;
	int N=0;
	for(int i=0; i<Nspe; i++)
	{
		N+=Natomspe[i];
	}
	/*Allocate the space accoeding to the atom number*/
	Aset->n_atoms =N;
	Aset->atom_storage	= SafeCalloc(N, sizeof(Atom));
	Aset->atom_array  	= SafeCalloc(N, sizeof(Atom*));
	Aset->neigh_storage = SafeCalloc(Nspe, sizeof(unsigned int));
	return Aset;


}

void ReadPoscar	(Lattice *s)
{
	double		lat_con		=0.0;
	char 		*filename;
	char        *name;
	char        *path;
	char 		buff[INPUT_BUFSIZE];
	FILE		*fp;


	StringClone(path, "./");
	StrCat(&path, s->title);
	StrCat(&path, "/");

	StringClone(name, "Poscar");

	filename = Fullpath(path, name);
	
	if(!FileExists(filename))
    {
        free(filename);
        printf("No Poscar file found.\n ");
        fflush(stdout);
        return;
    }

	fp = SafeFOpen(filename,"r");


//	fscanf(fp, "%[^\n]", title);
	if(fgets(buff,INPUT_BUFSIZE,fp)==NULL){
		puts(buff);
		printf("Poscar file is corrupted\n");
		return;
	};

	printf("Lattice %s is used for band structure calculation\n",buff);

	fscanf(fp,"%lf",&lat_con);
	
	ReadLatticeVector(s, fp, lat_con);
	ReadMatInfo(s, fp );
	ReadAtomNum(s, fp);

	AsetInitial(s->a_set, s->n_spe, s->natom_spe);

	ReadAtomPos(s, fp);
	fclose(fp);

}