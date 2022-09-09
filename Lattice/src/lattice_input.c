#include <stdio.h>  
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include "memory.h"
#include "util.h"
#include "reader.h"
#include "atom.h"

static inline void ReadLatticeVector( Lattice *s, Reader *fr, double lat_con)
{
	int	i,j;
	char *s0;
	for(i = 0; i < 3; i++){
		

		for (j = 0 ; j < 3; j++){
			s0 = ReaderGetEntry(fr, 2+i,j);
			s->a[i][j] = atof(s0)* lat_con;
		}// end j
	}//end i
}

static inline void ReadMatInfo( Lattice *s, Reader *fr )
{
	int    	i;
	char	*s0;	


	// count number of material name

	s->n_spe = ReaderGetLineLength(fr, 5);
	s->mat_name = SafeCalloc(s->n_spe, sizeof( char*));
	// read again and record the names
  	for(i=0; i<s->n_spe; i++)
  	{
		s0=ReaderGetEntry(fr, 5,i);
		StringClone(s->mat_name[i],s0);
  	}

}

static inline void ReadAtomNum(Lattice *s, Reader *fr)
{
	char	*s0;
	s->natom_spe = SafeCalloc(s->n_spe, sizeof(int));
	for(int i=0; i < s->n_spe; i++)
	{
		s0=ReaderGetEntry(fr, 6,i);
		s->natom_spe[i] = atoi(s0);
	}
}


static inline void ReadAtomPos(Lattice *s, Reader *fr)
{
	Atom *atom;
	char *type;
	char *s0;
	int  i,j,k;
	int	 N = 0;
	double tmp[3];
	double lv[3][3];

	for(i=0; i<3; i++)
	{
		lv[i][0]=s->a[i][0];
		lv[i][1]=s->a[i][1];
		lv[i][2]=s->a[i][2];
	}


	/* get the type of atom position*/
	type = ReaderGetEntry(fr, 7, 0);

	// get type of coordination of atoms position
	if (strcmp(type, "Direct"))
	{
		for(k=0, i=0; k<s->n_spe; k++)
		{
			N += s->natom_spe[k];
				for(; i < N; i++)
				{
					atom= s->a_set->atom_array[i];

					for (j = 0 ; j < 3; j++)
					{
						s0=ReaderGetEntry(fr, 8+i,j);
						// tmp[j]=atof(s0);
						atom->tau[j]=atof(s0);
					}

					atom->spe = k;
					// atom->tau[0] = lv[0][0]*tmp[0]+lv[1][0]*tmp[1]+lv[2][0]*tmp[2];
					// atom->tau[1] = lv[0][1]*tmp[0]+lv[1][1]*tmp[1]+lv[2][1]*tmp[2];
					// atom->tau[2] = lv[0][2]*tmp[0]+lv[1][2]*tmp[1]+lv[2][2]*tmp[2];
				}
		}
		FindNeighbor (4.0, s);
		for(k=0, i=0; k<s->n_spe; k++)
		{
			N += s->natom_spe[k];
				for(; i < N; i++)
				{
					atom= s->a_set->atom_array[i];

					for (j = 0 ; j < 3; j++)
					{
						s0=ReaderGetEntry(fr, 8+i,j);
						tmp[j]=atof(s0);
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
		for(k=0, i = 0; k<s->n_spe; k++)
		{
			N += s->natom_spe[k];
				for(; i < N; i++)
				{
					atom= s->a_set->atom_array[i];

					for (j = 0 ; j < 3; j++){
						s0=ReaderGetEntry(fr, 8+i,j);
						tmp[j]=atof(s0);
					}

					atom->spe = k;
					atom->tau[0]=tmp[0];
					atom->tau[1]=tmp[1];
					atom->tau[2]=tmp[2];
				}
		}

	}
	
}


static inline AtomStack *AsetInitial( unsigned int Nspe, unsigned int *Natomspe)
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
	Aset->neigh_storage = SafeCalloc(N, sizeof(unsigned int )*Nspe);
	for(int i =0; i<N; i++)
	{
		Aset->atom_array[i]=&(Aset->atom_storage[i]);
		Aset->atom_array[i]->n_neighbor = &(Aset->neigh_storage[i]);
	}
	return Aset;


}

static inline void ReadPoscar	(Lattice *s, char *title )
{
	double		lat_con		=0.0;
	char 		*filename;
	char        *name;
	char        *path;
	char 		*s0;
	Reader		*fr;


	StringClone(path, "./");
	StrCat(&path, title);
	StrCat(&path, "/");

	StringClone(name, "Poscar");

	filename = FullPath(path, name);
	
	if(!FileExists(filename))
    {
        free(filename);
        printf("No Poscar file found.\n ");
        fflush(stdout);
        return;
    }

	fr = ReaderReadFile(filename);
	s0 = ReaderGetEntry(fr,0,0);

    
	printf("Lattice %s is used for band structure calculation\n",s0);

	s0 = ReaderGetEntry(fr,1,0);
	lat_con = atof(s0);
	
	ReadLatticeVector(s, fr, lat_con);
	ReadMatInfo(s, fr );
	ReadAtomNum(s, fr);

	s->a_set=AsetInitial( s->n_spe, s->natom_spe);

	ReadAtomPos(s, fr);

}

Lattice *LatticeInitial ( char *filename )
{
	char 	*s0;
	Lattice *L;
	L = SafeCalloc(1,sizeof(Lattice));
	ReadPoscar(L, filename);
	
	
	return L;


}