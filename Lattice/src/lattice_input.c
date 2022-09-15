#include <stdio.h>  
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include "memory.h"
#include "util.h"
#include "reader.h"
#include "constants.h"
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

/**
 * @brief Read the atom position and finding the neighbor atoms
 * 
 * @param s the lattice pointer
 * @param fr the reader pointer of the input structure file
 * @param N  total number of atoms in the lattice
 */
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
		FindNeighbor (3, s);
		for(i=0; i<N; i++)
		{

				atom= s->a_set->atom_array[i];
				for (j = 0 ; j < 3; j++)
				{
					tmp[j]=atom->tau[j];
				}
				atom->tau[0] = lv[0][0]*tmp[0]+lv[1][0]*tmp[1]+lv[2][0]*tmp[2];
				atom->tau[1] = lv[0][1]*tmp[0]+lv[1][1]*tmp[1]+lv[2][1]*tmp[2];
				atom->tau[2] = lv[0][2]*tmp[0]+lv[1][2]*tmp[1]+lv[2][2]*tmp[2];

		 }
	}
	else
	{
		printf("To utilize the neighbor finding function, please use direct position of stoms\n");
		fflush(stdout);

		// for(k=0, i = 0; k<s->n_spe; k++)
		// {
		// 	N += s->natom_spe[k];
		// 		for(; i < N; i++)
		// 		{
		// 			atom= s->a_set->atom_array[i];

		// 			for (j = 0 ; j < 3; j++){
		// 				s0=ReaderGetEntry(fr, 8+i,j);
		// 				tmp[j]=atof(s0);
		// 			}

		// 			atom->spe = k;
		// 			atom->tau[0]=tmp[0];
		// 			atom->tau[1]=tmp[1];
		// 			atom->tau[2]=tmp[2];
		// 		}
		// }

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
		Aset->atom_array[i]->n_neighbor = Aset->neigh_storage+i*Nspe;
	}
	return Aset;
}

static inline void BuildRe( Lattice *s)
{
	int i,j;
	double vol;
	double a[3][3];
	double b[3][3];

	for(i=0; i<3; i++)
	{
		a[i][0]=s->a[i][0];
		a[i][1]=s->a[i][1];
		a[i][2]=s->a[i][2];
	}
	
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
			s->b[i][j]=b[i][j]*TWOPI/vol;  
		}
	}
}


static inline void ReadPoscar	(Lattice *s, char *title )
{
	int			N;
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


static inline void GVecInitial(Lattice *s, int N)
{
	int i,j,k;
	double *G_bulk;
  	s->G_vec =  SafeCalloc(N*N*N, sizeof(double *));
	G_bulk   =  SafeCalloc(3*N*N*N, sizeof(double) );
	for(i=0; i<N*N*N; i++)
	{
		s->G_vec[i] = G_bulk+i*3;
	}

}

Lattice *LatticeInitial ( char *filename, int K_max )
{
	Lattice *L;
	L = SafeCalloc(1,sizeof(Lattice));
	ReadPoscar(L, filename);
	BuildRe(L);
	GVecInitial(L,K_max);
	return L;


}