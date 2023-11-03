#include <stdio.h>  
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "epm.h"
#include "memory.h"
#include "util.h"
#include "reader.h"

static void ReadKpathName(Reader *fr,int i0,EPM *s);
static void ReadSimulation(Reader *fr,int i0,EPM *s);
static void ReadNThreads(Reader *fr,int i0,EPM *s);
static void ReadSimType(Reader *fr,int i0,EPM *s);
static void ReadType(Reader *fr,int i0,EPM *s);
static void ReadSimwave(Reader *fr,int i0,EPM *s);
static void ReadEcut(Reader *fr,int i0,EPM *s);

static void ReadKpath(EPM *s);
static void ReadKlist(EPM *s);
static void ReadK_sym( Reader *fr,EPM *s);

/*********
 * Read the input file
 * 1st line: [simulation name] //default value is current directory 
 * 2nd line: [simulation type] //default value is BAND
 * 3rd line: [potential type]  //default value is LOC
 * 4th line: [k_path file name]//default value is k_path.in
 * 5th line: [E_cut]           //in unit of Rydberg
 * 6th line: [n_threads]       //number of threads
 * 7th line: [simulation wave]       //save wave function if 1, don't if 0
 * ********/
void EPMReadInput(EPM *s,const char *filename)
{
    Reader *fr;
    int i0;
    char *s0;
    //

    if(!FileExists(filename))
    {
        free(filename);
        printf("No input file found.\n ");
        fflush(stdout);
        return;
    }
    fr = ReaderReadFile(filename);
    for(i0 = 0;i0 < ReaderGetFileLength(fr); i0++)
    {
        s0 = ReaderGetEntry(fr,i0,0);
        if(!strcmp(s0,SIMNAME)) ReadSimulation(fr,i0,s);
        else if(!strcmp(s0,SIMWAVE)) ReadSimwave(fr,i0,s);
        else if(!strcmp(s0,SIMTYPE)) ReadSimType(fr,i0,s);
        else if(!strcmp(s0,NTHREADS)) ReadNThreads(fr,i0,s);
        else if(!strcmp(s0,POTTYPE)) ReadType(fr,i0,s);
        else if(!strcmp(s0,KPATH)) ReadKpathName(fr,i0,s);
        else if(!strcmp(s0,ECUT)) ReadEcut(fr,i0,s);
    }
    // Read the statistics
    if(!strcmp(s->simtype,"BAND"))
        ReadKpath(s);
    else if(!strcmp(s->simtype,"KLIST"))
        ReadKlist(s);
    ReaderFree(fr);
}




static void ReadSimulation(Reader *fr,int i0,EPM *s)
{
    char *s0;
    //
    s0 = ReaderGetEntry(fr,i0,1);
    StringClone(s->simname,"./");//default value is current directory
    StringClone(s->simname,s0);//s->simname = s0;
    printf("Simulation name: %s \n",s->simname);
    fflush(stdout);
}

static void ReadNThreads(Reader *fr,int i0,EPM *s)
{
    char *s0;
    //
    s0 = ReaderGetEntry(fr,i0,1);
    s->n_threads = atoi(s0);
}


static void ReadSimType(Reader *fr,int i0,EPM *s)
{
    char *s0;
    //
    s0 = ReaderGetEntry(fr,i0,1);
    StringClone(s->simtype,"BAND");//default value is BAND
    StringClone(s->simtype,s0);//SimType = s0;
    printf("Simulation type: %s \n",s->simtype);
    fflush(stdout);
}

static void ReadSimwave(Reader *fr,int i0,EPM *s)
{
    char *s0;
    //
    s0 = ReaderGetEntry(fr,i0,1);
    s->simwave = atoi(s0);
    if(s->simwave)
        printf("Eigen energy and Wave function will be saved.\n");
    else
        printf("Only Eigen energy will be saved.\n");
    fflush(stdout);

}

static void ReadEcut(Reader *fr,int i0,EPM *s)
{
    char *s0;
    //
    s0 = ReaderGetEntry(fr,i0,1);
    s->Emax = atof(s0)*RYD; //s->Emax = convert s0 to Rydberg
    printf("Cut off energy: %f eV\n",s->Emax);
    fflush(stdout);
}

static void ReadType(Reader *fr,int i0,EPM *s)
{
    char *s0;
    //
    s0 = ReaderGetEntry(fr,i0,1);
    s->pottype=1;//default value is LOC
    if(!strcmp(s0,"SO"))
        s->pottype=2;
    printf("Potential  type: %s \n",s0);
    fflush(stdout);
}

static void ReadKpathName(Reader *fr,int i0,EPM *s)
{
    char *s0;
    char        *name;
	char        *path;
    //
    s0 = ReaderGetEntry(fr,i0,1);
    StringClone(path, s->simname);
	StrCat(&path, "/");

    StringClone(name,"k_path.in");//default value is current directory
    StringClone(name,s0);

    s->kpath_name = FullPath(path, name);//s->kpath = "simname/s0";
}

static void ReadK_sym( Reader *fr,EPM *s)
{
    int i;
    char *s0;
    for(i=0;i<=s->n_kpath;i++)
    {
        s0 = ReaderGetEntry(fr,i+2,0);
        s->k_symetry[i*3+0] = atof(s0);
        s0 = ReaderGetEntry(fr,i+2,1);
        s->k_symetry[i*3+1] = atof(s0);
        s0 = ReaderGetEntry(fr,i+2,2);
        s->k_symetry[i*3+2] = atof(s0);
    }
}

/******
 * Read the k_path file
 * 
 * 1st line: [k_path info]          //should be comments on the k path
 * 2nd line: [n_kpath] [n_k_p_path] //number of k path and number of k points per path
 * 3rd line: [kx] [ky] [kz]         //the first k point 
 * 4th line: [kx] [ky] [kz]         //the second k point
 * .
 * .
 * .
 *****/
static void ReadKpath(EPM *s)
{
    char *s0;
    char *filename=s->kpath_name;
    int linenum=0;
    Reader *fr;

    printf("K_path info:  ");
    
    fr = ReaderReadFile(filename);

    for (int i0=1; i0<ReaderGetLineLength(fr,linenum); i0++)
    {
        printf(" %s", ReaderGetEntry(fr,linenum,i0));
            
    }
    printf("\n");
    fflush(stdout);

	s0 = ReaderGetEntry(fr,1,0);
	s->n_kpath = atof(s0);
    s0 = ReaderGetEntry(fr,1,1);
    s->n_k_p_path = atof(s0);
    s->k_symetry = SafeCalloc(3*s->n_kpath+3,sizeof(double));
	
	ReadK_sym(fr, s);
    ReaderFree(fr);
}





static void ReadK_point( Reader *fr,EPM *s)
{
    int i;
    char *s0;
    for(i=0;i<s->n_kpoint;i++)
    {
        s0 = ReaderGetEntry(fr,i+2,0);
        s->k_list[i*3+0] = atof(s0);
        s0 = ReaderGetEntry(fr,i+2,1);
        s->k_list[i*3+1] = atof(s0);
        s0 = ReaderGetEntry(fr,i+2,2);
        s->k_list[i*3+2] = atof(s0);
    }
}
/******
 * Read the k_path file
 * 
 * 1st line: [k_path info]          //should be comments on the k path
 * 2nd line: [n_kpoint]             //number of k points
 * 3rd line: [kx] [ky] [kz]         //the first k point
 * 4th line: [kx] [ky] [kz]         //the second k point
 * .
 * .
 * .
 *****/
static void ReadKlist(EPM *s)
{
    char *s0;
    char *filename=s->kpath_name;
    int linenum=0;

    Reader *fr;
    fr = ReaderReadFile(filename);
	s0 = ReaderGetEntry(fr,0,0);

    printf("K_list info:  ");
    for (int i0=1; i0<ReaderGetLineLength(fr,linenum); i0++)
    {
        printf(" %s", ReaderGetEntry(fr,linenum,i0));
    }
    printf("\n");     
    fflush(stdout);

	s0 = ReaderGetEntry(fr,1,0);
	s->n_kpoint = atof(s0);
    
    s->k_list = SafeCalloc(3*s->n_kpoint,sizeof(double));
	
	ReadK_point(fr, s);
    ReaderFree(fr);
}