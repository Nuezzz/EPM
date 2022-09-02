#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "memc.h"
#include "error.h"
#include "memory.h"
#include "util.h"
#include "otmesher.h"
#include "scattering.h"
#include "memc.h"
#include "demc.h"
#include "particle.h"
#include "version.h"
//
//
double timer_interval(struct timespec start, struct timespec end);
void DFISEConvert(const char *basename);
//
//
void GMSHConvert(const char *basename);
//
//

//global machine precision epsilon values

static void Help()
{
    printf("To use the MC3D, please enter the following into the command line:\n");
    printf(">\tMC3D_Executable Tool Parameters\n");
    printf("Available Tools: HELP, GMSH, DFISE, TPMESHER, OTMESHER, SCATTERING, MEMC, DEMC, TEST\n");
    printf("\n");
    printf("Examples of commonly used tools\n");
    printf(">\tDEMC Example: ./MC3D_debug DEMC demc.in\n");
    printf(">\tMEMC Example: ./MC3D_debug MEMC memc.in\n");
    printf(">\tGMSH Example: ./MC3D_debug GMSH meshname\n");
    printf("\n");
    printf("For more information on a specific tool, please enter the following into the command line (under construction):\n");
    printf(">\tMC3D_Executable Tool\n");
}
//
//
int main(int argc, char **argv)
{
    int flag = 0;
	printf("MC3D Software Suite ");
    VersionPrint(stdout);
    printf("\n");
    fflush(stdout);
    ErrorStreamOpen("error.log");
    if(argc > 1)
    {
        if(argc == 3)
        {
            if (!strcmp(argv[1], "DFISE"))
            {
                printf("Tool: DFISE File Converter\n");
                fflush(stdout);
                DFISEConvert(argv[2]);
                flag = 1;
            }
            else if (!strcmp(argv[1], "GMSH"))
            {
                flag = 1;
                printf("Tool: GMSH Converter\n");
                fflush(stdout);
                GMSHConvert(argv[2]);
            }
            else if (!strcmp(argv[1], "TPMESHER"))
            {
                flag = 1;
                printf("Tool: Tensor Product Mesher (TPMESHER)\n");
                fflush(stdout);
                TPMesher(argv[2]);
            }
            else if (!strcmp(argv[1], "OTMESHER"))
            {
                flag = 1;
                printf("Tool: Octree Tetrahedral Mesher (OTMESHER)\n");
                fflush(stdout);
                OTMGenerate(argv[2]);
            }
            else if (!strcmp(argv[1], "SCATTERING"))
            {
                flag = 1;
                printf("Tool: Full-Band Scattering Rate Calculator (SCATTERING)\n");
                fflush(stdout);
                ScatteringRun(argv[2]);
            }
            else if (!strcmp(argv[1], "MEMC"))
            {
                flag = 1;
                printf("Tool: Material Ensemble Monte Carlo (MEMC)\n");
                fflush(stdout);
                memc(argv[2]);
            }
            else if (!strcmp(argv[1], "DEMC"))
            {
                flag = 1;
                printf("Tool: Device Ensemble Monte Carlo (DEMC)\n");
                fflush(stdout);
                demc(argv[2]);
            }
            else
            {
                printf("Tool keyword not recognized!\n");
                fflush(stdout);             
            }
        }
        if(argc == 2)
        {
            if (!strcmp(argv[1], "GMSH"))
            {
                flag = 1;
                printf("Help with GMSH Converter\n");
                printf("To use this tool, please enter the mesh name after MC3D_Executable GMSH in the command line\n");
                printf("For example, if the mesh name is 'mesh3d', with a symbolically linked MC3D_debug executable,\n");
                printf("the command << ./MC3D_debug GMSH mesh3d >> will require the files ./mesh3d/mesh3d.geo and ./mesh3d/mesh3d.in\n");
                fflush(stdout);
            }
            else if (!strcmp(argv[1], "HELP"))
            {
                printf("Tool: Not-So-Useful Helper\n");
                fflush(stdout);
                Help();
                flag = 1;
            }
            else if (!strcmp(argv[1], "TEST"))
            {
                flag = 1;
                printf("Tool: A Vanilla Testing Environment\n");
                fflush(stdout);
                printf("End of test run\n");
            }
            else
            {
                printf("Tool keyword not recognized or not yet implemented!\n");
                fflush(stdout);             
            }
        }
    }
    else
    {
        printf("No Tools were selected!\n");
    }
    if (!flag)
    {
        Help();
    }
    ErrorStreamClose();
    return 0;
}
