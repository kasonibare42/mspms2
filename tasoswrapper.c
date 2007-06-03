#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "vars.h"

extern "C" void initpotentialgrid_(int*, double*, double*, double*, int*, int*, int*, double*);
extern "C" void pass_grid_file_name_(char*, int*);
extern "C" void read_grids_(int* nspecies_yang);

int init_tasos_grid()
{
    int	ii, jj;
    const int datalen = 200;
    char buffer[200];
    double auc, buc, cuc, nanotuberadius;
    int	na, nb, nc;
    int	nspecies_yang;
    int tnuatoms;
    char szgrid[12][80]; // assuming only 12 types of grids
    FILE *fptasos;

    fptasos = fopen("tasos.mspms","r");

    sscanf(fgets(buffer,datalen,fptasos), "%lf %lf %lf %lf", &auc, &buc, &cuc, &nanotuberadius);
    sscanf(fgets(buffer,datalen,fptasos), "%d %d %d", &na, &nb, &nc);
    sscanf(fgets(buffer,datalen,fptasos), "%d", &nspecies_yang);
    sscanf(fgets(buffer,datalen,fptasos), "%d", &tnuatoms);
    for (ii=0;ii<tnuatoms;ii++)
       	sscanf(fgets(buffer,datalen,fptasos), "%s", szgrid+ii);

    initpotentialgrid_(&tnuatoms, &auc, &buc, &cuc, &na, &nb, &nc, &nanotuberadius);
    fprintf(stderr,"init potential ok\n");
    for (ii=0;ii<tnuatoms;ii++)
    {
	jj = ii + 1;
       	pass_grid_file_name_(szgrid[ii], &jj);
    }
    fprintf(stderr,"pass grid file name ok\n");
    read_grids_(&nspecies_yang);
    fprintf(stderr,"read grids ok\n");

    fclose(fptasos);

}
