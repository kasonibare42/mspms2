#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#define natom_max 9999

int usage()
{}

int main(int argc, char *argv[])
{
	int ii, jj, kk, ff;
	char *trjfile; // trajectory file
	char *outputfile; // output msd file
	int nframe;

	int natom;
	double xx[natom_max], yy[natom_max], zz[natom_max];



	FILE *fptrj, *fpouts;

	for (ii=1;ii<argc;ii++)
	{
		if (!strcmp(argv[ii],"-in"))
		{
			trjfile = argv[++ii];
			printf("%s\n",trjfile);
		}
		else if (!strcmp(argv[ii],"-out"))
			outputfile = argv[++ii];
		else if (!strcmp(argv[ii],"-natom"))
			natom = atoi(argv[++ii]);
		else if (!strcmp(argv[ii],"-nframe"))
			nframe = atoi(argv[++ii]);
		else if (!strcmp(argv[ii],"-h"))
		{
			usage();
			exit(0);
		}
		else
		{
			fprintf(stderr,"Error: Command-line argument '%s' not recognized.\n", argv[ii]);
			exit(1);
		}
	}

	assert(natom<natom_max);

	fptrj = fopen(trjfile, "rb");
	fpouts = fopen(outputfile,"w");
	for (ii=0;ii<nframe;ii++)
	{
		fread(xx,sizeof(double),natom,fptrj);
		fread(yy,sizeof(double),natom,fptrj);
		fread(zz,sizeof(double),natom,fptrj);

		fprintf(fpouts,"%d\n%d\n",natom,ii);
		for (jj=0;jj<natom;jj++) 
		{
			fprintf(fpouts,"Cu  %lf   %lf   %lf\n",xx[jj],yy[jj],zz[jj]);

		}
	}
	fclose(fptrj);
	fclose(fpouts);

}





