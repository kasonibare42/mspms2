#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#define NATOM_MAX	4000
#define NMOLE_MAX	1000
#define NSPECIE_MAX	3
#define SAMPLE_NATOM_MAX	200

#define _safealloc(pt,num,size)         pt=calloc(num,size); assert(pt!=NULL)
#define _safefree(pt)	if ((pt)!=NULL) {free(pt); (pt)=NULL;}


int usage()
{}

int main(int argc, char *argv[])
{
	int ii, jj, kk, ff;
	char *trjfile = NULL; // trajectory file
	char *outputfile = NULL; // output msd file
	int nframe;
	int iFrame, iAtom, iSampleAtom;
	int nspecie, nmole, natom;
	int nmole_per_specie[NSPECIE_MAX];
	int sample_natom_per_mole[NSPECIE_MAX];
	char sample_atomname[SAMPLE_NATOM_MAX][5];
	double sample_aw[SAMPLE_NATOM_MAX];
	double xx[NATOM_MAX], yy[NATOM_MAX], zz[NATOM_MAX];

	FILE *fptrj, *fpouts;

	for (ii=1;ii<argc;ii++)
	{
		if (!strcmp(argv[ii],"-in"))
		{
			trjfile = argv[++ii];
			printf("%s\n",trjfile);
		}
		else if (!strcmp(argv[ii],"-out"))
		{
			outputfile = argv[++ii];
		}
		else if (!strcmp(argv[ii],"-natom"))
		{
			natom = atoi(argv[++ii]);
		}
		else if (!strcmp(argv[ii],"-nframe"))
		{
			nframe = atoi(argv[++ii]);
		}
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
	if (!trjfile)
	{
		_safealloc(trjfile, 20, sizeof(char));
		strcpy(trjfile, "trj.mspms");
	}
	if (!outputfile)
	{
		_safealloc(outputfile, 20, sizeof(char));
		strcpy(outputfile, "trjmspms.xyz");
	}

	fptrj = fopen(trjfile, "rb");
	fpouts = fopen(outputfile,"w");
	fread(&nspecie, sizeof(int), 1, fptrj);
	fread(&nmole, sizeof(int), 1, fptrj);
	fread(&natom, sizeof(int), 1, fptrj);
	assert(natom<NATOM_MAX);
	iAtom = 0;
	for (ii=0;ii<nspecie;ii++)
	{
		fread(&nmole_per_specie[ii], sizeof(int), 1, fptrj);
		fread(&sample_natom_per_mole[ii], sizeof(int), 1, fptrj);
		for (jj=iAtom;jj<iAtom+sample_natom_per_mole[ii];jj++)
		{
			fread(sample_atomname[jj], sizeof(char), 5, fptrj);
			fread(&sample_aw[jj], sizeof(double), 1, fptrj);
		}
		iAtom += sample_natom_per_mole[ii];
	}

	iFrame = 0;
	while (!feof(fptrj) )
	{
		fread(nmole_per_specie, sizeof(int), nspecie, fptrj);
		natom = 0; 
		for (ii=0;ii<nspecie;ii++)
		{
			natom += nmole_per_specie[ii]*sample_natom_per_mole[ii];
		}
		assert(natom<NATOM_MAX);
		fread(xx,sizeof(double),natom,fptrj);
		fread(yy,sizeof(double),natom,fptrj);
		fread(zz,sizeof(double),natom,fptrj);

		fprintf(fpouts,"%d\n%d\n",natom,iFrame);
		iAtom = 0;
		iSampleAtom = 0;
		for (ii=0;ii<nspecie;ii++)
		{
			for (jj=0;jj<nmole_per_specie[ii];jj++)
			{
				// printf("sample_natom_per_mole[%d]=%d\n",ii,sample_natom_per_mole[ii]);
				for (kk=iSampleAtom; kk<iSampleAtom+sample_natom_per_mole[ii]; kk++)
				{
					// printf("ii=%d  jj=%d   kk=%d   iAtom=%d\n",ii,jj,kk,iAtom);
					fprintf(fpouts, "%s  %lf   %lf   %lf\n", sample_atomname[kk], xx[iAtom], yy[iAtom], zz[iAtom]);
					iAtom++;
				} 
			}
			iSampleAtom += sample_natom_per_mole[ii];
		}

		// printf("iFrame=%d\n", iFrame);
		iFrame++;
	}
	fclose(fptrj);
	fclose(fpouts);

}





