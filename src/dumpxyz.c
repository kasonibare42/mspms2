/**
 * Project: mspms2
 * File: dumpxyz.c
 * 
 * Copyright (C) 2008    Yang Wang <ywangd@gmail.com>
 * Created @ 15/04/2008
 * 
 * Description:
 * 
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <ctype.h>
#include "mspms2.h"

int snapshot()
{
	int ii;

	fpss = fopen(SNAPSHOT,"w");
	fprintf(fpss, "%d\n", natom);
	fprintf(fpss, "%d %lf %lf %lf\n", istep, boxlx, boxly, boxlz);
	for (ii=0; ii<natom; ii++)
	{
		fprintf(fpss, "%s  %lf  %lf  %lf\n", atomname[ii], xx[ii], yy[ii],
				zz[ii]);
	}
	fclose(fpss);

	return 0;
}

int trajectory()
{
	int ii, jj;

	// write out molecular information at the beginning of the trajectory file
	if (!nframe)
	{
		fwrite(&nspecie, sizeof(int), 1, fptrj);
		fwrite(&nmole, sizeof(int), 1, fptrj);
		fwrite(&natom, sizeof(int), 1, fptrj);
		for (ii=0; ii<nspecie; ii++)
		{
			fwrite(&nmole_per_specie[ii], sizeof(int), 1, fptrj);
			fwrite(&sample_natom_per_mole[ii], sizeof(int), 1, fptrj);
			for (jj=sample_mole_first_atom_idx[ii]; jj
					<sample_mole_last_atom_idx[ii]; jj++)
			{
				fwrite(sample_atomname[jj], sizeof(char), 5, fptrj);
				fwrite(&sample_aw[jj], sizeof(double), 1, fptrj);
			}
		}
	}
	nframe++;
	fwrite(nmole_per_specie, sizeof(int), nspecie, fptrj);
	fwrite(xx, sizeof(double), natom, fptrj);
	fwrite(yy, sizeof(double), natom, fptrj);
	fwrite(zz, sizeof(double), natom, fptrj);

	return 0;
}

