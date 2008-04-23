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
	int ii, iSpecie, iAtom;

	fpss = fopen(SNAPSHOT,"w");
	fprintf(fpss, "%d\n", natom);
	fprintf(fpss, "%d %lf %lf %lf\n", istep, boxlx, boxly, boxlz);
	for (ii=0; ii<natom; ii++)
	{
		get_specie_and_relative_atom_id(ii, &iSpecie, &iAtom);
		fprintf(fpss, "%s  %lf  %lf  %lf\n", sample_mole[iSpecie].atom_name[iAtom], xx[ii], yy[ii],
				zz[ii]);
	}
	fclose(fpss);

	return 0;
}

int trajectory()
{
	int ii, jj;

	// write out molecular information at the beginning of the trajectory file
	if (nframe==0)
	{
		fwrite(&nspecie, sizeof(int), 1, fptrj);
		fwrite(&nmole, sizeof(int), 1, fptrj);
		fwrite(&natom, sizeof(int), 1, fptrj);
		for (ii=0; ii<nspecie; ii++)
		{
			fwrite(&nmole_per_specie[ii], sizeof(int), 1, fptrj);
			fwrite(&sample_mole[ii].natom, sizeof(int), 1, fptrj);
			for (jj=0; jj<sample_mole[ii].natom; jj++)
			{
				fwrite(sample_mole[ii].atom_name[jj], sizeof(char), 5, fptrj);
				fwrite(&sample_mole[ii].aw[jj], sizeof(double), 1, fptrj);
			} // Sample atom loop
		} // Specie loop
	}
	nframe++;
	fwrite(nmole_per_specie, sizeof(int), nspecie, fptrj);
	fwrite(xx, sizeof(double), natom, fptrj);
	fwrite(yy, sizeof(double), natom, fptrj);
	fwrite(zz, sizeof(double), natom, fptrj);

	return 0;
}

