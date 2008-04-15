/**
 * Project: mspms2
 * File: loadsave.c
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

int saveit()
{
	fpsave = fopen(SAVE,"wb");

	fwrite(xx, sizeof(double), natom, fpsave);
	fwrite(yy, sizeof(double), natom, fpsave);
	fwrite(zz, sizeof(double), natom, fpsave);
	fwrite(vx, sizeof(double), natom, fpsave);
	fwrite(vy, sizeof(double), natom, fpsave);
	fwrite(vz, sizeof(double), natom, fpsave);

	fwrite(&qq, sizeof(double), 1, fpsave);
	fwrite(&ps, sizeof(double), 1, fpsave);
	fwrite(&gg, sizeof(double), 1, fpsave);
	fwrite(&ss, sizeof(double), 1, fpsave);
	fwrite(&qqs, sizeof(double), 1, fpsave);
	fwrite(&pss, sizeof(double), 1, fpsave);
	fwrite(&ggs, sizeof(double), 1, fpsave);
	fwrite(&sss, sizeof(double), 1, fpsave);

	fwrite(&vts, sizeof(double), 1, fpsave);
	fwrite(&rts, sizeof(double), 1, fpsave);
	fwrite(&vbs, sizeof(double), 1, fpsave);

	fwrite(&boxlx, sizeof(double), 1, fpsave);
	fwrite(&boxly, sizeof(double), 1, fpsave);
	fwrite(&boxlz, sizeof(double), 1, fpsave);
	fwrite(&boxv, sizeof(double), 1, fpsave);

	fwrite(&bEquilibrium, sizeof(bool), 1, fpsave);
	
	fwrite(&istep, sizeof(int), 1, fpsave);
	fwrite(counts, sizeof(int), NCOUNTS_MAX, fpsave);
	fwrite(accum, sizeof(double), NCOUNTS_MAX, fpsave);

	fclose(fpsave);

	return 0;
}

int loadit()
{
	int ii;

	fprintf(stderr,"loading from saved file...\n");
	fprintf(fpouts, "loading from saved file...\n");

	fpload = fopen(LOAD,"rb");

	fread(xx, sizeof(double), natom, fpload);
	fread(yy, sizeof(double), natom, fpload);
	fread(zz, sizeof(double), natom, fpload);
	fread(vx, sizeof(double), natom, fpload);
	fread(vy, sizeof(double), natom, fpload);
	fread(vz, sizeof(double), natom, fpload);

	// recaculate kinetic energy and temperature using the loaded velocities
	ukin = 0.0;
	for (ii=0; ii<natom; ii++)
	{
		ukin += aw[ii]*(vx[ii]*vx[ii]+vy[ii]*vy[ii]+vz[ii]*vz[ii]);
	}
	ukin = 0.5*ukin;
	tinst = 2.0*ukin/(RGAS*nfree);

	// calculate instant temperature
	fread(&qq, sizeof(double), 1, fpload);
	fread(&ps, sizeof(double), 1, fpload);
	fread(&gg, sizeof(double), 1, fpload);
	fread(&ss, sizeof(double), 1, fpload);
	fread(&qqs, sizeof(double), 1, fpload);
	fread(&pss, sizeof(double), 1, fpload);
	fread(&ggs, sizeof(double), 1, fpload);
	fread(&sss, sizeof(double), 1, fpload);

	fread(&vts, sizeof(double), 1, fpload);
	fread(&rts, sizeof(double), 1, fpload);
	fread(&vbs, sizeof(double), 1, fpload);

	fread(&boxlx, sizeof(double), 1, fpload);
	fread(&boxly, sizeof(double), 1, fpload);
	fread(&boxlz, sizeof(double), 1, fpload);
	fread(&boxv, sizeof(double), 1, fpload);

	// recalculate the thermostat energy
	utsbs = 0.5*Qbs*vbs*vbs + 0.5*Qts*vts*vts + (nfree+1)*RGAS*treq*rts + preq
			*boxv*PA_A3_TO_J_PER_MOL;
	if (isLJlrcOn)
	{
		// recalculate the long range corrections since the box size may be changed
		calculate_ljlrc();
	}

	// read counters and accumulators only if its a continue run
	if (fStart_option==CONTINUE)
	{
		fread(&bEquilibrium, sizeof(bool), 1, fpload);
		fread(&istep, sizeof(int), 1, fpload);
		nstep_start = istep + 1;
		fread(counts, sizeof(int), NCOUNTS_MAX, fpload);
		fread(accum, sizeof(double), NCOUNTS_MAX, fpload);
	}

	fclose(fpload);

	return 0;
}
