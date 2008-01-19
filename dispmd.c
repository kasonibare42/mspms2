#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <ctype.h>
#include "random.h"
#include "vars.h"
#include "funcs.h"

/**
 * \brief MD moves for HMC etc.
 */
int fnMDmove(int nStepMD, int (*pfnAlgorithm)())
{
	int ii;
	
	// initialize energy difference
	fDeltaU = 0.0;

	// velocity init and energy calcualtions have been done for the 1st step outside the loop
	if (istep!=nstep_start)
	{
		velinit();
		erfrc();
		rafrc(); // we may not need rafrc here??
	}

	// save the old states
	for (ii=0; ii<natom; ii++)
	{
		xx_old[ii] = xx[ii];
		yy_old[ii] = yy[ii];
		zz_old[ii] = zz[ii];
	}
	ukin_old = ukin;
	tinst_old = tinst;
	uinter_old = uinter;
	uintra_old = uintra;
	uvdw_old = uvdw;
	ubond_old = ubond;
	uangle_old = uangle;
	udih_old = udih;
	uimp_old = uimp;
	uewald_old = uewald;
	uwolf_old = uwolf;
	ucoulomb_old = ucoulomb;
	usflj_old = usflj;
	virial_inter_old = virial_inter;
	virial_intra_old = virial_intra;

	// MD moves
	// callback
	for (ii=0; ii<nStepMD; ii++)
	{
		pfnAlgorithm();
	}

	// calculate the energy difference
	fDeltaU = (uinter - uinter_old) + (uintra - uintra_old)
			+ (ukin - ukin_old);

	// Hamotonial difference
	dH = fDeltaU*rRgas/treq;

	icounter[20]++; // canonical moves

	// check if the move is accepted
	isAccept = 0;
	if (dH <= 0.0)
	{
		isAccept = 1;
	}
	else
	{
		ranmar(rndnum, 1);
		if (rndnum[0] < exp(-dH))
		{
			isAccept = 1;
		}
	}
	if (isAccept == 1)
	{
		icounter[21]++; // accepted canonical moves
	}
	else
	{
		// restore old values
		for (ii=0; ii<natom; ii++)
		{
			xx[ii] = xx_old[ii];
			yy[ii] = yy_old[ii];
			zz[ii] = zz_old[ii];
		}
		ukin = ukin_old;
		tinst = tinst_old;
		uinter = uinter_old;
		uintra = uintra_old;
		uvdw = uvdw_old;
		ubond = ubond_old;
		uangle = uangle_old;
		udih = udih_old;
		uimp = uimp_old;
		uewald = uewald_old;
		uwolf = uwolf_old;
		ucoulomb = ucoulomb_old;
		usflj = usflj_old;
		virial_inter = virial_inter_old;
		virial_intra = virial_intra_old;
	}

	return 0;
}
