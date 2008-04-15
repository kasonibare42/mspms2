#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <ctype.h>
#include "mspms2.h"

/**
 * \brief MD moves for HMC etc.
 */
int fnMDmove(int nStepMD, void(*pfnRezero)(), int (*pfnAlgorithm)())
{
	int ii;

	// rezero thermostat
	if (pfnRezero)
	{
		pfnRezero();
	}
	
	// initialize energy difference
	fDeltaU = 0.0;

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
	ushift_old = ushift;

	// MD moves
	// callback
	for (ii=0; ii<nStepMD; ii++)
	{
		pfnAlgorithm();
	}

	// calculate the energy difference
	fDeltaU = (uinter - uinter_old) + (uintra - uintra_old) + (ukin - ukin_old) + (ushift - ushift_old);  // kinetic energy should be included??

	// Hamotonial difference
	dH = fDeltaU*R_RGAS/treq;

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
		// printf("delta_U=%lf   dH=%lf   %lf < %lf ?\n",fDeltaU, dH, rndnum[0], exp(-dH));
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
		ushift = ushift_old;
	}

	return isAccept;
}
