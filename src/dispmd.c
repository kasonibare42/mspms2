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
bool md_move(int nStepMD)
{
	int ii;
	double du;

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
	for (ii=0; ii<nStepMD; ii++)
	{
		vver();
	}

	// calculate the energy difference
	du = (uinter - uinter_old) + (uintra - uintra_old) + (ukin - ukin_old) + (ushift - ushift_old);  // kinetic energy should be included??

	// Hamotonial difference
	dH = du/treq;

	counts[20]++; // canonical moves

	// Check if the move is accepted
	isAccept = false;
	if (dH <= 0.0)
	{
		isAccept = true;
	}
	else
	{
		ranmar(rndnum, 1); 
		if (rndnum[0] < exp(-dH))
		{
			isAccept = true;
		}
	}
	if (isAccept == true)
	{
		counts[21]++; // accepted canonical moves
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
