/**
 * Project: mspms2
 * File: bndfrc.c
 * 
 * Copyright (C) 2008    Yang Wang <ywangd@gmail.com>
 * Created @ Apr 22, 2008
 * 
 * Description:
 * 
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include "mspms2.h"

int bndfrc(int iSpecie, int iBond, int iabs, double *uij, double *virial_ij)
{
	int mm1, mm2, ii1, ii2;
	double rxij, ryij, rzij, fij, fxij, fyij, fzij;
	double rijsq, rij, delta_r, De, exp_term, one_minus_exp_term, rmax2;
	PSAMPLE_MOLECULE pSampleMole;

	pSampleMole = sample_mole + iSpecie;
	
	mm1 = pSampleMole->bnd_idx[iBond][0]; // relative atom id 1
	mm2 = pSampleMole->bnd_idx[iBond][1]; // relative atom id 2
	// Calculate absolute atom id
	ii1 = iabs + mm1; // absolute atom id 1
	ii2 = iabs + mm2; // absolute atom id 2
	// Calculate separation
	rxij = xx[ii1] - xx[ii2];
	ryij = yy[ii1] - yy[ii2];
	rzij = zz[ii1] - zz[ii2];
	rijsq = rxij*rxij + ryij*ryij + rzij*rzij;
	rij = sqrt(rijsq);

	switch (pSampleMole->bnd_type[iBond])
	{
	case BOND_NONE: // 0
		break;
	case BOND_HARMONIC: // 1
		delta_r = rij - pSampleMole->Req[iBond];
		*uij = 0.5*pSampleMole->Kb[iBond]*delta_r*delta_r;
		fij = -pSampleMole->Kb[iBond]*delta_r/rij;
		break;
	case BOND_MORSE: // 2
		delta_r = rij - pSampleMole->Req[iBond];
		De = pSampleMole->Kb[iBond];
		exp_term = exp(-pSampleMole->alpha[iBond]*delta_r);
		one_minus_exp_term = 1.0 - exp_term;
		*uij = De*one_minus_exp_term*one_minus_exp_term;
		fij = -2.0*De*pSampleMole->alpha[iBond]*one_minus_exp_term*exp_term/rij;
		break;
	case BOND_FENE: // 3
		rmax2 = pSampleMole->Req[iBond]*pSampleMole->Req[iBond];
		*uij = pSampleMole->Kb[iBond]*rmax2*log(1.0-rijsq/rmax2);
		fij = -2.0*pSampleMole->Kb[iBond]*rmax2/(rmax2-rijsq);
		break;
	default:
		printf("unknown bond type.\n");
		exit(1);
	} // switch for different bond type
	*virial_ij = fij*rijsq;
	fxij = fij*rxij;
	fyij = fij*ryij;
	fzij = fij*rzij;
	// force on atom 1
	fxs[ii1] += fxij;
	fys[ii1] += fyij;
	fzs[ii1] += fzij;
	// force on atom 2
	fxs[ii2] -= fxij;
	fys[ii2] -= fyij;
	fzs[ii2] -= fzij;
	
	return 0;
}


