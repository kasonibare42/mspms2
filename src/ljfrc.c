/**
 * Project: mspms2
 * File: ljfrc.c
 * 
 * Copyright (C) 2008    Yang Wang <ywangd@gmail.com>
 * Created @ 21/04/2008
 * 
 * Description:
 * 
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include "mspms2.h"

/**
 * Energy and force calculations for LJ interactions.
 * Output energy, raw force and shift energy
 */
int ljfrc(double rijsq, double sigmaij, double epsilonij, double *uij,
		double *fij, double *uijshift)
{
	double LJswitch;
	double r_rijsq, r6, r12;
	double sigmaij2, sigmaij6;

	// if switch potential for LJ is on, then calculate the switch
	if (isLJswitchOn)
	{
		if (rijsq<rcutonsq)
		{
			LJswitch = 1.0;
		}
		else
		{
			LJswitch = (rcutoffsq-rijsq)*(rcutoffsq-rijsq)*(rcutoffsq+2.0
					*rijsq-3.0*rcutonsq)/roff2_minus_ron2_cube;
		}
	}
	r_rijsq = 1.0/rijsq;
	sigmaij2 = sigmaij*sigmaij;
	r6 = sigmaij2*r_rijsq;
	r6 = r6*r6*r6;
	r12 = r6*r6;
	*uij = 4.0*epsilonij*(r12-r6);
	if (isLJswitchOn) // if switch is used
	{
		if (rijsq<rcutonsq)
		{
			*fij = 24.0*epsilonij*(2*r12-r6)*r_rijsq;
			// we do not need to modify energy because LJswitch is 1.0 here
		}
		else
		{
			*fij = 24.0*epsilonij*(2*r12-r6)*r_rijsq*LJswitch -*uij*12.0
					*(rcutoffsq-rijsq)*(rcutonsq-rijsq)/roff2_minus_ron2_cube;
			*uij *= LJswitch; // switch energy modifier
		}
		*uijshift = 0.0; // If switch is used, we do not need shift energy.
	}
	else
	{
		*fij = 24.0*epsilonij*(2*r12-r6)*r_rijsq;
		// Only when LJ swith is not used, we calculate shift energies shift energies.
		sigmaij6 = sigmaij2*sigmaij2*sigmaij2;
		*uijshift = epsilonij*sigmaij6*(sigmaij6*shift1-1.0); ///< The factor of 4.0 is in shift4
	}

	return 0;
}
