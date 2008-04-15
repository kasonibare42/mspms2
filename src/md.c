/*
 * Functions related to MD operation go here
 *
 * Written by Yang Wang 11-29-2007
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <ctype.h>
#include "mspms2.h"

/// Velocity verlet
int vver()
{
	int ii, ll;

	for (ii=0; ii<natom; ii++)
	{
		// the factor of 1.0e-5 is based on Angstrom (from the force)
		// and femto second (from delt)
		vx[ii] += (deltby2*fxl[ii]*1.0e-5/aw[ii]);
		vy[ii] += (deltby2*fyl[ii]*1.0e-5/aw[ii]);
		vz[ii] += (deltby2*fzl[ii]*1.0e-5/aw[ii]);
	}

	for (ll=0; ll<nstep_inner; ll++)
	{
		for (ii=0; ii<natom; ii++)
		{
			vx[ii] += (deltsby2*fxs[ii]*1.0e-5/aw[ii]);
			vy[ii] += (deltsby2*fys[ii]*1.0e-5/aw[ii]);
			vz[ii] += (deltsby2*fzs[ii]*1.0e-5/aw[ii]);
			xx[ii] = xx[ii] + delts*vx[ii]*1.0e-5;
			yy[ii] = yy[ii] + delts*vy[ii]*1.0e-5;
			zz[ii] = zz[ii] + delts*vz[ii]*1.0e-5;
		}

		// intra forces, short ranged
		rafrc();

		// compute the pseudo velocity at delts
		for (ii=0; ii<natom; ii++)
		{
			vx[ii] += (deltsby2*fxs[ii]*1.0e-5/aw[ii]);
			vy[ii] += (deltsby2*fys[ii]*1.0e-5/aw[ii]);
			vz[ii] += (deltsby2*fzs[ii]*1.0e-5/aw[ii]);
		}
	}

	// inter forces, long ranged
	erfrc();

	// use the new forces to calculate the new velocities at t+delt
	ukin = 0.0;
	for (ii=0; ii<natom; ii++)
	{
		vx[ii] += (deltby2*fxl[ii]*1.0e-5/aw[ii]);
		vy[ii] += (deltby2*fyl[ii]*1.0e-5/aw[ii]);
		vz[ii] += (deltby2*fzl[ii]*1.0e-5/aw[ii]);
		ukin += aw[ii]*(vx[ii]*vx[ii]+vy[ii]*vy[ii]+vz[ii]*vz[ii]);
	}
	ukin = ukin/2.0;

	// calculate instant temperature
	tinst = 2.0*ukin*R_RGAS/nfree;

	return 0;
}


