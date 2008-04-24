/**
 * Project: mspms2
 * File: md.c
 * 
 * Copyright (C) 2008    Yang Wang <ywangd@gmail.com>
 * Created @ Nov 29, 2007
 * Modified @ Apr 24, 2008
 * 
 * Description:
 *   Functions related to MD operation go here
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
	int ii, ll, iSpecie, iAtom;
	double r_mass;

	for (ii=0; ii<natom; ii++)
	{
		// Get specie ID and relative atom ID
		get_specie_and_relative_atom_id(ii, &iSpecie, &iAtom);
		r_mass = 1.0/sample_mole[iSpecie].aw[iAtom];
		vx[ii] += (deltby2*fxl[ii]*r_mass);
		vy[ii] += (deltby2*fyl[ii]*r_mass);
		vz[ii] += (deltby2*fzl[ii]*r_mass);
	}

	for (ll=0; ll<nstep_inner; ll++)
	{
		for (ii=0; ii<natom; ii++)
		{
			// Get specie ID and relative atom ID
			get_specie_and_relative_atom_id(ii, &iSpecie, &iAtom);
			r_mass = 1.0/sample_mole[iSpecie].aw[iAtom];
			vx[ii] += (deltsby2*fxs[ii]*r_mass);
			vy[ii] += (deltsby2*fys[ii]*r_mass);
			vz[ii] += (deltsby2*fzs[ii]*r_mass);
			xx[ii] = xx[ii] + delts*vx[ii];
			yy[ii] = yy[ii] + delts*vy[ii];
			zz[ii] = zz[ii] + delts*vz[ii];
		}
		// intra forces, short ranged
		frcshort();
		// compute the pseudo velocity at delts
		for (ii=0; ii<natom; ii++)
		{
			// Get specie ID and relative atom ID
			get_specie_and_relative_atom_id(ii, &iSpecie, &iAtom);
			r_mass = 1.0/sample_mole[iSpecie].aw[iAtom];
			vx[ii] += (deltsby2*fxs[ii]*r_mass);
			vy[ii] += (deltsby2*fys[ii]*r_mass);
			vz[ii] += (deltsby2*fzs[ii]*r_mass);
		}
	} // Inner step loop
	// inter forces, long ranged
	frclong();
	// use the new forces to calculate the new velocities at t+delt
	ukin = 0.0;
	for (ii=0; ii<natom; ii++)
	{
		// Get specie ID and relative atom ID
		get_specie_and_relative_atom_id(ii, &iSpecie, &iAtom);
		r_mass = 1.0/sample_mole[iSpecie].aw[iAtom];
		vx[ii] += (deltby2*fxl[ii]*r_mass);
		vy[ii] += (deltby2*fyl[ii]*r_mass);
		vz[ii] += (deltby2*fzl[ii]*r_mass);
		ukin += (vx[ii]*vx[ii]+vy[ii]*vy[ii]+vz[ii]*vz[ii])/r_mass;
	}
	ukin = ukin*0.5;
	// Calculate instant temperature
	tinst = 2.0*ukin/nfree;
	
	/*
	for (ii=0;ii<natom;ii++)
	{
		printf("%d: fxl=%lf, fyl=%lf, fzl=%lf | fxs=%lf, fys=%lf, fzs=%lf\n",
				ii, fxl[ii]*force_base, fyl[ii]*force_base, fzl[ii]*force_base, 
				fxs[ii]*force_base, fys[ii]*force_base, fzs[ii]*force_base);
	}
	for (ii=0;ii<natom;ii++)
	{
		printf("%d, vx=%lf, vy=%lf, vz=%lf\n", 
				ii, vx[ii]*velocity_base, vy[ii]*velocity_base, vz[ii]*velocity_base);
	}
	*/


	return 0;
}


