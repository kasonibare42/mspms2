/**
 * Project: mspms2
 * File: misc.c
 * 
 * Copyright (C) 2008    Yang Wang <ywangd@gmail.com>
 * Created @ 15/04/2008
 * 
 * Description:
 * 		Anything that does not have a perfect place to go goes here.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <ctype.h>
#include "mspms2.h"

int printit()
{
	fprintf(stdout,"%10d %10.4le %10.4le %10.4le %10.4le %10.4le %10.4le\n",
	istep,utot*epsilon_base*RGAS,upot*epsilon_base*RGAS,ukin*epsilon_base*RGAS,
	tinst*epsilon_base,pinst*pressure_base*1.0e5,boxv*sigma_base*sigma_base*sigma_base);

	fprintf(
			fplog,
			"%10d %10.4le %10.4le %10.4le %10.4le %10.4le %10.4le %10.4le %10.4le %10.4le ",
			istep, utot, upot, ukin, tinst, uinter, uintra, uvdw, ubond, uangle);

	fprintf(
			fplog,
			"%10.4le %10.4le %10.4le %10.4le %10.4le %10.4le %10.4le %10.4le %10.4le ",
			udih, uimp, uewald, usflj, unhts, unhtss, virial, virial_inter,
			virial_intra);

	fprintf(fplog, "%10.4le %10.4le %10.4le\n", utsbs, pinst, boxv);

	return 0;
}

/// Give the absolute atom id, calculate its specie id and relative atom id in the sample
/// molecule.
void get_specie_and_relative_atom_id(int abs_atom_id, int *specie_id, int *sample_atom_id)
{
	int ii, iAtom;
	
	iAtom = 0;
	for (ii=0; ii<nspecie; ii++)
	{
		if (abs_atom_id < iAtom+natom_per_specie[ii])
		{
			// We know this atom is for this specie
			*specie_id = ii;
			*sample_atom_id = (abs_atom_id-iAtom)%sample_mole[ii].natom;
			break;
		}
		iAtom += natom_per_specie[ii];
	}
}

/** 
 * Give the molecule coordination, calculate its center of mass 
 * and inner coordinates. iMole is the absolute molecule ID.
 */
int cal_com_and_efg_one(int iSpecie, int iMole, int iabs, int iAtom)
{
	int ii, jj;
	PSAMPLE_MOLECULE pSampleMole;
	
	pSampleMole = sample_mole + iSpecie;
	// Calculate the center of mass
	mole_xx[iMole] = 0.0;
	mole_yy[iMole] = 0.0;
	mole_zz[iMole] = 0.0;
	for (ii=0;ii<iAtom;ii++)
	{
		jj = iabs + ii;
		mole_xx[iMole] += xx[jj]*pSampleMole->aw[ii];
		mole_yy[iMole] += yy[jj]*pSampleMole->aw[ii];
		mole_zz[iMole] += zz[jj]*pSampleMole->aw[ii];
	}
	mole_xx[iMole] /= pSampleMole->mw;
	mole_yy[iMole] /= pSampleMole->mw;
	mole_zz[iMole] /= pSampleMole->mw;
	
	// Calculate the inner coordinates
	for (ii=0;ii<iAtom;ii++)
	{
		jj = iabs + ii;
		ex[jj] = xx[jj] - mole_xx[iMole];
		fy[jj] = yy[jj] - mole_yy[iMole];
		gz[jj] = zz[jj] - mole_zz[iMole];
	}
	
	return 0;
}

/**
 * Reconstruct the molecule from the position of center of mass.
 * The new coordinates are saved to ex, fy, gz. So, the inner
 * coordinates are no longer kept after this calculation.
 */
int reconstruct_from_com_one(int iMole, int iabs, int iAtom)
{
	int ii;

	for (ii=iabs;ii<iabs+iAtom;ii++)
	{
		ex[ii] = ex[ii] + mole_xx[iMole];
		fy[ii] = fy[ii] + mole_yy[iMole];
		gz[ii] = gz[ii] + mole_zz[iMole];
	}

	return 0;
}


double deriv_inc_gamma(double x, void *params)
{
	return gsl_sf_gamma_inc(0.0, x);
}


int calculate_ljlrc()
{
	int mm, nn;

	// calculate the lrc for the real system
	uljlrc = 0.0;
	pljlrc = 0.0;
	for (mm=0; mm<nspecie; mm++)
	{
		for (nn=0; nn<nspecie; nn++)
		{
			uljlrc += (uljlrc_term[mm][nn]*nmole_per_specie[mm]*nmole_per_specie[nn]/boxv);
			pljlrc += (pljlrc_term[mm][nn]*nmole_per_specie[mm]*nmole_per_specie[nn]/boxv/boxv);
		}
	}

	return 0;
}

