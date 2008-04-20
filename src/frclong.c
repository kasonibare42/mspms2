/**
 * Project: mspms2
 * File: frclong.c
 * 
 * Copyright (C) 2008    Yang Wang <ywangd@gmail.com>
 * Created @ Apr 18, 2008
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

/// Calcualte interaction between atom ii and jj, exclude those from excluding list.
/**
 * @param iStartSpecie is the starting specie for the loop.
 * @param iStartMole is the starting molecule for the loop.
 * @param jStartSpecie is the ending specie for the loop.
 * @param jStartMole is the ending molecule for the loop.\n
 * 
 * These four parameters are used to decide what are the lower and upper limits for ii, jj loop.
 * The value of -1 means all possible values for this parameter. A -1 of the specie overides
 * whatever value for the molecule parameter. \n
 * There are only server possible combination of the parameters.\n
 * 1) (-1,X, -1,X) means loop through all molecules.\n
 * 2) (X,X, -1,X) means loop for a specific molecule and all other molecules.\n
 * 3) (X,X, Y,Y) means calculations between two molecules.\n
 */
int frclong()
{
	int isum_atom, jsum_atom;
	int isum_mole, jsum_mole;
	PSAMPLE_MOLECULE pSampleMole_i, pSampleMole_j;
	int iStartMole, jStartMole;
	int iEndMole, jEndMole;
	int iSpecie, jSpecie;
	int iMole, jMole;
	int iAtom, jAtom;
	int iabs, jabs; // absolute atom id 
	int ii, jj, mm, nn;
	double xxi, yyi, zzi;
	double rxij, ryij, rzij, rijsq;
	double sigmaij, epsilonij, chargeij;

	isum_atom = 0; // Count atoms for iSpecie
	isum_mole = 0; // Count molecules for iSpecie
	// ---------------------- Start of Species ------------------------------------------------------------------------------
	for (iSpecie=0; iSpecie<nspecie; iSpecie++) // iSpecie
	{
		pSampleMole_i = sample_mole + iSpecie;
		iAtom = sample_mole[iSpecie].natom; // Number of atoms per molecule for iSpecie
		jsum_atom = isum_atom;
		jsum_mole = isum_mole; // Count molecules for jSpecie
		for (jSpecie=iSpecie; jSpecie<nspecie; jSpecie++) // jSpecie
		{
			pSampleMole_j = sample_mole + jSpecie;
			jAtom = sample_mole[jSpecie].natom; // Number of atoms per molecules for jSpecie
			// Assign the correct lower and upper bound of i molecule and j moleclue
			if (iSpecie==jSpecie) // If within the same specie, we should have i=(0,n-1) and j=(i+1,n)
			{
				iStartMole = isum_mole;
				iEndMole = isum_mole + nmole_per_specie[iSpecie] - 1;
				jStartMole = DYNAMIC_ID;
				jEndMole = jsum_mole + nmole_per_specie[jSpecie];
			}
			else // If cross species, we have i=(0-n) and j=(n-n+m)
			{
				iStartMole = isum_mole;
				iEndMole = isum_mole + nmole_per_specie[iSpecie];
				jStartMole = jsum_mole;
				jEndMole = jsum_mole + nmole_per_specie[jSpecie];
			}
			iabs = isum_atom; // Set the first absolute atom id for i atom
			// ---------------------- Start of Molecules ----------------------------------------------------------------------------
			for (iMole=iStartMole; iMole<iEndMole; iMole++) // molecules in iSpecie
			{
				// Set the first absolute atom id for j atom dynamically!!
				jabs = jsum_atom + (jStartMole==DYNAMIC_ID ? jAtom : 0);
				// Dynamically assign the j start molecule
				for (jMole=(jStartMole==DYNAMIC_ID ? iMole+1 : jStartMole); jMole
						<jEndMole; jMole++) // molecules in jSpecie
				{
					// ---------------------- Start of Atoms --------------------------------------------------------------------------------
					for (mm=0; mm<iAtom; mm++) // atoms in iMole (relative position)
					{
						if (pSampleMole_i->ghost_type[mm]==GHOST_FULL) // Skip ghost aotm
						{
							continue;
						}
						ii = iabs + mm; // Set the absolute atom id for i atom
						xxi = xx[ii];
						yyi = yy[ii];
						zzi = zz[ii];
						for (nn=0; nn<jAtom; nn++) // atoms in jMole (relative position)
						{
							if (pSampleMole_j->ghost_type[nn]==GHOST_FULL)
							{
								continue;
							}
							jj = jabs + nn; // Set the absolute atom id for j atom
							// printf("Specie: %d, mole: %d, atom: %d, abs aid=%d; || ", iSpecie, iMole, mm, iabs+mm);
							// printf("Specie: %d, mole: %d, atom: %d, abs aid=%d\n", jSpecie, jMole, nn, jabs+nn);
							rxij = xxi - xx[jj];
							ryij = yyi - yy[jj];
							rzij = zzi - zz[jj];
							// minimum image convention
							rxij = rxij - boxlx*rint(rxij/boxlx);
							ryij = ryij - boxly*rint(ryij/boxly);
							rzij = rzij - boxlz*rint(rzij/boxlz);
							rijsq = rxij*rxij + ryij*ryij + rzij*rzij;
							// Lennard-Jones
							if (pSampleMole_i->ghost_type[mm]!=GHOST_LJ&&pSampleMole_j->ghost_type[nn]!=GHOST_LJ&&rijsq<rcutoffsq)
							{
								sigmaij = 0.5*(pSampleMole_i->sigma[mm]+pSampleMole_j->sigma[nn]);
								epsilonij = sqrt(pSampleMole_i->epsilon[mm]*pSampleMole_j->epsilon[nn]);
								// STOPPED HERE!
								// ljfrc();
							}
							
							
							
							

						} // End of nn  
					} // End of mm
					// ---------------------- End of Atoms ----------------------------------------------------------------------------------
					// printf("%d || %d\n",iMole, jMole);
					jabs += jAtom;
				} // End of jMole 
				iabs += iAtom;
			} // End of iMole
			// ---------------------- End of Molecules ------------------------------------------------------------------------------
			jsum_atom += jAtom*nmole_per_specie[jSpecie];
			jsum_mole += nmole_per_specie[jSpecie]; // Count molecule for jSpecie
		} // End of jSpecie
		isum_atom += iAtom*nmole_per_specie[iSpecie];
		isum_mole += nmole_per_specie[iSpecie]; // Count molecule for iSpecie
	} // End of iSpecie
	// ---------------------- End of Species --------------------------------------------------------------------------------

	exit(1);

}


int ljfrc(double rijsq, double sigmaij, double epsilonij, double *uij, double *fij)
{
	double LJswitch;
	double r_rijsq, r_r6, r_r12, r_r12_minus_r_r6;
	double uij_vdw_temp, uij_vdw;
	double sigmaij6;
	
	// if switch potential for LJ is on, then calculate the switch
	if (isLJswitchOn)
	{
		if (rijsq<rcutonsq)
		{
			LJswitch = 1.0;
		}
		else
		{
			LJswitch = (rcutoffsq-rijsq)*(rcutoffsq-rijsq)
					*(rcutoffsq+2.0*rijsq-3.0*rcutonsq)
					/roff2_minus_ron2_cube;
		}
	}
	r_rijsq = sigmaij*sigmaij/rijsq;
	r_r6 = r_rijsq*r_rijsq*r_rijsq;
	r_r12 = r_r6*r_r6;
	r_r12_minus_r_r6 = r_r12 - r_r6;
	uij_vdw_temp = epsilonij*r_r12_minus_r_r6; // still need *4.0
	if (isLJswitchOn) // if switch is used
	{
		uij_vdw += uij_vdw_temp*LJswitch; // still need 4.0
	}
	else
	{
		uij_vdw += uij_vdw_temp; // still need *4.0
	}
	// force calculations
	if (isLJswitchOn)
	{
		if (rijsq<rcutonsq)
		{
			*fij = 24.0*epsilonij*(r_r12_minus_r_r6+r_r12)/rijsq;
		}
		else
		{
			*fij = 24.0*epsilonij*(r_r12_minus_r_r6+r_r12)/rijsq
					*LJswitch // 4.0 for the real energy
					-4.0*uij_vdw_temp*12.0*(rcutoffsq-rijsq)
							*(rcutonsq-rijsq)
							/roff2_minus_ron2_cube;
		}
	}
	else
	{
		*fij = 24.0*epsilonij*(r_r12_minus_r_r6+r_r12)/rijsq;

		// only when LJ swith is not used, we calculate shift energies
		// shift energies
		sigmaij6 = sigmaij*sigmaij*sigmaij*sigmaij*sigmaij*sigmaij;
		gUShiftSession += epsilonij*sigmaij6*(sigmaij6*shift1 -1.0); // still need constant
	}
	
	return 0;
}




