/**
 * Project: mspms2
 * File: elecfrc.c
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

int ewald_real_frc(double rijsq, double chargeij, double *uij, double *fij)
{
	double rij;
	double uij_real_temp;
	double temp1, temp2;

	rij = sqrt(rijsq);

	uij_real_temp = coulomb_prefactor*chargeij/rij;
	*uij = uij_real_temp*erfc(kappa*rij); // real part ewald energy, still need 1/4*pi*epsilon0
	temp1 = kappa*rij;
	temp2 = uij_real_temp*erfc(temp1);
	// real part force calculation
	*fij = (temp2+uij_real_temp*2.0*temp1*exp(-temp1*temp1)/sqrt(pi))/rijsq;

	return 0;
}

int ewald_excl_frc(double rijsq, double chargeij, double *uij, double *fij)
{
	double rij;

	rij = sqrt(rijsq);
	*uij = coulomb_prefactor*chargeij/rij;
	// The negative sign for the fij is because these forces are meant to be
	// substracted from the total forces!
	*fij = -*uij/rijsq; ///< NOTE: negative sign

	return 0;
}

int wolf_real_frc(double rijsq, double chargeij, double *uij, double *fij)
{
	double rij;

	rij = sqrt(rijsq);

	*uij = coulomb_prefactor*chargeij*(erfc(kappa*rij)/rij + wolfvcon1
			+ wolfvcon2*(rij-rcutoffelec));
	*fij = coulomb_prefactor*(chargeij/rij)*(erfc(kappa*rij)/rijsq + wolffcon1
			*exp(-(kappa*rij)*(kappa*rij))/rij + wolffcon2);

	return 0;
}

/**
 * The coulomb interactions between 2 molecules (cross molecules).
 * The coulomb interactions intra-molecular is calculated elsewhere.
 * @iMole is the absolute molecule ID
 */
int xmole_coulomb_frc(int iSpecie, int iMole, int iabs, int iAtom, int jSpecie,
		int jMole, int jabs, int jAtom)
{
	int ii, jj;
	int mm, nn;
	PSAMPLE_MOLECULE pSampleMole_i, pSampleMole_j;
	double rxcm, rycm, rzcm, rcmsq, rcm; // distance between molecular center of mass
	double uij, fij;
	double rxij, ryij, rzij, rijsq, rij; // distance between atoms
	double fxij, fyij, fzij;

	pSampleMole_i = sample_mole + iSpecie;
	pSampleMole_j = sample_mole + jSpecie;
	// first calculate center of mass of the molecules
	cal_com_and_efg_one(iSpecie, iMole, iabs, iAtom);
	cal_com_and_efg_one(jSpecie, jMole, jabs, jAtom);

	// calculate the distance between molecular center of mass
	rxcm = mole_xx[iMole] - mole_xx[jMole];
	rycm = mole_yy[iMole] - mole_yy[jMole];
	rzcm = mole_zz[iMole] - mole_zz[jMole];
	// minimum image convention
	rxcm = rxcm - boxlx*rint(rxcm/boxlx);
	rycm = rycm - boxly*rint(rycm/boxly);
	rzcm = rzcm - boxlz*rint(rzcm/boxlz);
	rcmsq = rxcm*rxcm + rycm*rycm + rzcm*rzcm;
	/** 
	 * The cut off check is based on the molecular center of mass.
	 * There could be different types of the check.
	 * E.g. (1) based on distances of any sites between 2 molecules. 
	 * (2) based on distances of any groups between 2 molecules.
	 */
	// cutoff check
	if (rcmsq < rcutoffelecsq)
	{
		// Distance between two center of mass
		rcm = sqrt(rcmsq);
		// first atom
		for (mm=0; mm<iAtom; mm++)
		{
			ii = iabs + mm; // Get the absolute id of the atom
			if (pSampleMole_i->ghost_type[mm]==GHOST_FULL) // ghost check
			{
				continue;
			}
			// second atom
			for (nn=0; nn<jAtom; nn++)
			{
				jj = jabs + nn; // Get the absolute id of the atom
				if (pSampleMole_j->ghost_type[nn]==GHOST_FULL) // ghost check
				{
					continue;
				}
				// i         j
				//  \       /
				//   m-----n
				//  
				//  (i<-j) = (i<-m) + (m<-n) + (n<-j)
				// where x<-y = x - y, x and y are vectors.
				// Calculate the distance between atoms
				// based on the minimum image convention'd centers of mass
				rxij = ex[ii] + rxcm - ex[jj];
				ryij = fy[ii] + rycm - fy[jj];
				rzij = gz[ii] + rzcm - gz[jj];
				rijsq = rxij*rxij + ryij*ryij + rzij*rzij;
				rij = sqrt(rijsq);
				uij = coulomb_prefactor*pSampleMole_i->charge[mm]
						*pSampleMole_j->charge[nn]/rij;
				ucoulomb += uij;
				virial_inter += uij; // Coulomb virial equal to the energy
				fij = uij/rijsq;
				fxij = fij*rxij;
				fyij = fij*ryij;
				fzij = fij*rzij;
				// force on atom ii
				fxl[ii] += fxij;
				fyl[ii] += fyij;
				fzl[ii] += fzij;
				// force on atom jj
				fxl[jj] -= fxij;
				fyl[jj] -= fyij;
				fzl[jj] -= fzij;
			} // second atom nn
		} // first atom mm
	} // cut off check

	return 0;
}

/**
 * The coulomb interactions between 2 atoms (cross atoms).
 */
int xatom_coulomb_frc(double rijsq, double chargeij, double *uij, double *fij)
{
	int ii, jj;
	double rij;

	rij = sqrt(rijsq);

	*uij = coulomb_prefactor*chargeij/rij;
	*fij = *uij/rijsq;

	return 0;
}

