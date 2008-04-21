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
	double fxi, fyi, fzi;
	double rxij, ryij, rzij, rijsq;
	double sigmaij, epsilonij, chargeij;
	double uij, fij, uijshift;
	double fxij, fyij, fzij;

	uvdw = 0.0;
	ureal = 0.0;
	uexcl = 0.0;
	ucoulomb = 0.0;
	virial_inter = 0.0;
	for (ii=0; ii<natom; ii++)
	{
		fxl[ii] = fyl[ii] = fzl[ii] = 0.0;
	}

	// Start of inter-molecule interactions --------------------------------------------------------------------------------------
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
						fxi = fxl[ii];
						fyi = fyl[ii];
						fzi = fzl[ii];
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
							if (pSampleMole_i->ghost_type[mm]!=GHOST_LJ
									&&pSampleMole_j->ghost_type[nn]!=GHOST_LJ
									&&rijsq<rcutoffsq)
							{
								sigmaij = 0.5*(pSampleMole_i->sigma[mm]
										+pSampleMole_j->sigma[nn]);
								epsilonij = sqrt(pSampleMole_i->epsilon[mm]
										*pSampleMole_j->epsilon[nn]);
								ljfrc(rijsq, sigmaij, epsilonij, &uij, &fij,
										&uijshift);
								uvdw += uij;
								virial_inter += fij*rijsq;
								fxij = fij*rxij;
								fyij = fij*ryij;
								fzij = fij*rzij;
								// forces on atom ii
								fxi += fxij;
								fyi += fyij;
								fzi += fzij;
								// forces on atom jj
								fxl[jj] -= fxij;
								fyl[jj] -= fyij;
								fzl[jj] -= fzij;
							} // End of LJ calculations

							// Electrostatic interactions
							if (iChargeType!=ELECTROSTATIC_NONE)
							{
								chargeij = pSampleMole_i->charge[mm]
										*pSampleMole_j->charge[nn];
								if (iChargeType==ELECTROSTATIC_EWALD) // Ewald summation
								{
									if (rijsq<rcutoffelecsq)
									{
										ewald_real_frc(rijsq, chargeij, &uij,
												&fij);
										ureal += uij;
										virial_inter += fij*rijsq;
										fxij = fij*rxij;
										fyij = fij*ryij;
										fzij = fij*rzij;
										// forces on atom ii
										fxi += fxij;
										fyi += fyij;
										fzi += fzij;
										// forces on atom jj
										fxl[jj] -= fxij;
										fyl[jj] -= fyij;
										fzl[jj] -= fzij;
									} // End of cutoff check
								}
								else if (iChargeType==ELECTROSTATIC_WOLF) // Wolf method
								{
									if (rijsq<rcutoffelecsq)
									{
										wolf_real_frc(rijsq, chargeij, &uij,
												&fij);
										ureal += uij;
										virial_inter += fij*rijsq;
										fxij = fij*rxij;
										fyij = fij*ryij;
										fzij = fij*rzij;
										// forces on atom ii
										fxi += fxij;
										fyi += fyij;
										fzi += fzij;
										// forces on atom jj
										fxl[jj] -= fxij;
										fyl[jj] -= fyij;
										fzl[jj] -= fzij;
									} // End of cutoff check
								}
							} // End of eletrostaic interaction type check
						} // End of nn atom (jj)
						fxl[ii] = fxi;
						fyl[ii] = fyi;
						fzl[ii] = fzi;
					} // End of mm atom (ii)
					// ---------------------- End of Atoms ----------------------------------------------------------------------------------
					// printf("%d || %d\n",iMole, jMole);
					if (iChargeType==ELECTROSTATIC_SIMPLE_COULOMB)
					{
						xmole_coulomb_frc(iSpecie, iMole, iabs, iAtom, jSpecie,
								jMole, jabs, jAtom);
					}
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
	// End of inter-molecule interactions ----------------------------------------------------------------------------------------

	int iBond, iAngle, iDih, iNbp;
	int rxij_old, ryij_old, rzij_old;
	isum_atom = 0;
	isum_mole = 0;
	// Start of inter-molecule interactions ----------------------------------------------------------------------------------------
	for (iSpecie=0; iSpecie<nspecie; iSpecie++) // specie loop
	{
		pSampleMole_i = sample_mole + iSpecie; // Get the sample molecule
		for (iMole=0; iMole<nmole_per_specie[iSpecie]; iMole++) // molecule loop
		{
			for (iBond=0; iBond<pSampleMole_i->nbond; iBond++) // Bond loop
			{
				mm = pSampleMole_i->bnd_idx[iBond][0]; // relative atom id 1
				nn = pSampleMole_i->bnd_idx[iBond][1]; // relative atom id 2
				ii = isum_atom + mm; // absolute atom id 1
				jj = isum_atom + nn; // absolute atom id 2
				// If full ghost atom, skip it.
				if (pSampleMole_i->ghost_type[mm]==GHOST_FULL
						|| pSampleMole_i->ghost_type[nn]==GHOST_FULL)
				{
					continue;
				}
				rxij = xx[ii] - xx[jj];
				ryij = yy[ii] - yy[jj];
				rzij = zz[ii] - zz[jj];
				// Lennard-Jones interactions
				// Check if 1,2 distance < 1,2' distance. If it is true, no LJ needed.
				// Otherwise, calculate LJ between 1 and 2'.
				if (fabs(rxij)>boxlx*0.5 || fabs(ryij)>boxly*0.5 || fabs(rzij)
						>boxlz*0.5)
				{
					fprintf(stderr,"Warning: long 1,2 bond ending pair %d-%d found...\n",ii,jj);
					fprintf(
							fpouts,
							"Warning: long 1,2 bond ending pair %d-%d found...\n",
							ii, jj);
					// Save the values before MIC
					rxij_old = rxij;
					ryij_old = ryij;
					rzij_old = rzij;
					// Minimum image connvention
					rxij = rxij - boxlx*rint(rxij/boxlx);
					ryij = ryij - boxly*rint(ryij/boxly);
					rzij = rzij - boxlz*rint(rzij/boxlz);
					rijsq = rxij*rxij + ryij*ryij + rzij*rzij;
					if (pSampleMole_i->ghost_type[mm]!=GHOST_LJ
							&& pSampleMole_i->ghost_type[nn]!=GHOST_LJ && rijsq
							<rcutoffsq)
					{
						sigmaij = 0.5*(pSampleMole_i->sigma[mm]
								+pSampleMole_i->sigma[nn]);
						epsilonij = sqrt(pSampleMole_i->epsilon[mm]
								*pSampleMole_i->epsilon[nn]);
						ljfrc(rijsq, sigmaij, epsilonij, &uij, &fij, &uijshift);
						uvdw += uij;
						virial_inter += fij*rijsq;
						fxij = fij*rxij;
						fyij = fij*ryij;
						fzij = fij*rzij;
						// forces on atom ii
						fxl[ii] += fxij;
						fyl[ii] += fyij;
						fzl[ii] += fzij;
						// forces on atom jj
						fxl[jj] -= fxij;
						fyl[jj] -= fyij;
						fzl[jj] -= fzij;
					} // End of ghost atom and rijsq cutoff check
					// Restore the old values for future Ewald interactions
					rxij = rxij_old;
					ryij = ryij_old;
					rzij = rzij_old;
				} // End LJ of 1-2, 1-2' check

				// Electrostatic interactions
				if (iChargeType!=ELECTROSTATIC_NONE)
				{
					chargeij = pSampleMole_i->charge[mm]
							*pSampleMole_i->charge[nn];
					if (iChargeType==ELECTROSTATIC_EWALD) // Ewald summation
					{
						// Excluding part of ewald summation
						rijsq = rxij*rxij + ryij*ryij + rzij*rzij;
						ewald_excl_frc(rijsq, chargeij, &uij, &fij);
						uexcl += uij;
						virial_inter += fij*rijsq;
						fxij = fij*rxij;
						fyij = fij*ryij;
						fzij = fij*rzij;
						// forces on atom ii
						fxl[ii] += fxij;
						fyl[ii] += fyij;
						fzl[ii] += fzij;
						// forces on atom jj
						fxl[jj] -= fxij;
						fyl[jj] -= fyij;
						fzl[jj] -= fzij;
						// Real part of Ewald summation
						rxij = rxij - boxlx*rint(rxij/boxlx);
						ryij = ryij - boxly*rint(ryij/boxly);
						rzij = rzij - boxlz*rint(rzij/boxlz);
						rijsq = rxij*rxij + ryij*ryij + rzij*rzij;
						if (rijsq<rcutoffelecsq)
						{
							ewald_real_frc(rijsq, chargeij, &uij, &fij);
							ureal += uij;
							virial_inter += fij*rijsq;
							fxij = fij*rxij;
							fyij = fij*ryij;
							fzij = fij*rzij;
							// forces on atom ii
							fxl[ii] += fxij;
							fyl[ii] += fyij;
							fzl[ii] += fzij;
							// forces on atom jj
							fxl[jj] -= fxij;
							fyl[jj] -= fyij;
							fzl[jj] -= fzij;
						} // Cutoff check
					}
					else if (iChargeType==ELECTROSTATIC_WOLF) // Wolf method
					{
						// Excluding part of Wolf method
						rijsq = rxij*rxij + ryij*ryij + rzij*rzij;
						ewald_excl_frc(rijsq, chargeij, &uij, &fij);
						uexcl += uij;
						virial_inter += fij*rijsq;
						fxij = fij*rxij;
						fyij = fij*ryij;
						fzij = fij*rzij;
						// forces on atom ii
						fxl[ii] += fxij;
						fyl[ii] += fyij;
						fzl[ii] += fzij;
						// forces on atom jj
						fxl[jj] -= fxij;
						fyl[jj] -= fyij;
						fzl[jj] -= fzij;
						// Real part of Wolf method
						rxij = rxij - boxlx*rint(rxij/boxlx);
						ryij = ryij - boxly*rint(ryij/boxly);
						rzij = rzij - boxlz*rint(rzij/boxlz);
						rijsq = rxij*rxij + ryij*ryij + rzij*rzij;
						if (rijsq<rcutoffelecsq)
						{
							wolf_real_frc(rijsq, chargeij, &uij, &fij);
							ureal += uij;
							virial_inter += fij*rijsq;
							fxij = fij*rxij;
							fyij = fij*ryij;
							fzij = fij*rzij;
							// forces on atom ii
							fxl[ii] += fxij;
							fyl[ii] += fyij;
							fzl[ii] += fzij;
							// forces on atom jj
							fxl[jj] -= fxij;
							fyl[jj] -= fyij;
							fzl[jj] -= fzij;
						} // Cutoff check
					}
				} // End of electrostatic interaction type check
			} // End of bond loop

			for (iAngle=0; iAngle<pSampleMole_i->nangle; iAngle++) // Angle loop
			{

			}

			isum_atom += pSampleMole_i->natom;
			isum_mole++;
		} // End of molecule loop
	} // End of specie loop


	printf("uvdw = %lf, %lf, %lf\n", uvdw*epsilon_base*RGAS, epsilon_base, RGAS);
	printf("ureal = %lf\n", ureal*epsilon_base*RGAS);
	printf("ucoulomb = %lf\n", ucoulomb*epsilon_base*RGAS);
	exit(1);

}

