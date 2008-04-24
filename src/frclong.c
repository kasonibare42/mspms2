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

/// Calcualte soft force interactions, including LJ and electrostatic forces.
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
	ushift = 0.0;
	uelec = 0.0;
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
								ushift += uijshift;
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
	bool isSameMole;
	isum_atom = 0;
	isum_mole = 0;
	// Start of intra-molecule interactions ----------------------------------------------------------------------------------------
	for (iSpecie=0; iSpecie<nspecie; iSpecie++) // specie loop
	{
		pSampleMole_i = sample_mole + iSpecie; // Get the sample molecule
		for (iMole=isum_mole; iMole<isum_mole+nmole_per_specie[iSpecie]; iMole++) // molecule loop
		{
			// Start of Bond ----------------------------------------------------------------------------
			for (iBond=0; iBond<pSampleMole_i->nbond; iBond++) // Bond loop
			{
				mm = pSampleMole_i->bnd_idx[iBond][0]; // relative atom id 1
				nn = pSampleMole_i->bnd_idx[iBond][1]; // relative atom id 2
				// Calculate absolute atom id
				ii = isum_atom + mm; // absolute atom id 1
				jj = isum_atom + nn; // absolute atom id 2
				// Calculate separation
				rxij = xx[ii] - xx[jj];
				ryij = yy[ii] - yy[jj];
				rzij = zz[ii] - zz[jj];
				// Check if 1-2 is in the same molecule.
				// 	If they are in the same molecule
				// 		NO LJ. Calculate charge (excl + real(mic))
				// 	If they are NOT in the same molecule
				// 		Check type first!!!
				// 		Calculate necssary LJ and charge (excl + real(mic)) same as above.
				// Check if any separation is out of half box, which means minimum image convetion is required.
				// Check if 1,2 distance > 1,2' distance. If it is true, 1-2, and 1-2' are in the different molecule.
				// Otherwise, they are in the same molecule
				if (fabs(rxij)>boxlx*0.5 || fabs(ryij)>boxly*0.5 || fabs(rzij)
						>boxlz*0.5)
				{
					isSameMole = false;
					// If full ghost atom, skip it.  
					if (pSampleMole_i->ghost_type[mm]==GHOST_FULL
							|| pSampleMole_i->ghost_type[nn]==GHOST_FULL)
					{
						continue;
					}
				}
				else
				{
					isSameMole = true;
				}
				// Lennard-Jones interactions
				if (isSameMole == false)
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
						// Excluding part of ewald summation.
						// This does not need minimum image convention.
						rijsq = rxij*rxij + ryij*ryij + rzij*rzij;
						if (rijsq<rcutoffelecsq)
						{
							ewald_excl_frc(rijsq, chargeij, &uij, &fij);
							uexcl += uij;
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
						}
						// Real part of Ewald summation.
						// This requires minimum image convention
						rxij = rxij - boxlx*rint(rxij/boxlx);
						ryij = ryij - boxly*rint(ryij/boxly);
						rzij = rzij - boxlz*rint(rzij/boxlz);
						rijsq = rxij*rxij + ryij*ryij + rzij*rzij;
						if (rijsq<rcutoffelecsq)
						{
							ewald_real_frc(rijsq, chargeij, &uij, &fij);
							ureal += uij;
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
						if (rijsq<rcutoffelecsq)
						{
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
						}
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
			// End of Bond ------------------------------------------------------------------------------

			// Start of Angle ---------------------------------------------------------------------------
			for (iAngle=0; iAngle<pSampleMole_i->nangle; iAngle++) // Angle loop
			{
				if (pSampleMole_i->isAngle_unique[iAngle]==false) // Only calculate the unique 1,3 pair
				{
					continue;
				}
				mm = pSampleMole_i->agl_idx[iAngle][0]; // relative atom id 1
				nn = pSampleMole_i->agl_idx[iAngle][2]; // relative atom id 2
				// Calculate absolute atom id
				ii = isum_atom + mm; // absolute atom id 1
				jj = isum_atom + nn; // absolute atom id 2
				// Calculate separation
				rxij = xx[ii] - xx[jj];
				ryij = yy[ii] - yy[jj];
				rzij = zz[ii] - zz[jj];
				// Check if 1-3 belong to the same molecule.
				// 	If they are:
				// 		NO LJ! Calculate charge (excl + real(mic)). No ghost type check!
				// 	If they are NOT
				// 		Check ghost type first!!
				// 		Calculate necessary LJ and charge (excl + real(mic)) same as above.
				// Check if any separation is out of half box, which means minimum image convention
				// is needed. Check if 1,3 distance > 1,3' distance. If it is true, they are belong 
				// to different molecule.
				// Otherwise, They are in the same molecule
				if (fabs(rxij)>boxlx*0.5 || fabs(ryij)>boxly*0.5 || fabs(rzij)
						>boxlz*0.5)
				{
					isSameMole = false; // Not in same molecule
					// If full ghost atom, skip it.  
					if (pSampleMole_i->ghost_type[mm]==GHOST_FULL
							|| pSampleMole_i->ghost_type[nn]==GHOST_FULL)
					{
						continue;
					}
				}
				else
				{
					isSameMole = true; // in Same molecule
				}
				// Lennard-Jones interactions
				if (isSameMole == false) // If 1-3' are in different molecules, calculate LJ as required
				{
					fprintf(stderr,"Warning: long 1,3 angle ending pair %d-%d found...\n",ii,jj);
					fprintf(
							fpouts,
							"Warning: long 1,3 angle ending pair %d-%d found...\n",
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
				} // End LJ of 1-3, 1-3' check

				// Electrostatic interactions
				if (iChargeType!=ELECTROSTATIC_NONE)
				{
					chargeij = pSampleMole_i->charge[mm]
							*pSampleMole_i->charge[nn];
					if (iChargeType==ELECTROSTATIC_EWALD) // Ewald summation
					{
						// Excluding part of ewald summation.
						// No minimum image convetion is needed for the excluding part.
						rijsq = rxij*rxij + ryij*ryij + rzij*rzij;
						if (rijsq<rcutoffelecsq)
						{
							ewald_excl_frc(rijsq, chargeij, &uij, &fij);
							uexcl += uij;
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
						}
						// Real part of Ewald summation
						// This requires minimum image convention.
						rxij = rxij - boxlx*rint(rxij/boxlx);
						ryij = ryij - boxly*rint(ryij/boxly);
						rzij = rzij - boxlz*rint(rzij/boxlz);
						rijsq = rxij*rxij + ryij*ryij + rzij*rzij;
						if (rijsq<rcutoffelecsq)
						{
							ewald_real_frc(rijsq, chargeij, &uij, &fij);
							ureal += uij;
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
						if (rijsq<rcutoffelecsq)
						{
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
						}
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
			} // End of Angle loop
			// End of Angle ----------------------------------------------------------------------------

			// Start of Dihedral ------------------------------------------------------------------------
			for (iDih=0; iDih<pSampleMole_i->ndih; iDih++)
			{
				if (pSampleMole_i->isDih_unique[iDih]==false) // Skip it if it is NOT unique
				{
					continue;
				}
				mm = pSampleMole_i->dih_idx[iDih][0]; // relative atom id 1
				nn = pSampleMole_i->dih_idx[iDih][3]; // relative atom id 2
				ii = isum_atom + mm; // absolute atom id 1
				jj = isum_atom + nn; // absolute atom id 2
				// Calculate separation
				rxij = xx[ii] - xx[jj];
				ryij = yy[ii] - yy[jj];
				rzij = zz[ii] - zz[jj];
				// Check if 1-4 belong to the same molecule, which means 1-4 and 1-4' are the same pair.
				// If they are the same pair:
				// 	No type check are needed!! Calculate LJ (with f0) and charge (ONLY real (mic)). 
				// If they are NOT the same pair: (two different molecules)
				// 	Check ghost type!!!
				// 	Calculate LJ (withOUT f0) and charge (ONLY real (mic)).
				// Check if 1,4 distance > 1,4' distance. If it is true, they are not in the same molecule.
				// Otherwise, 1,4 and 1,4' belong to the same molecule.
				if (fabs(rxij)>boxlx*0.5 || fabs(ryij)>boxly*0.5 || fabs(rzij)
						>boxlz*0.5)
				{
					isSameMole == false;
					// If full ghost atom, skip it.  
					if (pSampleMole_i->ghost_type[mm]==GHOST_FULL
							|| pSampleMole_i->ghost_type[nn]==GHOST_FULL)
					{
						continue;
					}
				}
				else
				{
					isSameMole == true;
				}
				// Lennard-Jones interactions are needed for 1,4 no matter if
				// 1-4 and 1-4' are in the same molecule or NOT.
				if (isSameMole == false) // If 1-4' are in different molecules
				{
					fprintf(stderr,"Warning: long 1,4 dihedral ending pair %d-%d found...\n",ii,jj);
					fprintf(
							fpouts,
							"Warning: long 1,4 dihedral ending pair %d-%d found...\n",
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
							<rcutoffsq) // we need ghost check for two molecules
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
				}
				else // If 1-4' are in the same molecule
				{
					// We do not need MIC since we already know 1-4 and 1-4' are the same.
					// isSameMole equal to true here!
					rijsq = rxij*rxij + ryij*ryij + rzij*rzij;
					if (rijsq <rcutoffsq) // we do NOT need ghost check for the same molecule
					{
						sigmaij = 0.5*(pSampleMole_i->sigma[mm]
								+pSampleMole_i->sigma[nn]);
						epsilonij = f0*sqrt(pSampleMole_i->epsilon[mm]
								*pSampleMole_i->epsilon[nn]); // Apply f0 for 1-4 LJ
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
				}

				// Electrostatic interactions
				if (iChargeType!=ELECTROSTATIC_NONE)
				{
					chargeij = pSampleMole_i->charge[mm]
							*pSampleMole_i->charge[nn];
					if (iChargeType==ELECTROSTATIC_EWALD) // Ewald summation
					{
						// Real part of Ewald summation
						rxij = rxij - boxlx*rint(rxij/boxlx);
						ryij = ryij - boxly*rint(ryij/boxly);
						rzij = rzij - boxlz*rint(rzij/boxlz);
						rijsq = rxij*rxij + ryij*ryij + rzij*rzij;
						if (rijsq<rcutoffelecsq)
						{
							ewald_real_frc(rijsq, chargeij, &uij, &fij);
							ureal += uij;
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
					else if (iChargeType==ELECTROSTATIC_SIMPLE_COULOMB) // Simple coulomb
					{
						rxij = rxij - boxlx*rint(rxij/boxlx);
						ryij = ryij - boxly*rint(ryij/boxly);
						rzij = rzij - boxlz*rint(rzij/boxlz);
						rijsq = rxij*rxij + ryij*ryij + rzij*rzij;
						if (rijsq<rcutoffelecsq)
						{
							xatom_coulomb_frc(rijsq, chargeij, &uij, &fij);
							ucoulomb += uij;
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
						}// End of cutoff check
					}
				} // End of electrostatic interaction type check
			} // End of dihedral loop
			// End of Dihedral --------------------------------------------------------------------------

			// Start of non-bonded --------------------------------------------------------------------------
			for (iNbp=0; iNbp<pSampleMole_i->nnbp; iNbp++)
			{
				mm = pSampleMole_i->nbp_idx[iDih][0]; // relative atom id 1
				nn = pSampleMole_i->nbp_idx[iDih][1]; // relative atom id 2
				ii = isum_atom + mm; // absolute atom id 1
				jj = isum_atom + nn; // absolute atom id 2
				// Calculate separation
				rxij = xx[ii] - xx[jj];
				ryij = yy[ii] - yy[jj];
				rzij = zz[ii] - zz[jj];
				// Check if the nbp belongs to the same molecule, which means it is the same pair as
				// its minimum image convention
				// If they are the same pair:
				// 	No type check are needed!! Calculate LJ and charge (ONLY real (mic)). 
				// If they are NOT the same pair: (two different molecules)
				// 	Check ghost type!!!
				// 	Calculate LJ and charge (ONLY real (mic)).
				// Check if real distance > image distance. If it is true, they are not in the same molecule.
				// Otherwise, real and image belong to the same molecule.
				if (fabs(rxij)>boxlx*0.5 || fabs(ryij)>boxly*0.5 || fabs(rzij)
						>boxlz*0.5)
				{
					isSameMole == false;
					// If full ghost atom, skip it.  
					if (pSampleMole_i->ghost_type[mm]==GHOST_FULL
							|| pSampleMole_i->ghost_type[nn]==GHOST_FULL)
					{
						continue;
					}
				}
				else
				{
					isSameMole == true;
				}
				// Lennard-Jones interactions are needed for nbp no matter if
				// the real and image are in the same molecule or NOT.
				if (isSameMole == false) // If image are in different molecules
				{
					fprintf(stderr,"Warning: long non-bonded pair ending pair %d-%d found...\n",ii,jj);
					fprintf(
							fpouts,
							"Warning: long non-bonded pair ending pair %d-%d found...\n",
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
							<rcutoffsq) // we need ghost check for two molecules
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
				}
				else // If real and image are in the same molecule
				{
					// We do not need MIC since we already know 1-4 and 1-4' are the same.
					// isSameMole equal to true here!
					rijsq = rxij*rxij + ryij*ryij + rzij*rzij;
					if (rijsq <rcutoffsq) // we do NOT need ghost check for the same molecule
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
				}
				// Electrostatic interactions
				if (iChargeType!=ELECTROSTATIC_NONE)
				{
					chargeij = pSampleMole_i->charge[mm]
							*pSampleMole_i->charge[nn];
					if (iChargeType==ELECTROSTATIC_EWALD) // Ewald summation
					{
						// Real part of Ewald summation
						rxij = rxij - boxlx*rint(rxij/boxlx);
						ryij = ryij - boxly*rint(ryij/boxly);
						rzij = rzij - boxlz*rint(rzij/boxlz);
						rijsq = rxij*rxij + ryij*ryij + rzij*rzij;
						if (rijsq<rcutoffelecsq)
						{
							ewald_real_frc(rijsq, chargeij, &uij, &fij);
							ureal += uij;
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
					else if (iChargeType==ELECTROSTATIC_SIMPLE_COULOMB) // Simple coulomb
					{
						rxij = rxij - boxlx*rint(rxij/boxlx);
						ryij = ryij - boxly*rint(ryij/boxly);
						rzij = rzij - boxlz*rint(rzij/boxlz);
						rijsq = rxij*rxij + ryij*ryij + rzij*rzij;
						if (rijsq<rcutoffelecsq)
						{
							xatom_coulomb_frc(rijsq, chargeij, &uij, &fij);
							ucoulomb += uij;
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
						}// End of cutoff check
					}
				} // End of electrostatic interaction type check
			} // End of non-bonded pair loop
			// End of non-bonded ----------------------------------------------------------------------------			
			isum_atom += pSampleMole_i->natom;
		} // End of molecule loop 
		isum_mole += nmole_per_specie[iSpecie];
	} // End of specie loop
	
	// Factor for shift energy
	ushift *= shift4;

	// More Ewald and Wolf calculations
	if (iChargeType==ELECTROSTATIC_EWALD)
	{
		int kx, ky, kz, ksq;
		double rkx, rky, rkz, rksq;
		double kvec, sr, si, t, chargei;
		ufourier = 0.0;
		uself = 0.0;
		// Fourier summation of energy
		for (kx=-KMAXX; kx<=KMAXX; kx++) // NOTE: <=
		{
			for (ky=-KMAXY; ky<=KMAXY; ky++) // <= 
			{
				for (kz=-KMAXZ; kz<=KMAXZ; kz++) // <=
				{
					ksq = kx*kx + ky*ky + kz*kz;
					if (ksq<=KSQMAX && ksq!=0)
					{
						rkx = TWOPI_LX*kx;
						rky = TWOPI_LY*ky;
						rkz = TWOPI_LZ*kz;
						rksq = rkx*rkx + rky*rky + rkz*rkz;
						kvec = exp(-Bfactor_ewald*rksq)/rksq;
						// calculate |rho(k)|^2 = sr*sr + si*si
						sr = 0.0;
						si = 0.0;
						// energy
						for (ii=0; ii<natom; ii++)
						{
							// the x,y,z coordinates dont have to be PBC'd
							// before doing the calculation since the existence of the TWOPI_L
							// factor will make the cos, sin functions have the same results
							// with periodical system
							t = rkx*xx[ii] + rky*yy[ii] + rkz*zz[ii];
							// Get the relative atom id
							get_specie_and_relative_atom_id(ii, &iSpecie,
									&iAtom);
							chargei = sample_mole[iSpecie].charge[iAtom];
							sr = sr + chargei*cos(t);
							si = si + chargei*sin(t);
						}
						ufourier += kvec*(sr*sr+si*si);
						// forces
						for (ii=0; ii<natom; ii++)
						{
							t = rkx*xx[ii] + rky*yy[ii] + rkz*zz[ii];
							// Get the relative atom id 
							get_specie_and_relative_atom_id(ii, &iSpecie,
									&iAtom);
							chargei = sample_mole[iSpecie].charge[iAtom];
							fij = 2.0*chargei*(sr*sin(t)-si*cos(t))
									*Vfactor_ewald*kvec*coulomb_prefactor;
							fxl[ii] += fij*rkx;
							fyl[ii] += fij*rky;
							fzl[ii] += fij*rkz;
						}
					} // if (ksq<=KSQMAX && ksq!=0)
				} // KMAXZ
			} // KMAXY
		} // KMAXX
		// Total fourier energy part of ewald
		ufourier *= Vfactor_ewald*coulomb_prefactor;

		// Self interaction corrections, it is constant, no forces.
		for (ii=0; ii<natom; ii++)
		{
			// Get the relative atom id
			get_specie_and_relative_atom_id(ii, &iSpecie, &iAtom);
			chargei = sample_mole[iSpecie].charge[iAtom];
			uself += chargei*chargei;
		}
		uself *= coulomb_prefactor*sqrt(kappa*kappa/pi);

		// Total 3D ewald energy with tinfoil boundary condition
		uewald = ureal + ufourier -uself - uexcl;

		// Calculate additional energy and forces for vaccum boundary condition
		// for ewald summation.
		if (fEwald_BC==EWALD_BC_VACUUM)
		{
			double xxpri, yypri, zzpri; // PBC coords into primal box 
			double qrx, qry, qrz;
			// the atoms have to be grouped adjacent to their respective molecular
			// centers of mass before performing this sum
			// Because any molecules straddling a periodic boundary will cause
			// the total dipole moment of the cell to be greatly exaggerated.
			// -------------------------
			isum_atom = 0;
			isum_mole = 0;
			for (iSpecie=0; iSpecie<nspecie; iSpecie++)
			{
				pSampleMole_i = sample_mole + iSpecie;
				iAtom = pSampleMole_i->natom;
				for (iMole=isum_mole; iMole<isum_mole+nmole_per_specie[iSpecie]; iMole++)
				{
					// calculate the center of mass and relative positions
					cal_com_and_efg_one(iSpecie, iMole, isum_atom, iAtom);
					// calcuate the new molecular center of mass using PBC
					mole_xx[iMole] -= boxlx*rint(mole_xx[iMole]/boxlx);
					mole_yy[iMole] -= boxly*rint(mole_yy[iMole]/boxly);
					mole_zz[iMole] -= boxlz*rint(mole_zz[iMole]/boxlz);
					// calculate the new positions to the PBC'd center of mass
					reconstruct_from_com_one(iMole, isum_atom, iAtom);
					isum_atom += iAtom;
				}
				isum_mole += nmole_per_specie[iSpecie];
			}
			qrx = qry = qrz = 0.0;
			for (ii=0; ii<natom; ii++)
			{
				get_specie_and_relative_atom_id(ii, &iSpecie, &iAtom);
				// use the reconstructed coordinates for this calculations
				// the reconstructed coordinates make sure that the
				// molecular center of mass are in the primal simulation
				// box and the atoms in one molecule are grouped together
				xxpri = ex[ii];
				yypri = fy[ii];
				zzpri = gz[ii];
				chargei = sample_mole[iSpecie].charge[iAtom];
				qrx = qrx + chargei*xxpri;
				qry = qry + chargei*yypri;
				qrz = qrz + chargei*zzpri;
			}
			// energy
			uvacuum = qrx*qrx + qry*qry + qrz*qrz; // still need constant
			uvacuum *= twopi_over_3v*coulomb_prefactor;
			// forces
			for (ii=0; ii<natom; ii++)
			{
				get_specie_and_relative_atom_id(ii, &iSpecie, &iAtom);
				chargei = sample_mole[iSpecie].charge[iAtom];
				fij = -2.0*twopi_over_3v*chargei*coulomb_prefactor;
				fxl[ii] += fij*qrx;
				fyl[ii] += fij*qry;
				fzl[ii] += fij*qrz;
			}
			// TODO: Does this term have additional contribution to the pressure?
			// Or is it already included in the ewald??
			uewald += uvacuum;
		}
		// The virial contribution from ewald summation is the same as its energy
		virial_inter += uewald;
	}
	else if (iChargeType==ELECTROSTATIC_WOLF)
	{
		double chargei;
		uself = 0.0;
		for (ii=0; ii<natom; ii++)
		{
			get_specie_and_relative_atom_id(ii, &iSpecie, &iAtom);
			chargei = sample_mole[iSpecie].charge[iAtom];
			uself += chargei*chargei;
		}
		uself *= coulomb_prefactor*wolfvcon3;
		// Total Wolf energy
		uwolf = ureal - uself - uexcl;
	}

	// Start of Solid-Fluid interactions --------------------------------------------------
	double tforce[3];
	// reset energy
	usflj = 0.0;
	// Virial for solid-fluid with solid fixed is not well defined
	if (iSF_type==SF_NANOTUBE_HYPERGEO)
	{
		for (ii=0; ii<natom; ii++)
		{
			// Get the specie ID and relative atom ID
			get_specie_and_relative_atom_id(ii, &iSpecie, &iAtom);
			if (sample_mole[iSpecie].ghost_type[iAtom]==GHOST_NONE) // If NOT ghost
			{
				cal_sf_hypergeo(ii, iSpecie, iAtom, &uij, &fij);
				usflj += uij;
			} // ghost check
		} // Atom loop
	}
	else if (iSF_type==SF_NANOTUBE_ATOM_EXPLICIT)
	{
		for (ii=0; ii<natom; ii++)
		{
			// Get the specie ID and relative atom ID
			get_specie_and_relative_atom_id(ii, &iSpecie, &iAtom);
			if (sample_mole[iSpecie].ghost_type[iAtom]==GHOST_NONE) // If Not ghost atom
			{
				cal_sf_atom_explicit(ii, iSpecie, iAtom, &uij, &fij);
				usflj += uij;
			} // ghost check
		} // Atom loop
	}
	else if (iSF_type==SF_NANOTUBE_TASOS)
	{
		for (ii=0; ii<natom; ii++)
		{
			// Get the specie ID and relative atom ID
			get_specie_and_relative_atom_id(ii, &iSpecie, &iAtom);
			if (sample_mole[iSpecie].ghost_type[iAtom]==GHOST_NONE) // If Not ghost atom
			{
				// call Tasos's code for force calculations
				call_tasos_forces(sample_mole[iSpecie].interp_type[iAtom], xx[ii], yy[ii], zz[ii], &uij, tforce);
				usflj += uij;
				fxl[ii] += tforce[0];
				fyl[ii] += tforce[1];
				fzl[ii] += tforce[2];
			} // ghost check
		} // natom loop
	}
	else if (iSF_type==SF_NANOTUBE_MY_INTERP)
	{
		for (ii=0; ii<natom; ii++)
		{
			// Get the specie ID and relative atom ID
			get_specie_and_relative_atom_id(ii, &iSpecie, &iAtom);
			if (sample_mole[iSpecie].ghost_type[iAtom]==GHOST_NONE) // If Not ghost atom
			{
				get_values_from_grid(xx[ii], yy[ii], zz[ii], sample_mole[iSpecie].interp_type[iAtom], &uij, tforce);
				usflj += uij*R_RGAS/epsilon_base;
				fxl[ii] += tforce[0]*R_RGAS*sigma_base/epsilon_base;
				fyl[ii] += tforce[1]*R_RGAS*sigma_base/epsilon_base;
				fzl[ii] += tforce[2]*R_RGAS*sigma_base/epsilon_base;
			} // ghost check
		} // Atom loop
	}
	// End of Solid-Fluid interactions ----------------------------------------------------
	
	uotherff = 0.0;
	// Calculate metal cluster energy if necessary.
	if (iExternal_FF_type == FF_DFT_METAL_CLUSTER)
	{
		fnMetalClusterFF(); // Currently, only for the whole system. No single calculation is allowed!
	}

	// Set the inter-molecular energy
	uinter = uvdw + uelec + usflj + uotherff;

	// printf("uvdw = %lf, %lf, %lf\n", uvdw*epsilon_base*RGAS, epsilon_base, RGAS);
	// printf("ureal=%lf, uexcl=%lf\n", ureal*epsilon_base*RGAS, uexcl*epsilon_base*RGAS);
	// printf("ufourier=%lf, uself=%lf\n", ufourier*epsilon_base*RGAS, uself*epsilon_base*RGAS);
	// printf("uvacuum=%lf\n", uvacuum*epsilon_base*RGAS);
	// printf("uelec=%lf\n", uelec*epsilon_base*RGAS);
	// printf("uinter=%lf\n",uinter*epsilon_base*RGAS);
	// printf("virial_inter=%lf\n", virial_inter*epsilon_base*RGAS);

}

