/**
 * Project: mspms2
 * File: frcshort.c
 * 
 * Copyright (C) 2008    Yang Wang <ywangd@gmail.com>
 * Created @ Apr 22, 2008
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

/// Calcualte hard force interactions, bond, angle, dih etc.
int frcshort()
{
	/*
	int ii, isum_atom, isum_mole;
	PSAMPLE_MOLECULE pSampleMole;
	int iSpecie, iMole, iAtom, iabs;
	double uij, virial_ij;
	int iBond, iAngle, iDih, iImp;
	
	ubond = 0.0;
	uangle = 0.0;
	udih = 0.0;
	uimp = 0.0;
	virial_intra = 0.0;
	for (ii=0; ii<natom; ii++)
	{
		fxs[ii] = fys[ii] = fzs[ii] = 0.0;
	}
	
	isum_atom = 0;
	isum_mole = 0;
	// Start of intra-molecule interactions ----------------------------------------------------------------------------------------
	for (iSpecie=0; iSpecie<nspecie; iSpecie++) // specie loop
	{
		pSampleMole = sample_mole + iSpecie; // Get the sample molecule
		for (iMole=isum_mole; iMole<isum_mole+nmole_per_specie[iSpecie]; iMole++) // molecule loop
		{
			// Start of Bond ----------------------------------------------------------------------------
			for (iBond=0; iBond<pSampleMole->nbond; iBond++) // Bond loop
			{
				bndfrc(iSpecie, iBond, isum_atom, &uij, &virial_ij);
				ubond += uij;
				virial_intra += virial_ij;
			} // End of bond loop
			// End of Bond ------------------------------------------------------------------------------
			
			// Start of Angle ---------------------------------------------------------------------------
			for (iAngle=0;iAngle<pSampleMole->nangle;iAngle++) // Angle loop
			{
				aglfrc(iSpecie, iAngle, isum_atom, &uij, &virial_ij);
				uangle += uij;
				virial_intra += virial_ij;
			}
			// End of Angle -----------------------------------------------------------------------------
			
			// Start of Dihedral ------------------------------------------------------------------------
			for (iDih=0;iDih<pSampleMole->ndih;iDih++)
			{
				dihfrc(iSpecie, iDih, isum_atom, &uij, &virial_ij);
				udih += uij;
				virial_intra += virial_ij;
			}
			// End of Dihedral --------------------------------------------------------------------------
			
			// Start of improper ---------------------------------------------------------------------------
			for (iImp=0;iImp<pSampleMole->nimp;iImp++)
			{
				impfrc(iSpecie, iImp, isum_atom, &uij, &virial_ij);
				uimp += uij;
				virial_intra += virial_ij;
			}
			// End of improper ----------------------------------------------------------------------------
			isum_atom += pSampleMole->natom;
		} // End of molecule loop
		isum_mole += nmole_per_specie[iSpecie];
	} // End of specie loop
	
	// Set the intra-molecular energy
	uintra = ubond + uangle + udih + uimp;
	
	
	// printf("ubond=%lf\n",ubond*epsilon_base*RGAS);
	// printf("uangle=%lf\n",uangle*epsilon_base*RGAS);
	// printf("udih=%lf\n",udih*epsilon_base*RGAS);
	// printf("uimp=%lf\n",uimp*epsilon_base*RGAS);
	 
	 */
	
	return 0;
}



