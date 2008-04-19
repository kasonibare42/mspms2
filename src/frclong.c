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
	int isum_mole, jsum_mole;
	int iStartMole, jStartMole;
	int iEndMole, jEndMole;
	int iSpecie, jSpecie;
	int iMole, jMole;
	int iAtom, jAtom;
	int ii, jj, kk, ll, mm, nn;

	int count = 0;
	
	// STOPPED HERE!!
	
	printf("%d, %d, %d\n",nspecie, nmole, natom);
	
	isum_mole = 0; // Count molecules for iSpecie
	// ---------------------- Start of Species ------------------------------------------------------------------------------
	for (iSpecie=0; iSpecie<nspecie; iSpecie++) // iSpecie
	{
		iAtom = sample_natom_per_mole[iSpecie]; // Number of atoms per molecule for iSpecie
		jsum_mole = isum_mole; // Count molecules for jSpecie
		for (jSpecie=iSpecie; jSpecie<nspecie; jSpecie++) // jSpecie
		{
			jAtom = sample_natom_per_mole[jSpecie]; // Number of atoms per molecules for jSpecie
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
			// ---------------------- Start of Molecules ----------------------------------------------------------------------------
			for (iMole=iStartMole; iMole<iEndMole; iMole++) // molecules in iSpecie
			{  
				// Dynamically assign the j start molecule
				for (jMole=(jStartMole==DYNAMIC_ID?iMole+1:jStartMole); jMole<jEndMole; jMole++) // molecules in jSpecie
				{ 
					// ---------------------- Start of Atoms --------------------------------------------------------------------------------
					for (mm=0; mm<iAtom; mm++) // atoms in iMole (relative position)
					{
						for (nn=0; nn<jAtom; nn++) // atoms in jMole (relative position)
						{
							// printf( "Specie: %d, mole: %d, atom: %d; || Specie: %d, mole: %d, atom: %d\n", iSpecie, iMole, mm, jSpecie, jMole, nn);
							
							count++;
							// printf("%d - %d\t\t||    %d - %d\t%d\n",ii,jj,kk,ll,count);
							

							
							
							 
						} // End of nn  
					} // End of mm
					// ---------------------- End of Atoms ----------------------------------------------------------------------------------
					// printf("%d || %d\n",iMole, jMole);
				} // End of jMole 
			} // End of iMole
			// ---------------------- End of Molecules ------------------------------------------------------------------------------
			jsum_mole += nmole_per_specie[jSpecie]; // Count molecule for jSpecie
		} // End of jSpecie
		isum_mole += nmole_per_specie[iSpecie]; // Count molecule for iSpecie
	} // End of iSpecie
	// ---------------------- End of Species --------------------------------------------------------------------------------


}

