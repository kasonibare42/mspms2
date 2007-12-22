#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <ctype.h>
#include "random.h"
#include "vars.h"
#include "funcs.h"

int GetNextVacancy(int iSpecie)
{
	return (specie_first_vacancy_idx[iSpecie]==-1) ? nmole
			: specie_first_vacancy_idx[iSpecie];
}

int SetNextVacancy(int iSpecie, int index)
{
	int ii;

	if (specie_first_vacancy_idx[ii] != -1 && index
			> specie_first_vacancy_idx[iSpecie])
	{
		return 0;
	}
	else
	{
		specie_first_vacancy_idx[ii] = index;
	}

	return 0;
}

/// Build the inserted molecule using Sample template
int BulidInsertedMole()
{
	// copy efg coordinates for atoms
	// generate new molecule xyz coordinates
	// ranmar(rndnum, 3);
	// xxnew = pBox->length*(rndnum[0]-0.5);
	// yynew = pBox->width*(rndnum[1]-0.5);
	// zznew = pBox->height*(rndnum[2]-0.5);

	// choose orientation and calculate efg2xyz according to orientation if multi-atom molecule
	/*
	 if (pMoleInsert->natom > 1)
	 {
	 normvec(ex, ey, ez, fx, fy, fz, gx, gy, gz);
	 pMoleInsert->cal_all_atom_efg2xyz_ortn(ex, ey, ez, fx, fy, fz, gx,
	 gy, gz);
	 }
	 else
	 {
	 pMoleInsert->cal_all_atom_efg2xyz();
	 }
	 */

	return 0;
}

int fnInsDelMole()
{
	int ii, jj;
	int iSpecieSelected;
	int iMoleSelected;
	double pcreate, pkill;

	ranmar(rndnum, 2);

	// Find out which specie to insert/delete
	for (ii=0; ii<nspecie; ii++)
	{
		if (rndnum[0] < probability_to_be_selected[ii])
		{
			iSpecieSelected = ii;
			break;
		}
		else
		{
			rndnum[0] = rndnum[0] - probability_to_be_selected[ii];
		}
	}

	// insert or delete
	if (rndnum[1] < probability_to_insert[iSpecieSelected]) // insert
	{
		icounter[26]++; // counter of insertions

		// get the vacancy where the molecule can be inserted into
		iMoleSelected = GetNextVacancy(iSpecieSelected);
		// set the status of the molecule to normal
		mole_status[iMoleSelected] = MOLE_STATUS_NORMAL;

		// use the Sample of this specie to build up the insert molecule
		BulidInsertedMole();

		// Calculate the energy of the inserted molecule
		// del_upoter = 0.0;
		// del_ulj = 0.0;
		// del_ucs = 0.0;
		// del_plj = 0.0;

		// long range correction part
		// del_upoter = del_upoter + del_uljlrc;

		// Accept probability
		// pcreate = exp(-del_upoter*_R_RGAS/tset)*zact*pBox->volume/(pSpecieSelect->nmole+1.0);

		// acceptance?
		ranmar(rndnum, 1);
		if (rndnum[0] < pcreate) // accept
		{
			// number of molecules should be smaller than the NMOLE_MAX

			// accepted insertions
			// accepted insertions for the specie

			// update number of molecule
			// number of molecule in the specie
			// density
			// ideal pressure
			// degree of freedom
			// energies and others
		}
		else // reject
		{
		}
	}
	else // delete
	{
		// counter of reject
		// count of reject for this specie

		// cannot be empty specie
		if (nmole == 0) // empty
		{
		}
		else
		{
			ranmar(rndnum, 1);
			// which molecule to delete

			// calculate the energy for the deleted molecule
			// del_upoter = 0.0;
			// del_ulj = 0.0;
			// del_ucs = 0.0;
			// del_plj = 0.0;

			// long range correction part
			// del_upoter = -del_upoter + del_uljlrc;

			// acceptance probability
			// pkill = exp(-del_upoter*_R_RGAS/tset)*pSpecieSelect->nmole/zact/pBox->volume;

			ranmar(rndnum, 1);
			if (rndnum[0] < pkill) // accept

			{
				// accepted deletions
				// accepted deletions for the specie

				// turn the molecule status of the deleted molecule to vacancy

				// update number of molecule, atoms
				// number of molecule, atom for this specie
				// densities
				// ideal pressure
				// sanity check for empty box: pBox->uljlrc = 0.0;
				// LJ lrc changes
				// degree of freedom
			}
			else // reject
			{
			}
		} // if not empty

	} // if (rndnum[1] < prob_insert)

	return 0;
}
