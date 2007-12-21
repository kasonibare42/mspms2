#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <ctype.h>
#include "vars.h"

int GetNextVacancy(int iSpecie)
{
	return (vacancy_idx[iSpecie]==-1) ? nmole : vacancy_idx[iSpecie];
}

int SetNextVacancy(int iSpecie, int index)
{
	if (vacancy_idx[ii] != -1 && index > vacancy_idx[iSpecie])
	{
		return 0;
	}
	else
	{
		vacancy_idx[ii] = index;
	}
}

int fnInsDelMole()
{
	int ii, jj;
	int iSpecieSelected;
	int iMoleSelected;

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

		// use the molecule next to the last molecule as the inserted molecule
		pMoleInsert = pSpecieSelect->mole_list + pSpecieSelect->nmole;
		pSampleMole = pMoleInsert->pSample_mole;
		for (ii=0; ii<pMoleInsert->natom; ii++)
		{
			pAtomFromInsert = pMoleInsert->atom_list + ii;
			pAtomFromSample = pSampleMole->atom_list + ii;
			pAtomFromInsert->ee = pAtomFromSample->ee;
			pAtomFromInsert->ff = pAtomFromSample->ff;
			pAtomFromInsert->gg = pAtomFromSample->gg;
		}
		ranmar(rndnum, 3);
		xxnew = pBox->length*(rndnum[0]-0.5);
		yynew = pBox->width*(rndnum[1]-0.5);
		zznew = pBox->height*(rndnum[2]-0.5);
		*(pMoleInsert->pXx) = xxnew;
		*(pMoleInsert->pYy) = yynew;
		*(pMoleInsert->pZz) = zznew;
		// choose orientation and calculate efg2xyz according to orientation if multi-atom molecule
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
		// calculate the energy of the inserted molecule
		del_upoter = 0.0;
		del_ulj = 0.0;
		del_ucs = 0.0;
		del_plj = 0.0;
		for (ii=0; ii<pBox->nspecie; ii++)
		{
			pSpecie = pBox->specie_list + ii;
			for (jj=0; jj<pSpecie->nmole; jj++)
			{
				pMole = pSpecie->mole_list + jj;
				ljfrc(pMoleInsert, pMole, pBox, ulj, ucs, plj, false);
				del_upoter += ulj;
				del_ulj += ulj;
				del_ucs += ucs;
				del_plj += plj;
			}
		}
		// calculate the solid-fluid energy of the inserted molecule if nanotube exists
		if (pBox->nnanotube)
		{
			sfljfrc(pMoleInsert, pBox, del_usflj, del_usfljcs, del_psflj, false);
			del_upoter += del_usflj;
		}
		// long range correction part
		del_uljlrc = 0.0;
		del_pljlrc = 0.0;
		if (input_data.bLrc)
			deln_ljlrc(pBox, del_uljlrc, del_pljlrc, idSpecieSelect, 1); // calculate the change of ljlrc
		del_upoter = del_upoter + del_uljlrc;
		pcreate = exp(-del_upoter*_R_RGAS/tset)*zact*pBox->volume
				/(pSpecieSelect->nmole+1.0);
		ranmar(rndnum, 1);
		if (rndnum[0] < pcreate) // accept
		{
			if (pSpecieSelect->nmole == pSpecieSelect->nmole_max)
			{
				std::cout << "ERROR: Too many molecules in specie " << idSpecieSelect << ".\n";
				exit(1);
			}
			pBox->counter[7]++; // accepted insertions
			pSpecieSelect->counter[1]++; // accepted insertions for the specie
			pSpecieSelect->nmole += 1;
			pSpecieSelect->natom += pMoleInsert->natom;
			box_nmole_old = pBox->nmole;
			pBox->nmole = box_nmole_new = box_nmole_old + 1;
			pBox->natom += pMoleInsert->natom;
			pBox->weight += pMoleInsert->weight;
			pBox->rho_atom = pBox->natom/pBox->volume;
			pBox->rho_mole = pBox->nmole/pBox->volume;
			pBox->pideal = pBox->rho_atom*pBox->temperature*_KB_OVER_ANG3;
			pBox->ulj += del_ulj;
			pBox->uljcs += del_ucs;
			pBox->plj += del_plj;
			if (pBox->nnanotube)
			{
				pBox->usflj += del_usflj;
				pBox->usfljcs += del_usfljcs;
				pBox->psflj += del_psflj;
			}
			pBox->uljlrc = pBox->uljlrc*box_nmole_old/box_nmole_new + del_uljlrc/box_nmole_new; // del_uljlrc is for the whole system
			pBox->pljlrc += del_pljlrc;
			pBox->rfree = 1.0/(3.0*pBox->natom-3.0);
		}
		else // reject

		{}
	}
	else // delete

	{
		pBox->counter[8]++;
		pSpecieSelect->counter[2]++;
		if (pSpecieSelect->nmole == 0) // empty

		{}
		else
		{
			ranmar(rndnum, 1);
			idMoleDelete = int(rndnum[0]*pSpecieSelect->nmole);
			pMoleDelete = pSpecieSelect->mole_list + idMoleDelete;
			// calculate the energy for the deleted molecule
			del_upoter = 0.0;
			del_ulj = 0.0;
			del_ucs = 0.0;
			del_plj = 0.0;
			for (ii=0;ii<pBox->nspecie;ii++)
			{
				pSpecie = pBox->specie_list + ii;
				for (jj=0;jj<pSpecie->nmole;jj++)
				{
					if (ii!=idSpecieSelect || jj!=idMoleDelete)
					{
						pMole = pSpecie->mole_list + jj;
						ljfrc(pMoleDelete, pMole, pBox, ulj, ucs, plj, false);
						del_upoter += ulj;
						del_ulj += ulj;
						del_ucs += ucs;
						del_plj += plj;
					}
				}
			}
			// calculate the solid-fluid energy of the inserted molecule if nanotube exists
			if (pBox->nnanotube)
			{
				sfljfrc(pMoleDelete, pBox, del_usflj, del_usfljcs, del_psflj, false);
				del_upoter += del_usflj;
			}
			// long range correction part
			del_uljlrc = 0.0;
			del_pljlrc = 0.0;
			if (input_data.bLrc)
			deln_ljlrc(pBox, del_uljlrc, del_pljlrc, idSpecieSelect, -1);
			del_upoter = -del_upoter + del_uljlrc;
			pkill = exp(-del_upoter*_R_RGAS/tset)*pSpecieSelect->nmole/zact/pBox->volume;
			ranmar(rndnum, 1);
			if (rndnum[0] < pkill) // accept

			{
				pBox->counter[9]++; // accepted deletions
				pSpecieSelect->counter[3]++; // accepted deletions for the specie
				// get the address of the last molecule
				pMoleLast = pSpecieSelect->mole_list + (pSpecieSelect->nmole-1);
				// copy the last molecule to the positon of the deleted molecule
				*(pMoleDelete->pXx) = *(pMoleLast->pXx);
				*(pMoleDelete->pYy) = *(pMoleLast->pYy);
				*(pMoleDelete->pZz) = *(pMoleLast->pZz);
				for (ii=0;ii<pMoleLast->natom;ii++)
				{
					pAtomFromLast = pMoleLast->atom_list + ii;
					pAtomFromDelete = pMoleDelete->atom_list + ii;
					pAtomFromDelete->xx = pAtomFromLast->xx;
					pAtomFromDelete->yy = pAtomFromLast->yy;
					pAtomFromDelete->zz = pAtomFromLast->zz;
					pAtomFromDelete->xdisp = pAtomFromLast->xdisp;
					pAtomFromDelete->ydisp = pAtomFromLast->ydisp;
					pAtomFromDelete->zdisp = pAtomFromLast->zdisp;
					// do not need to copy the velocities and forces
					// since v and f will be generated or calculated
					// before a new dispmd cycle for a HMC run
				}
				pSpecieSelect->nmole -= 1;
				pSpecieSelect->natom -= pMoleDelete->natom;
				box_nmole_old = pBox->nmole;
				pBox->nmole = box_nmole_new = box_nmole_old - 1;
				pBox->natom -= pMoleDelete->natom;
				pBox->weight -= pMoleDelete->weight;
				pBox->rho_atom = pBox->natom/pBox->volume;
				pBox->rho_mole = pBox->nmole/pBox->volume;
				pBox->pideal = pBox->rho_atom*pBox->temperature*_KB_OVER_ANG3;
				pBox->ulj -= del_ulj;
				pBox->uljcs -= del_ucs;
				pBox->plj -= del_plj;
				if (pBox->nnanotube)
				{
					pBox->usflj -= del_usflj;
					pBox->usfljcs -= del_usfljcs;
					pBox->psflj -= del_psflj;
				}
				if (box_nmole_new==0) // sanity check for empty box
				pBox->uljlrc = 0.0;
				else
				pBox->uljlrc = pBox->uljlrc*box_nmole_old/box_nmole_new + del_uljlrc/box_nmole_new;
				pBox->pljlrc += del_pljlrc;
				pBox->rfree = 1.0/(3.0*pBox->natom-3.0);
			}
			else // reject

			{}
		} // if not empty

	} // if (rndnum[1] < prob_insert)

	return 0;
}
