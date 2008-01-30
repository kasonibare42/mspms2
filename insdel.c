#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <ctype.h>
#include "random.h"
#include "vars.h"
#include "funcs.h"

int normvec()
{
	int bRegenerate;
	double rndnum[3];
	double module, rmodulesqrt, p1, p2;

	// Gram-Schmidt vector orthonormalization process
	do
	{
		bRegenerate = false;
		// pick up an orientation
		ranmar(rndnum, 3);
		fOrientE_x = rndnum[0] - 0.5;
		fOrientE_y = rndnum[1] - 0.5;
		fOrientE_z = rndnum[2] - 0.5;
		ranmar(rndnum, 3);
		fOrientF_x = rndnum[0] - 0.5;
		fOrientF_y = rndnum[1] - 0.5;
		fOrientF_z = rndnum[2] - 0.5;
		ranmar(rndnum, 3);
		fOrientG_x = rndnum[0] - 0.5;
		fOrientG_y = rndnum[1] - 0.5;
		fOrientG_z = rndnum[2] - 0.5;
		module = fOrientE_x*fOrientE_x+fOrientE_y*fOrientE_y+fOrientE_z
				*fOrientE_z;
		if (module < 10.0e-3)
		{
			bRegenerate = true;
		}
		else
		{
			rmodulesqrt = 1.0/sqrt(module);
			fOrientE_x *= rmodulesqrt;
			fOrientE_y *= rmodulesqrt;
			fOrientE_z *= rmodulesqrt;
			module = 0.0;
			p1 = fOrientE_x*fOrientF_x + fOrientE_y*fOrientF_y + fOrientE_z
					*fOrientF_z;
			fOrientF_x = fOrientF_x - fOrientE_x*p1;
			fOrientF_y = fOrientF_y - fOrientE_y*p1;
			fOrientF_z = fOrientF_z - fOrientE_z*p1;
			module = fOrientF_x*fOrientF_x+fOrientF_y*fOrientF_y+fOrientF_z
					*fOrientF_z;
			if (module<10.0e-3)
			{
				bRegenerate = true;
			}
			else
			{
				rmodulesqrt = 1.0/sqrt(module);
				fOrientF_x *= rmodulesqrt;
				fOrientF_y *= rmodulesqrt;
				fOrientF_z *= rmodulesqrt;
				module = 0.0;
				p1 = fOrientE_x*fOrientG_x + fOrientE_y*fOrientG_y + fOrientE_z
						*fOrientG_z;
				p2 = fOrientF_x*fOrientG_x + fOrientF_y*fOrientG_y + fOrientF_z
						*fOrientG_z;
				fOrientG_x = fOrientG_x - fOrientE_x*p1 - fOrientF_x*p2;
				fOrientG_y = fOrientG_y - fOrientE_y*p1 - fOrientF_y*p2;
				fOrientG_z = fOrientG_z - fOrientE_z*p1 - fOrientF_z*p2;
				module = fOrientG_x*fOrientG_x + fOrientG_y*fOrientG_y
						+ fOrientG_z*fOrientG_z;
				if (module < 10.0e-3)
				{
					bRegenerate = true;
				}
				else
				{
					rmodulesqrt = 1.0/sqrt(module);
					fOrientG_x *= rmodulesqrt;
					fOrientG_y *= rmodulesqrt;
					fOrientG_z *= rmodulesqrt;
				}
			}
		}
	} while (bRegenerate);

	return 0;
}

int cal_all_atom_efg2xyz_ortn(int iSpecieSelected, int iMoleSelected,
		double xxnew, double yynew, double zznew)
{
	int ii;
	int idFirstRealAtom;
	int idRealAtom;

	idFirstRealAtom = mole_first_atom_idx[iMoleSelected];
	for (ii=sample_mole_first_atom_idx[iSpecieSelected]; ii
			<sample_mole_last_atom_idx[iSpecieSelected]; ii++)
	{
		idRealAtom = idFirstRealAtom + ii;
		// e2x
		xx[idRealAtom] = fOrientE_x*sample_ee[ii] + xxnew;
		xx[idRealAtom] += fOrientF_x*sample_ff[ii];
		xx[idRealAtom] += fOrientG_x*sample_gg[ii];
		// f2y
		yy[idRealAtom] = fOrientE_y*sample_ee[ii] + yynew;
		yy[idRealAtom] += fOrientF_y*sample_ff[ii];
		yy[idRealAtom] += fOrientG_y*sample_gg[ii];
		// g2z
		zz[idRealAtom] = fOrientE_z*sample_ee[ii] + zznew;
		zz[idRealAtom] += fOrientF_z*sample_ff[ii];
		zz[idRealAtom] += fOrientG_z*sample_gg[ii];
	}

	return 0;
}

int cal_all_atom_efg2xyz(int iSpecieSelected, int iMoleSelected, double xxnew,
		double yynew, double zznew)
{
	// always use (1,0,0), (0,1,0), (0,0,1) as the orientation vectors.
	// so the calculations are reduced as follows
	int ii;
	int idFirstRealAtom;
	int idRealAtom;

	idFirstRealAtom = mole_first_atom_idx[iMoleSelected];
	for (ii=sample_mole_first_atom_idx[iSpecieSelected]; ii
			<sample_mole_last_atom_idx[iSpecieSelected]; ii++)
	{
		idRealAtom = idFirstRealAtom + ii;

		xx[idRealAtom] = sample_ee[ii] + xxnew;
		yy[idRealAtom] = sample_ff[ii] + yynew;
		zz[idRealAtom] = sample_gg[ii] + zznew;
	}

	return 0;
}

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
int BulidInsertedMole(int iSpecieSelected, int iMoleSelected)
{
	double xxnew, yynew, zznew;

	// Generate the new xyz coordinates for the molecule 
	// Default coordinates of the molecule is the first atom of the molecule
	ranmar(rndnum, 3);
	xxnew = boxlx*(rndnum[0]-0.5);
	yynew = boxly*(rndnum[1]-0.5);
	zznew = boxlz*(rndnum[2]-0.5);

	// copy efg coordinates for atoms

	// choose orientation and calculate efg2xyz according to orientation if multi-atom molecule
	if (sample_natom_per_mole[iSpecieSelected] > 1)
	{
		// If the multiple atom molecule, we need to generate a random orientation for it
		normvec();
		cal_all_atom_efg2xyz_ortn(iSpecieSelected, iMoleSelected, xxnew, yynew,
				zznew);
	}
	else // single atom, spherical molecule does not need orientation
	{
		cal_all_atom_efg2xyz(iSpecieSelected, iMoleSelected, xxnew, yynew,
				zznew);
	}

	return 0;
}

/// Initialize the inserted molecule if it is NOT initialized
int InitInsertedMole(int iSpecieSelected, int iMoleSelected)
{
	int ii, jj, kk;
	int iAtom, iMole, iBond, iAngle, iDih, iImp, iNbp;
	int idMoleInSpecie;

	// init the real atom, bond, angle, dihedral, improper, non-bonded list
	iAtom = idAtomUninit;
	iMole = iMoleSelected;
	iBond = idBondUninit;
	iAngle = idAngleUninit;
	iDih = idDihUninit;
	iImp = idImpUninit;
	iNbp = idNbpUninit;
	
	// specie_first_atom_idx[iSpecieSelected] = iAtom;
	// specie_first_mole_idx[iSpecieSelected] = iMole;
	idMoleInSpecie = nmole_per_specie[iSpecieSelected];
	nmole_per_specie[iSpecieSelected]++;

	mole2specie[iMole] = iSpecieSelected;
	mole_status[iMole] = MOLE_STATUS_NORMAL;
	iPhysicalMoleIDFromMetaID[iMole] = iMole;
	iPhysicalMoleIDFromMetaIDinSpecie[iSpecieSelected][idMoleInSpecie] = iMole;

	// atom
	mole_first_atom_idx[iMole] = iAtom;
	for (kk=sample_mole_first_atom_idx[iSpecieSelected]; kk
			<sample_mole_last_atom_idx[iSpecieSelected]; kk++)
	{
		atom2mole[iAtom] = iMole;
		strcpy(atomname[iAtom], sample_atomname[kk]);
		aw[iAtom] = sample_aw[kk];
		epsilon[iAtom] = sample_epsilon[kk];
		sigma[iAtom] = sample_sigma[kk];
		charge[iAtom] = sample_charge[kk];
		isghost[iAtom] = sample_isghost[kk];
		tasostype[iAtom] = sample_tasostype[kk];
		iAtom++;
	}
	mole_last_atom_idx[iMole] = iAtom;

	// bond
	mole_first_bond_idx[iMole] = iBond;
	for (kk=sample_mole_first_bond_idx[iSpecieSelected]; kk
			<sample_mole_last_bond_idx[iSpecieSelected]; kk++)
	{
		bond_idx[iBond][0] = sample_bond_idx[kk][0] + idMoleInSpecie
			*sample_natom_per_mole[iSpecieSelected];
		bond_idx[iBond][1] = sample_bond_idx[kk][1] + idMoleInSpecie
			*sample_natom_per_mole[iSpecieSelected];
		bond_type[iBond] = sample_bond_type[kk];
		Kb[iBond] = sample_Kb[kk];
		Req[iBond] = sample_Req[kk];
		alpha[iBond] = sample_alpha[kk];
		iBond++;
	}
	mole_last_bond_idx[iMole] = iBond;

	// angle
	mole_first_angle_idx[iMole] = iAngle;
	for (kk=sample_mole_first_angle_idx[iSpecieSelected]; kk
			<sample_mole_last_angle_idx[iSpecieSelected]; kk++)
	{
		angle_idx[iAngle][0] = sample_angle_idx[kk][0] + idMoleInSpecie
			*sample_natom_per_mole[iSpecieSelected];
		angle_idx[iAngle][1] = sample_angle_idx[kk][1] + idMoleInSpecie
			*sample_natom_per_mole[iSpecieSelected];
		angle_idx[iAngle][2] = sample_angle_idx[kk][2] + idMoleInSpecie
			*sample_natom_per_mole[iSpecieSelected];
		angle_type[iAngle] = sample_angle_type[kk];
		Ktheta[iAngle] = sample_Ktheta[kk];
		Thetaeq[iAngle] = sample_Thetaeq[kk];
		agl_para_3[iAngle] = sample_agl_para_3[kk];
		agl_para_4[iAngle] = sample_agl_para_4[kk];
		agl_para_5[iAngle] = sample_agl_para_5[kk];
		iAngle++;
	}
	mole_last_angle_idx[iMole] = iAngle;

	// dihedral
	mole_first_dih_idx[iMole] = iDih;
	for (kk=sample_mole_first_dih_idx[iSpecieSelected]; kk
			<sample_mole_last_dih_idx[iSpecieSelected]; kk++)
	{
		dih_idx[iDih][0] = sample_dih_idx[kk][0] + idMoleInSpecie
			*sample_natom_per_mole[iSpecieSelected];
		dih_idx[iDih][1] = sample_dih_idx[kk][1] + idMoleInSpecie
			*sample_natom_per_mole[iSpecieSelected];
		dih_idx[iDih][2] = sample_dih_idx[kk][2] + idMoleInSpecie
			*sample_natom_per_mole[iSpecieSelected];
		dih_idx[iDih][3] = sample_dih_idx[kk][3] + idMoleInSpecie
			*sample_natom_per_mole[iSpecieSelected];
		dih_type[iDih] = sample_dih_type[kk];
		c1[iDih] = sample_c1[kk];
		c2[iDih] = sample_c2[kk];
		c3[iDih] = sample_c3[kk];
		c4[iDih] = sample_c4[kk];
		iDih++;
	}
	mole_last_dih_idx[iMole] = iDih;

	// improper
	mole_first_imp_idx[iMole] = iImp;
	for (kk=sample_mole_first_imp_idx[iSpecieSelected]; kk
			<sample_mole_last_imp_idx[iSpecieSelected]; kk++)
	{
		imp_idx[iImp][0] = sample_imp_idx[kk][0] + idMoleInSpecie
			*sample_natom_per_mole[iSpecieSelected];
		imp_idx[iImp][1] = sample_imp_idx[kk][1] + idMoleInSpecie
			*sample_natom_per_mole[iSpecieSelected];
		imp_idx[iImp][2] = sample_imp_idx[kk][2] + idMoleInSpecie
			*sample_natom_per_mole[iSpecieSelected];
		imp_idx[iImp][3] = sample_imp_idx[kk][3] + idMoleInSpecie
			*sample_natom_per_mole[iSpecieSelected];
		imp_type[iImp] = sample_imp_type[kk];
		komega[iImp] = sample_komega[kk];
		omega0[iImp] = sample_omega0[kk];
		iImp++;
	}
	mole_last_imp_idx[iMole] = iImp;

	// non-bonded
	mole_first_nbp_idx[iMole] = iNbp;
	for (kk=sample_mole_first_nbp_idx[iSpecieSelected]; kk
			<sample_mole_last_nbp_idx[iSpecieSelected]; kk++)
	{
		nbp_idx[iNbp][0] = sample_nbp_idx[kk][0] + idMoleInSpecie
			*sample_natom_per_mole[iSpecieSelected];
		nbp_idx[iNbp][1] = sample_nbp_idx[kk][1] + idMoleInSpecie
			*sample_natom_per_mole[iSpecieSelected];
		iNbp++;
	}
	mole_last_nbp_idx[iMole] = iNbp;

	specie_last_atom_idx[iSpecieSelected] = iAtom;
	specie_last_mole_idx[iSpecieSelected] = iMole;
	
	// track the uninit lists
	idAtomUninit = iAtom;
	idBondUninit = iBond;
	idAngleUninit = iAngle;
	idDihUninit = iDih;
	idImpUninit = iImp;
	idNbpUninit = iNbp;

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

		// use the Sample of this specie to build up the insert molecule
		BulidInsertedMole(iSpecieSelected, iMoleSelected);

		// Initialize bond, angle, dih, imp, nbp if the molecule status is un-initialized
		if (mole_status[iMoleSelected] == MOLE_STATUS_UNINIT)
		{
			InitInsertedMole(iSpecieSelected, iMoleSelected);
		}

		// set the status of the molecule to normal
		mole_status[iMoleSelected] = MOLE_STATUS_NORMAL;

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
