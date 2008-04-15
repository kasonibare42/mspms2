#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <ctype.h>
#include <time.h>
#include "mspms2.h"

int velinit()
{
	int ii;
	double px, py, pz;
	double stdvtmp, stdv;
	double totalmass;
	double scaling;

	// if at zero temperature, everything should be zero
	if (treq == 0.0)
	{
		for (ii=0; ii<natom; ii++)
		{
			vx[ii] = vy[ii] = vz[ii] = 0.0;
			ukin = tinst = 0.0;
		}
		return 0;
	}

	totalmass = 0.0;
	px = py = pz = 0.0;
	// mv^2=kT, when m is kg/mol, the equation becomes to mv^2=RT
	// p = sqrt(RTm), v = sqrt(RT/m)
	stdvtmp = sqrt(RGAS*treq);

	// Temperature is K. R is J/K/mol. m is kg/mol. v is m/s
	for (ii=0; ii<natom; ii++)
	{
		stdv = stdvtmp/sqrt(aw[ii]);
		vx[ii] = stdv*gaussran();
		vy[ii] = stdv*gaussran();
		vz[ii] = stdv*gaussran();
		px += vx[ii]*aw[ii];
		py += vy[ii]*aw[ii];
		pz += vz[ii]*aw[ii];
		totalmass += aw[ii];
	}
	// zero the momentum
	px /= totalmass;
	py /= totalmass;
	pz /= totalmass;
	for (ii=0; ii<natom; ii++)
	{
		vx[ii] -= px;
		vy[ii] -= py;
		vz[ii] -= pz;
	}
	// rescale velocity for required temperature
	ukin = 0.0;
	for (ii=0; ii<natom; ii++)
	{
		ukin += aw[ii]*(vx[ii]*vx[ii]+vy[ii]*vy[ii]+vz[ii]*vz[ii]);
	}
	ukin = 0.5*ukin;
	tinst = 2.0*ukin/(RGAS*nfree);
	scaling = sqrt(treq/tinst);
	for (ii=0; ii<natom; ii++)
	{
		vx[ii] *= scaling;
		vy[ii] *= scaling;
		vz[ii] *= scaling;
	}
	// recalculate the kinetic energy and instantaneous temperature
	// should be exactly the set tempature
	ukin = 0.0;
	for (ii=0; ii<natom; ii++)
	{
		ukin += aw[ii]*(vx[ii]*vx[ii]+vy[ii]*vy[ii]+vz[ii]*vz[ii]);
	}
	ukin = 0.5*ukin;
	tinst = 2.0*ukin/(RGAS*nfree);

	return 0;
}


/// Read in electrostatic parametes and initialize related variables
int fnInitCharge()
{
	int ii;
	const int datalen = 200;
	char buffer[200];
	char keyword[100];

	fprintf(stderr,"Reading input data for Electrostatic interactions...\n");
	fprintf(fpouts, "Reading input data for Electrostatic interactions...\n");

	// re-open input file to read extra data section
	fpins = fopen(INPUT,"r");

	while (fgets(buffer, datalen, fpins)!=NULL)
	{
		sscanf(buffer, "%s", keyword);
		for (ii=0; ii<strlen(keyword); ii++)
		{
			keyword[ii] = toupper(keyword[ii]);
		}
		if (!strcmp(keyword, "ELECTROSTATIC"))
		{
			fprintf(stderr,"Data section for electrostatics found...\n");
			fprintf(fpouts, "Data section for electrostatics found...\n");

			// preset all electrostatic methods to false
			// then turn on them according to the input
			isEwaldOn = isWolfOn = isSimpleCoulomb = 0;
			if (iChargeType == ELECTROSTATIC_EWALD)
			{
				isEwaldOn = 1;
				sscanf(fgets(buffer, datalen, fpins), "%lf", &kappa);
				sscanf(fgets(buffer, datalen, fpins), "%d %d %d %d", &KMAXX,
						&KMAXY, &KMAXZ, &KSQMAX);
				sscanf(fgets(buffer, datalen, fpins), "%d %d", &fEwald_BC,
						&fEwald_Dim);

				// Initialization
				// set ewald parameters
				kappasq = kappa*kappa;
				Bfactor_ewald = 1.0/(4.0*kappa*kappa);
				Vfactor_ewald = 2.0*pi/(boxlx*boxly*boxlz);
				TWOPI_LX = 2.0*pi/boxlx;
				TWOPI_LY = 2.0*pi/boxly;
				TWOPI_LZ = 2.0*pi/boxlz;
				// 1D ewald constant
				twopi_over_3v = 2.0*pi/3.0/boxlx/boxly/boxlz;
			}
			else if (iChargeType == ELECTROSTATIC_WOLF)
			{
				isWolfOn = 1;
				sscanf(fgets(buffer, datalen, fpins), "%lf", &kappa);

				// set wolf parameters
				wolfvcon1 = -erfc(kappa*rcutoffelec)/rcutoffelec;
				wolfvcon2 = erfc(kappa*rcutoffelec)/rcutoffelecsq + 2.0*kappa
						*exp(-(kappa *rcutoffelec)*(kappa*rcutoffelec))
						/(sqrt(pi) *rcutoffelec);
				wolffcon1 = 2.0*kappa/sqrt(pi);
				wolffcon2 = -wolfvcon2;
			}
			else if (iChargeType == ELECTROSTATIC_SIMPLE_COULOMB)
			{
				isSimpleCoulomb = 1;
			}

			fclose(fpins);
			return 0;
		} // if keyword found
	} // read through lines
	fprintf(stderr,"Error: data for electrostatics not found.\n");
	fprintf(fpouts, "Error: data for electrostatics not found.\n");
	fclose(fpins);
	exit(1);
}

/// Read in Silvera-Goldman parametes and initialize related variables
int fnInitSG()
{
	int ii;
	const int datalen = 200;
	char buffer[200];
	char keyword[100];

	fprintf(stderr,"Reading input data for Silvera-Goldman potential");
	fprintf(fpouts, "Reading input data for Silvera-Goldman potential");

	// re-open input file to read extra data section
	fpins = fopen(INPUT,"r");

	while (fgets(buffer, datalen, fpins)!=NULL)
	{
		sscanf(buffer, "%s", keyword);
		for (ii=0; ii<strlen(keyword); ii++)
		{
			keyword[ii] = toupper(keyword[ii]);
		}
		if (!strcmp(keyword, "SG"))
		{
			fprintf(stderr,"Data section for Slivera-Goldman found...\n");
			fprintf(fpouts, "Data section for Slivera-Goldman found...\n");

			sscanf(fgets(buffer, datalen, fpins), "%lf %lf %lf", &sg_alpha,
					&sg_beta, &sg_gama);
			sscanf(fgets(buffer, datalen, fpins), "%lf %lf %lf %lf", &sg_c6,
					&sg_c8, &sg_c9, &sg_c10);
			sscanf(fgets(buffer, datalen, fpins), "%lf", &sg_rc);

			/*
			 * * Calculate the constant part of spring constant: 
			 * K = m * P * (kB T)^2 / hbar^2
			 * the dimensions are mass/time^2 = energy/length^2
			 * 
			 * where P is the number of beads per chain
			 * 		m is the mass of each bead
			 */
			/* To use the unit of J/mol/m^2 for the spring constant
			 * the formula should be revised as follows
			 *  spring  = mass*natom*(R*T)^2/(hbar*avogadro)^2
			 * So the units are (kg/mol)*((J/mol/K)*K)^2/(J*s/mol)^2 = J/mol/m^2
			 * 
			 * The mass of the system (system_mass) is calculated in init_vars()
			 */

			for (ii=0; ii<nspecie; ii++)
			{
				spring[ii] = sample_mw[ii]*(RGAS*treq)*(RGAS*treq)/HBAR_AVOGADRO_SQ;
			}

			double rc1, rc2, rc3, rc5, rc9, rc10;
			rc1 = 1.0/rcutoff;
			rc2 = rc1*rc1;
			rc3 = rc2*rc1;
			rc5 = rc2*rc2;
			rc9 = rc3*rc3*rc3;
			rc10 = rc5*rc5;

			sgshift=(HARTREE_TO_J_PER_MOL*exp(sg_alpha-sg_beta*rcutoff-sg_gama*rcutoffsq)
			-(sg_c6*rc3*rc3+sg_c8*rc5*rc3+sg_c10*rc10-sg_c9*rc9));

			fclose(fpins);
			return 0;
		} // if keyword found
	} // read through lines
	fprintf(stderr,"Error: data for Slivera-Goldman not found.\n");
	fprintf(fpouts, "Error: data for Slivera-Goldman not found.\n");
	fclose(fpins);
	exit(1);
}

/**
 * \brief Initialize the excluding pair list.
 */
int make_exclude_list()
{
	int ii, jj, kk;
	int nexcllist;
	int iStartAtom, iEndAtom;
	int iStartBond, iEndBond;
	int iStartAngle, iEndAngle;
	int iStartDihedral, iEndDihedral;
	int iStartNbp, iEndNbp;
	int iatom;

	fprintf(fpouts, "making exclude list... ");

	// the exclude list uses the relative index
	nexcllist = 0;
	for (ii=0; ii<nspecie; ii++)
	{
		// get the start and end atom index for the specie sample
		// They are the start and end atom index for the first molecule of the specie
		iStartAtom = sample_mole_first_atom_idx[ii];
		iEndAtom = sample_mole_last_atom_idx[ii];
		iatom = 0;
		for (jj=iStartAtom; jj<iEndAtom; jj++)
		{
			pointexcl_atom[ii][iatom] = nexcllist;
			// self exclusion
			excllist[nexcllist] = jj - jj;
			nexcllist++;
			// exclude bonded atoms
			// get the first bond index of the first molecule of the specie
			iStartBond = mole_first_bond_idx[specie_first_mole_idx[ii]];
			iEndBond = iStartBond + sample_nbond_per_mole[ii];
			for (kk=iStartBond; kk<iEndBond; kk++)
			{
				if (bond_idx[kk][0] == jj)
				{
					excllist[nexcllist] = bond_idx[kk][1] - jj;
					nexcllist++;
				}
				else if (bond_idx[kk][1] == jj)
				{
					excllist[nexcllist] = bond_idx[kk][0] - jj;
					nexcllist++;
				}
				assert(nexcllist<EXCLUDE_LIST_MAX);
			}
			// exclude angled atoms
			// get first angle index for the first molecule of the specie
			iStartAngle = mole_first_angle_idx[specie_first_mole_idx[ii]];
			iEndAngle = iStartAngle + sample_nangle_per_mole[ii];
			for (kk=iStartAngle; kk<iEndAngle; kk++)
			{
				if (angle_idx[kk][0] == jj)
				{
					excllist[nexcllist] = angle_idx[kk][2] - jj;
					nexcllist++;
				}
				else if (angle_idx[kk][2] == jj)
				{
					excllist[nexcllist] = angle_idx[kk][0] - jj;
					nexcllist++;
				}
				assert(nexcllist<EXCLUDE_LIST_MAX);
			}
			// exclude dihedraled atoms
			// get first dihedral index for the first molecule of the specie
			iStartDihedral = mole_first_dih_idx[specie_first_mole_idx[ii]];
			iEndDihedral = iStartDihedral + sample_ndih_per_mole[ii];
			for (kk=iStartDihedral; kk<iEndDihedral; kk++)
			{
				if (dih_idx[kk][0] == jj)
				{
					excllist[nexcllist] = dih_idx[kk][3] - jj;
					nexcllist++;
				}
				else if (dih_idx[kk][3] == jj)
				{
					excllist[nexcllist] = dih_idx[kk][0] - jj;
					nexcllist++;
				}
				assert(nexcllist<EXCLUDE_LIST_MAX);
			}
			// exclude non-bonded pair atoms
			// get first non-bonded pair index for the first molecule of the specie
			iStartNbp = mole_first_nbp_idx[specie_first_mole_idx[ii]];
			iEndNbp = iStartNbp + sample_nnbp_per_mole[ii];
			for (kk=iStartNbp; kk<iEndNbp; kk++)
			{
				if (nbp_idx[kk][0] == jj)
				{
					excllist[nexcllist] = nbp_idx[kk][1] - jj;
					nexcllist++;
				}
				else if (nbp_idx[kk][1] == jj)
				{
					excllist[nexcllist] = nbp_idx[kk][0] - jj;
					nexcllist++;
				}
				assert(nexcllist<EXCLUDE_LIST_MAX);
			}
			iatom++;
		}
		pointexcl_atom[ii][sample_natom_per_mole[ii]] = nexcllist;
	}
	fprintf(fpouts, "%d excluding pairs\n", nexcllist);

	return 0;
}

int InitReplicateSamples()
{
	int ii, jj, kk;
	int iAtom, iMole, iBond, iAngle, iDih, iImp, iNbp;

	// initialize the properties for maximal number of molecules
	for (ii=0; ii<NMOLE_MAX; ii++)
	{
		mole_status[ii] = MOLE_STATUS_UNINIT;
	}

	// init the real atom, bond, angle, dihedral, improper, non-bonded list
	iAtom = 0;
	iMole = 0;
	iBond = 0;
	iAngle = 0;
	iDih = 0;
	iImp = 0;
	iNbp = 0;
	for (ii=0; ii<nspecie; ii++)
	{
		specie_first_atom_idx[ii] = iAtom;
		specie_first_mole_idx[ii] = iMole;
		for (jj=0; jj<nmole_per_specie[ii]; jj++)
		{
			mole2specie[iMole] = ii;
			mole_status[iMole] = MOLE_STATUS_NORMAL;
			iMIDS_to_PMID[ii][jj] = iMole;

			// atom
			mole_first_atom_idx[iMole] = iAtom;
			for (kk=sample_mole_first_atom_idx[ii]; kk
					<sample_mole_last_atom_idx[ii]; kk++)
			{
				atom2mole[iAtom] = iMole;
				strcpy(atomname[iAtom], sample_atomname[kk]);
				aw[iAtom] = sample_aw[kk];
				epsilon[iAtom] = sample_epsilon[kk];
				sigma[iAtom] = sample_sigma[kk];
				charge[iAtom] = sample_charge[kk];
				ghost_type[iAtom] = sample_ghost_type[kk];
				tasostype[iAtom] = sample_tasostype[kk];
				iAtom++;
			}
			mole_last_atom_idx[iMole] = iAtom;

			// bond
			mole_first_bond_idx[iMole] = iBond;
			for (kk=sample_mole_first_bond_idx[ii]; kk
					<sample_mole_last_bond_idx[ii]; kk++)
			{
				bnd2mole[iBond] = iMole;
				bond_idx[iBond][0] = sample_bond_idx[kk][0] + jj
						*sample_natom_per_mole[ii];
				bond_idx[iBond][1] = sample_bond_idx[kk][1] + jj
						*sample_natom_per_mole[ii];
				bond_type[iBond] = sample_bond_type[kk];
				Kb[iBond] = sample_Kb[kk];
				Req[iBond] = sample_Req[kk];
				alpha[iBond] = sample_alpha[kk];
				iBond++;
			}
			mole_last_bond_idx[iMole] = iBond;

			// angle
			mole_first_angle_idx[iMole] = iAngle;
			for (kk=sample_mole_first_angle_idx[ii]; kk
					<sample_mole_last_angle_idx[ii]; kk++)
			{
				agl2mole[iAngle] = iMole;
				angle_idx[iAngle][0] = sample_angle_idx[kk][0] + jj
						*sample_natom_per_mole[ii];
				angle_idx[iAngle][1] = sample_angle_idx[kk][1] + jj
						*sample_natom_per_mole[ii];
				angle_idx[iAngle][2] = sample_angle_idx[kk][2] + jj
						*sample_natom_per_mole[ii];
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
			for (kk=sample_mole_first_dih_idx[ii]; kk
					<sample_mole_last_dih_idx[ii]; kk++)
			{
				dih2mole[iDih] = iMole;
				dih_idx[iDih][0] = sample_dih_idx[kk][0] + jj
						*sample_natom_per_mole[ii];
				dih_idx[iDih][1] = sample_dih_idx[kk][1] + jj
						*sample_natom_per_mole[ii];
				dih_idx[iDih][2] = sample_dih_idx[kk][2] + jj
						*sample_natom_per_mole[ii];
				dih_idx[iDih][3] = sample_dih_idx[kk][3] + jj
						*sample_natom_per_mole[ii];
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
			for (kk=sample_mole_first_imp_idx[ii]; kk
					<sample_mole_last_imp_idx[ii]; kk++)
			{
				imp2mole[iImp] = iMole;
				imp_idx[iImp][0] = sample_imp_idx[kk][0] + jj
						*sample_natom_per_mole[ii];
				imp_idx[iImp][1] = sample_imp_idx[kk][1] + jj
						*sample_natom_per_mole[ii];
				imp_idx[iImp][2] = sample_imp_idx[kk][2] + jj
						*sample_natom_per_mole[ii];
				imp_idx[iImp][3] = sample_imp_idx[kk][3] + jj
						*sample_natom_per_mole[ii];
				imp_type[iImp] = sample_imp_type[kk];
				komega[iImp] = sample_komega[kk];
				omega0[iImp] = sample_omega0[kk];
				iImp++;
			}
			mole_last_imp_idx[iMole] = iImp;

			// non-bonded
			mole_first_nbp_idx[iMole] = iNbp;
			for (kk=sample_mole_first_nbp_idx[ii]; kk
					<sample_mole_last_nbp_idx[ii]; kk++)
			{
				nbp2mole[iNbp] = iMole;
				nbp_idx[iNbp][0] = sample_nbp_idx[kk][0] + jj
						*sample_natom_per_mole[ii];
				nbp_idx[iNbp][1] = sample_nbp_idx[kk][1] + jj
						*sample_natom_per_mole[ii];
				iNbp++;
			}
			mole_last_nbp_idx[iMole] = iNbp;

			// molecule increament
			iMole++;
		}
		specie_last_atom_idx[ii] = iAtom;
		specie_last_mole_idx[ii] = iMole;
	}

	return 0;
}

/**
 * Check if any dihedral, angle, bond share the same ending pairs.
 */
int CheckUniques()
{
	int ii, jj;
	int iAngle, iDih;
	int FirstMoleIDofSpecieii;
	int FirstAngleIDofTheMole;
	int LastAngleIDofTheMole;
	int FirstDihIDofTheMole;
	int LastDihIDofTheMole;

	// check unique for dihedrals
	// this is for possible ring structures where 1,4 atoms can form multiple dihedrals
	// e.g. 1-2-3-4
	//       \5-6/
	// 1234 and 1564
	// The 14 pair should only be calculated once for energy/force
	// Thats the unique check for
	for (ii=0; ii<ndih; ii++)
	{
		isDih_unique[ii] = true;
	}
	for (ii=0; ii<ndih-1; ii++)
	{
		if (isDih_unique[ii])
		{
			for (jj=ii+1; jj<ndih; jj++)
			{
				if (isDih_unique[jj])
				{
					if (dih_idx[ii][0]==dih_idx[jj][0] && dih_idx[ii][3]
							==dih_idx[jj][3])
					{
						isDih_unique[jj] = false;
						fprintf(
								fpouts,
								"Dihedral %d and dihedral %d have the same ending pairs.\n",
								ii, jj);
					}
					else if (dih_idx[ii][0]==dih_idx[jj][3] && dih_idx[ii][3]
							==dih_idx[jj][0])
					{
						isDih_unique[jj] = false;
						fprintf(
								fpouts,
								"Dihedral %d and dihedral %d have the same ending pairs.\n",
								ii, jj);
					}
				}
			}
		}
	}

	// following codes make sure 14 and 13 do not share same ending pairs
	// this is also for ring kind structures
	for (ii=0; ii<ndih; ii++)
	{
		if (isDih_unique[ii])
		{
			for (jj=0; jj<nangle; jj++)
			{
				if (dih_idx[ii][0]==angle_idx[jj][0] && dih_idx[ii][3]
						==angle_idx[jj][2])
				{
					isDih_unique[jj] = false;
					fprintf(
							fpouts,
							"Dihedral %d and angle %d have the same ending pairs.\n",
							ii, jj);
				}
				else if (dih_idx[ii][0]==angle_idx[jj][2] && dih_idx[ii][3]
						==angle_idx[jj][0])
				{
					isDih_unique[jj] = false;
					fprintf(
							fpouts,
							"Dihedral %d and angle %d have the same ending pairs.\n",
							ii, jj);
				}
			}
		}
	}

	// following codes make sure 14 and 12 do not share the same ending pairs
	for (ii=0; ii<ndih; ii++)
	{
		if (isDih_unique[ii])
		{
			for (jj=0; jj<nbond; jj++)
			{
				if (dih_idx[ii][0]==bond_idx[jj][0] && dih_idx[ii][3]
						==bond_idx[jj][1])
				{
					isDih_unique[jj] = false;
					fprintf(
							fpouts,
							"Dihedral %d and bond %d have the same ending pairs.\n",
							ii, jj);
				}
				else if (dih_idx[ii][0]==bond_idx[jj][1] && dih_idx[ii][3]
						==bond_idx[jj][0])
				{
					isDih_unique[jj] = false;
					fprintf(
							fpouts,
							"Dihedral %d and bond %d have the same ending pairs.\n",
							ii, jj);
				}
			}
		}
	}

	// check unique for angles
	// see above comments for dihedrals
	for (ii=0; ii<nangle; ii++)
	{
		isAngle_unique[ii] = true;
	}
	for (ii=0; ii<nangle-1; ii++)
	{
		if (isAngle_unique[ii])
		{
			for (jj=ii+1; jj<nangle; jj++)
			{
				if (isAngle_unique[jj])
				{
					if (angle_idx[ii][0]==angle_idx[jj][0] && angle_idx[ii][2]
							==angle_idx[jj][2])
					{
						isAngle_unique[jj] = false;
						fprintf(
								fpouts,
								"Angle %d and angle %d have the same ending pairs.\n",
								ii, jj);
					}
					else if (angle_idx[ii][0]==angle_idx[jj][2]
							&& angle_idx[ii][2]==angle_idx[jj][0])
					{
						isAngle_unique[jj] = false;
						fprintf(
								fpouts,
								"Angle %d and angle %d have the same ending pairs.\n",
								ii, jj);
					}
				}
			}
		}
	} // end of checking unique angles

	// following codes make sure 13 and 12 do not share the same ending pairs
	for (ii=0; ii<nangle; ii++)
	{
		if (isAngle_unique[ii])
		{
			for (jj=0; jj<nbond; jj++)
			{
				if (angle_idx[ii][0]==bond_idx[jj][0] && angle_idx[ii][2]
						==bond_idx[jj][1])
				{
					isAngle_unique[jj] = false;
					fprintf(
							fpouts,
							"Angle %d and bond %d have the same ending pairs.\n",
							ii, jj);
				}
				else if (angle_idx[ii][0]==bond_idx[jj][1] && angle_idx[ii][2]
						==bond_idx[jj][0])
				{
					isAngle_unique[jj] = false;
					fprintf(
							fpouts,
							"Angle %d and bond %d have the same ending pairs.\n",
							ii, jj);
				}
			}
		}
	}

	// Set the angle, dihedral unique for Samples using the information from real lists.
	// This is the reverse way to what we have done to populate the real lists using Samples.
	// It is because we do not want to make many changes to the existing code of checking uniques for real lists.
	// So we assign the uniques for Samples using the info from the real lists. 
	iAngle = 0;
	iDih = 0;
	for (ii=0; ii<nspecie; ii++)
	{
		// Get ID of the first molecule for this specie
		FirstMoleIDofSpecieii = iMIDS_to_PMID[ii][0]; // molecule's physical id of the first molecule of specie ii

		// unique angles
		// Get the first and last Angle ID for this molecule
		FirstAngleIDofTheMole = mole_first_angle_idx[FirstMoleIDofSpecieii]; // this molecule's first bond index
		LastAngleIDofTheMole = mole_last_angle_idx[FirstMoleIDofSpecieii]; // this molecule's last bond index
		for (jj=FirstAngleIDofTheMole; jj<LastAngleIDofTheMole; jj++)
		{
			sample_isAngle_unique[iAngle] = isAngle_unique[jj];
			iAngle++;
		}

		// unique dihedrals
		// Get the first and last dihedral ID for this molecule
		FirstDihIDofTheMole = mole_first_dih_idx[FirstMoleIDofSpecieii];
		LastDihIDofTheMole = mole_last_dih_idx[FirstMoleIDofSpecieii];
		for (jj=FirstDihIDofTheMole; jj<LastDihIDofTheMole; jj++)
		{
			sample_isDih_unique[iDih] = isDih_unique[jj];
			iDih++;
		}
	}

	return 0;
}

int InitLJlrcCommonTerms()
{
	// Calculate the common terms for LJ lrc
	int ii, jj;
	int mm, nn;
	double uljlrc_term1, uljlrc_term2;
	double pljlrc_term1, pljlrc_term2;
	int atomid_1, atomid_2;
	double sigmaij, epsilonij;
	double temp1, temp2, temp3;
	// set lrc to zero even long range correction is not needed just in case
	uljlrc = 0.0;
	pljlrc = 0.0;
	uljlrc_term1 = (8.0/9.0)*pi*pow(rcutoff, -9.0);
	uljlrc_term2 = -(8.0/3.0)*pi*pow(rcutoff, -3.0);
	pljlrc_term1 = (32.0/9.0)*pi*pow(rcutoff, -9.0);
	pljlrc_term2 = -(16.0/3.0)*pi*pow(rcutoff, -3.0);
	// calculate long range correction terms for single molecules
	for (mm=0; mm<nspecie; mm++)
	{
		for (nn=0; nn<nspecie; nn++)
		{
			uljlrc_term[mm][nn] = 0.0;
			pljlrc_term[mm][nn] = 0.0;
			// loop through the atoms in one molecule of specie mm
			for (ii=0; ii<sample_natom_per_mole[mm]; ii++)
			{
				// use the first molecule of one specie to do the calculation
				atomid_1 = sample_mole_first_atom_idx[mm] + ii;
				// loop through all the atoms in one molecule of specie nn
				for (jj=0; jj<sample_natom_per_mole[nn]; jj++)
				{
					atomid_2 = sample_mole_first_atom_idx[nn] + ii;
					sigmaij = 0.5*(sample_sigma[atomid_1]
							+sample_sigma[atomid_2]);
					epsilonij = sqrt(sample_epsilon[atomid_1]
							*sample_epsilon[atomid_2]);
					temp1 = pow(sigmaij, 9.0)*uljlrc_term1 + pow(sigmaij, 3.0)
							*uljlrc_term2;
					temp2 = epsilonij*pow(sigmaij, 3.0);
					temp3 = pow(sigmaij, 9.0)*pljlrc_term1 + pow(sigmaij, 3.0)
							*pljlrc_term2;
					uljlrc_term[mm][nn] += temp1*temp2;
					pljlrc_term[mm][nn] += temp3*temp2;

				} // through all atom in one molecule of specie 2
			} // through all atom in one molecule of specie 1

		} // loop through specie 2
	} // loop through speice 1

	return 0;
}

/// Initiate variables 
/**
 * Initialize the real atom, bond, angle, dihedral, improper, non-bonded list
 * using Samples.
 * Initialize the readin variables.
 * Readin additional data section according to the simulation and ensemble
 * type.
 * Check the uniqueness for dihedrals and angles.
 * Calculate the long range corrections.
 */
int init_vars()
{
	int ii, jj, kk;
	int iAtom, iMole, iBond, iAngle, iDih, iImp, iNbp;

	fprintf(stderr,"Initializing variables...\n");
	fprintf(fpouts, "Initializing variables...\n");

	/// Set the un-initialized ID for atom, bond, angle, dih, imp, nbp lists
	idAtomUninit = natom;
	idBondUninit = nbond;
	idAngleUninit = nangle;
	idDihUninit = ndih;
	idImpUninit = nimp;
	idNbpUninit = nnbp;

	/// Set the history maximum records
	natom_hist_max = natom;
	nbond_hist_max = nbond;
	nangle_hist_max = nangle;
	ndih_hist_max = ndih;
	nimp_hist_max = nimp;
	nnbp_hist_max = nnbp;

	nmole_hist_max = nmole;

	/// Initialize the real atom, bond, angle, dihedral, improper, nbp lists using Samples
	InitReplicateSamples();

	/// Initiate file variables, LOG, TRAJECTORY\n
	/// Output file is initialized already at the very beginning of the run.
	fplog = fopen(LOG,"w");
	fptrj = fopen(TRAJECTORY,"wb");

	/// Calculate molecule weight for the real list.
	system_mass = 0.0;
	for (ii=0; ii<nmole; ii++)
	{
		mw[ii] = 0.0;
		for (jj=mole_first_atom_idx[ii]; jj<mole_last_atom_idx[ii]; jj++)
		{
			mw[ii] += aw[jj];
		}
		system_mass += mw[ii];
	}
	/// Zero the number of frames in trajectory file.
	nframe = 0;

	/// Set the starting step to 1, will be changed by load it if it is continue run.
	/// istep is used for printit, the first print should be at step zero.
	istep = 0;
	nstep_start = 1;
	nstep_end = (nstep_eq>0)?nstep_eq:nstep;

	// zero the counts and accums
	rezero();
	
	/// set the counter for equilibrium
	/// it will decrease during the run
	counts[11] = nstep_eq;

	/// initialize random number generator
	rmarin(ij, jk);

	/// calculate the degree of freedom
	nfree = 3*natom - nconstraint;

	// cutoff related
	rcutoffsq = rcutoff*rcutoff;
	rcutoffelecsq = rcutoffelec*rcutoffelec;
	rcutonsq = rcuton*rcuton;
	roff2_minus_ron2_cube = (rcutoffsq-rcutonsq)*(rcutoffsq-rcutonsq)
			*(rcutoffsq-rcutonsq);

	// shift energies, use rcutoff
	shift1 = pow((1/rcutoff), 6.0);
	shift4 = 4*shift1;

	// volume calculation
	boxv = boxlx*boxly*boxlz;

	// delt related
	deltby2 = delt/2.0;
	delts = delt/nstep_inner;
	deltsby2 = delts/2.0;

	dt_outer2 = deltby2;
	dt_outer4 = delt/4.0;
	dt_outer8 = delt/8.0;

	if (what_simulation == MOLECULAR_DYNAMICS)
	{
	}
	else if (what_simulation == HYBRID_MONTE_CARLO) // initialize HMC input data
	{
		init_hmc();
	}
	else if (what_simulation == SIMULATED_ANNEALING)
	{
		init_siman();
	}

	// initialize thermostat/baron stat input data
	if (what_ensemble == NPT)
	{
		init_npt_respa();
	}
	else if (what_ensemble == NVT)
	{
		init_nvt();
	}

	// initialize velocities for needed simulations
	if (what_simulation==MOLECULAR_DYNAMICS || what_simulation==HYBRID_MONTE_CARLO || what_simulation
			==SIMULATED_ANNEALING)
	{
		fprintf(stderr, "initializing velocities...\n");
		fprintf(fpouts, "initializing velocities...\n");
		velinit();
	}

	// If Solid-fluid interaction is required,
	// initiliaze the related variables
	if (sf_type==SF_NANOTUBE_HYPERGEO)
	{
		init_sf_hypergeo();
	}
	else if (sf_type==SF_NANOTUBE_ATOM_EXPLICIT)
	{
		init_sf_atom_explicit();
	}
	else if (sf_type==SF_NANOTUBE_TASOS)
	{
		init_tasos_grid();
	}
	else if (sf_type==SF_NANOTUBE_MY_INTERP)
	{
		init_my_interp();
	}

	// Read in electrostatic parametes and initialize if needed
	if (iChargeType != ELECTROSTATIC_NONE)
	{
		fnInitCharge();
	}

	// Read in SG potential parameters if required
	if (iInterMolePotType == INTER_MOLE_SG)
	{
		fnInitSG();
	}
	// initilize spring constant no matter what, since S.G. potential is not
	// necessarily part of the ring polymer approach, even L.J can be used
	// as inter-polymer bead interactions. But spring constant is always needed.
	for (ii=0; ii<nspecie; ii++)
	{
		spring[ii] = sample_mw[ii]*(RGAS*treq)*(RGAS*treq)/HBAR_AVOGADRO_SQ;
	}


	// Check if any bond, angle, dihedral share the same ending pairs
	CheckUniques();

	// make the exlcuding pair list
	make_exclude_list();

	// Initialize the common terms and calculate LJ lrc if it is requested
	if (isLJlrcOn)
	{
		fprintf(stderr,"calculating LJ long range corrections...\n");
		fprintf(fpouts, "calculating LJ long range corrections...\n");
		// Initialize common terms
		InitLJlrcCommonTerms();
		// calculate the total lj lrc
		calculate_ljlrc();
	}
	
	// if not new run, load from old file
	if (fStart_option!=NEW)
	{
		loadit();
	}

	return 0;
}
