/**
 * Project: mspms2
 * File: readins.c
 * 
 * Copyright (C) 2008    Yang Wang <ywangd@gmail.com>
 * Created @ 15/04/2008
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
#include <time.h>
#include "mspms2.h"

/* Read in input and config files */
int readins()
{
	int ii, jj;
	char buffer[STRING_LENGTH];
	int atomid;
	int iAtom, iBond, iAngle, iDih, iImp, iNbp;

	fprintf(stderr,"Reading input file...\n");
	fprintf(fpouts, "Reading input file...\n");
	/* read input file */
	fpins = fopen(INPUT,"r");
	fgets(title, STRING_LENGTH, fpins);
	sscanf(fgets(buffer, STRING_LENGTH, fpins), "%d", &ij);
	sscanf(fgets(buffer, STRING_LENGTH, fpins), "%d", &jk);
	sscanf(fgets(buffer, STRING_LENGTH, fpins), "%lf", &sigma_base);
	sscanf(fgets(buffer, STRING_LENGTH, fpins), "%lf", &epsilon_base);
	sscanf(fgets(buffer, STRING_LENGTH, fpins), "%lf", &mass_base);
	// Calculate reduced unit for time step, unit is fs.
	time_base = sqrt(mass_base*1.0e-27/(epsilon_base*BOLTZMAN_CONSTANT))*sigma_base*1.0e5;
	// Calculate reduced unit for pressure, unit is bar.
	pressure_base = BOLTZMAN_CONSTANT*epsilon_base*1.0e25/(sigma_base*sigma_base*sigma_base);
	sscanf(fgets(buffer, STRING_LENGTH, fpins), "%lf", &treq);
	treq /= epsilon_base;
	sscanf(fgets(buffer, STRING_LENGTH, fpins), "%lf", &preq);
	preq /= pressure_base;
	sscanf(fgets(buffer, STRING_LENGTH, fpins), "%lf", &boxlx);
	sscanf(fgets(buffer, STRING_LENGTH, fpins), "%lf", &boxly);
	sscanf(fgets(buffer, STRING_LENGTH, fpins), "%lf", &boxlz);
	boxlx /= sigma_base;
	boxly /= sigma_base;
	boxlz /= sigma_base;
	sscanf(fgets(buffer, STRING_LENGTH, fpins), "%lf", &rcuton);
	sscanf(fgets(buffer, STRING_LENGTH, fpins), "%lf", &rcutoff);
	sscanf(fgets(buffer, STRING_LENGTH, fpins), "%lf", &rcutoffelec);
	rcuton /= sigma_base;
	rcutoff /= sigma_base;
	rcutoffelec /= sigma_base;
	sscanf(fgets(buffer, STRING_LENGTH, fpins), "%d", &isLJlrcOn);
	sscanf(fgets(buffer, STRING_LENGTH, fpins), "%d", &isLJswitchOn);
	sscanf(fgets(buffer, STRING_LENGTH, fpins), "%d", &iStart_option);
	sscanf(fgets(buffer, STRING_LENGTH, fpins), "%d", &nstep);
	sscanf(fgets(buffer, STRING_LENGTH, fpins), "%d", &nstep_eq);
	// Set the equilibrium run flag to true if number of equilibrium run is greater than 0
	bEquilibrium = (nstep_eq>0) ? true : false;
	sscanf(fgets(buffer, STRING_LENGTH, fpins), "%d", &nstep_ave);
	sscanf(fgets(buffer, STRING_LENGTH, fpins), "%d", &nstep_print);
	sscanf(fgets(buffer, STRING_LENGTH, fpins), "%d", &nstep_save);
	sscanf(fgets(buffer, STRING_LENGTH, fpins), "%d", &nstep_ss);
	sscanf(fgets(buffer, STRING_LENGTH, fpins), "%d", &nstep_trj);
	sscanf(fgets(buffer, STRING_LENGTH, fpins), "%lf", &delt);
	delt /= time_base;
	sscanf(fgets(buffer, STRING_LENGTH, fpins), "%d", &nstep_inner);
	sscanf(fgets(buffer, STRING_LENGTH, fpins), "%lf", &f0);
	sscanf(fgets(buffer, STRING_LENGTH, fpins), "%d", &what_simulation);
	sscanf(fgets(buffer, STRING_LENGTH, fpins), "%d", &what_ensemble);
	sscanf(fgets(buffer, STRING_LENGTH, fpins), "%d", &nconstraint);
	sscanf(fgets(buffer, STRING_LENGTH, fpins), "%d", &iInterMolePotType);
	sscanf(fgets(buffer, STRING_LENGTH, fpins), "%d", &iChargeType);
	sscanf(fgets(buffer, STRING_LENGTH, fpins), "%d", &iSF_type);
	sscanf(fgets(buffer, STRING_LENGTH, fpins), "%d", &iExternal_FF_type);
	fclose(fpins);
	
	iAtom = 0;
	iBond = 0;
	iAngle = 0;
	iDih = 0;
	iImp = 0;
	iNbp = 0;
	natom = 0;
	nmole = 0;
	nbond = 0;
	nangle = 0;
	ndih = 0;
	nimp = 0;
	nnbp = 0;
	
	fprintf(stderr,"Reading cfg file...\n");
	fprintf(fpouts, "Reading cfg file...\n");
	// Read config file 
	fpcfg = fopen(CONFIG,"r");
	fscanf(fpcfg, "%s\n", sysname);
	sscanf(fgets(buffer, STRING_LENGTH, fpcfg), "%d", &nspecie);
	assert(nspecie<NSPECIE_MAX);
	for (ii=0; ii<nspecie; ii++)
	{
		fscanf(fpcfg, "%s\n", szSpecieName[ii]);
		sscanf(fgets(buffer, STRING_LENGTH, fpcfg), "%d", &nmole_per_specie[ii]);
		
		// Read atom information
		sscanf(fgets(buffer, STRING_LENGTH, fpcfg), "%d", &sample_natom_per_mole[ii]);
		natom_per_specie[ii] = sample_natom_per_mole[ii]*nmole_per_specie[ii];
		natom += natom_per_specie[ii]; // Total number of atoms
		nmole += nmole_per_specie[ii]; // Total number of molecules
		assert(natom<NATOM_MAX);
		assert(nmole<NMOLE_MAX);
		sample_mole_first_atom_idx[ii] = iAtom; // First atom of this sample molecule
		iAtom += sample_natom_per_mole[ii]; // Number of atoms for this sample molecule
		sample_mole_last_atom_idx[ii] = iAtom; // Last atom of this sample molecule
		// Set the weight of the molecule
		sample_mw[ii] = 0.0;
		// Read detailed atom information
		// Input parameters have units of Angstrom, K, Kg/mol
		for (jj=sample_mole_first_atom_idx[ii];jj<sample_mole_last_atom_idx[ii];jj++)
		{
			fscanf(fpcfg,
					"%d %s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d\n",
					&atomid, sample_atomname[jj], &sample_xx[jj],
					&sample_yy[jj], &sample_zz[jj], &sample_ee[jj],
					&sample_ff[jj], &sample_gg[jj], &sample_aw[jj],
					&sample_epsilon[jj], &sample_sigma[jj], &sample_charge[jj],
					&sample_ghost_type[jj], &sample_tasos_type[jj]);
			// Reduce units
			sample_xx[jj] /=sigma_base;
			sample_yy[jj] /=sigma_base;
			sample_zz[jj] /=sigma_base;
			sample_ee[jj] /=sigma_base;
			sample_ff[jj] /=sigma_base;
			sample_gg[jj] /=sigma_base;
			sample_sigma[jj] /= sigma_base;
			sample_epsilon[jj] /= epsilon_base;
			// Calculate the reduced mass
			sample_aw[jj] = (sample_aw[jj]/AVOGADRO)/(mass_base*1.0e-27);
			sample_mw[ii] += sample_aw[jj];
		}
		
		// Read bond information
		sscanf(fgets(buffer, STRING_LENGTH, fpcfg), "%d", &sample_nbond_per_mole[ii]);
		nbond += sample_nbond_per_mole[ii]*nmole_per_specie[ii]; // Total number of bonds
		assert(nbond<=NBOND_MAX);
		sample_mole_first_bond_idx[ii] = iBond; // First bond of this sample molecule
		iBond += sample_nbond_per_mole[ii]; // Number of bonds of this sample molecule
		sample_mole_last_bond_idx[ii] = iBond; // Last bond of this sample molecule
		// Read detailed bond information
		// Input parameters have units of Angstrom, K
		for (jj=sample_mole_first_bond_idx[ii];jj<sample_mole_last_bond_idx[ii];jj++)
		{
			fscanf(fpcfg, "%d %d %d %lf %lf %lf\n", &sample_bond_idx[jj][0],
					&sample_bond_idx[jj][1], &sample_bond_type[jj],
					&sample_Kb[jj], &sample_Req[jj], &sample_alpha[jj]);
			if (sample_bond_type[jj]==BOND_MORSE)
			{ 
				sample_Kb[jj] /= epsilon_base;
			}
			else
			{
				sample_Kb[jj] = sample_Kb[jj]*sigma_base*sigma_base/epsilon_base;
			}
			sample_Req[jj] /= sigma_base;
			sample_alpha[jj] /= sigma_base;
		}

		// Read angle information
		sscanf(fgets(buffer, STRING_LENGTH, fpcfg), "%d", &sample_nangle_per_mole[ii]);
		nangle += sample_nangle_per_mole[ii]*nmole_per_specie[ii]; // Total number of angles
		assert(nangle<=NANGLE_MAX);
		sample_mole_first_angle_idx[ii] = iAngle; // First angle of this sample molecule
		iAngle += sample_nangle_per_mole[ii]; // Number of angles of this sample molecule
		sample_mole_last_angle_idx[ii] = iAngle; // Last angle of this sample molecule
		// Read detailed angle information
		for (jj=sample_mole_first_angle_idx[ii]; jj<sample_mole_last_angle_idx[ii]; jj++)
		{
			fscanf(fpcfg, "%d %d %d %d %lf %lf %lf %lf %lf\n",
					&sample_angle_idx[jj][0], &sample_angle_idx[jj][1],
					&sample_angle_idx[jj][2], &sample_angle_type[jj],
					&sample_Ktheta[jj], &sample_Thetaeq[jj],
					&sample_agl_para_3[jj], &sample_agl_para_4[jj],
					&sample_agl_para_5[jj]); // these 3 parameters only for TRwater
			if (sample_angle_type[jj]==ANGLE_TR_WATER)
			{
				sample_Ktheta[jj] = sample_Ktheta[jj]*sigma_base*sigma_base/epsilon_base;
				sample_Thetaeq[jj] = sample_Thetaeq[jj]*sigma_base*sigma_base/epsilon_base;
				sample_agl_para_3[jj] = sample_agl_para_3[jj]*sigma_base*sigma_base/epsilon_base;
				sample_agl_para_4[jj] /= sigma_base;
				sample_agl_para_5[jj] /= sigma_base;
			}
			else
			{
				sample_Ktheta[jj] /= epsilon_base;
			}
		}

		// Read in dihedral information
		sscanf(fgets(buffer, STRING_LENGTH, fpcfg), "%d", &sample_ndih_per_mole[ii]);
		ndih += sample_ndih_per_mole[ii]*nmole_per_specie[ii]; // Total number of dihedrals
		assert(ndih<=NDIH_MAX);
		sample_mole_first_dih_idx[ii] = iDih; // First dihedral of this sample molecule
		iDih += sample_ndih_per_mole[ii]; // Number of dihedrals in this sample molecule
		sample_mole_last_dih_idx[ii] = iDih; // Last dihedral of this sample molecule
		// Readin detailed dihedral information
		for (jj=sample_mole_first_dih_idx[ii]; jj<sample_mole_last_dih_idx[ii]; jj++)
		{
			fscanf(fpcfg, "%d %d %d %d %d %lf %lf %lf %lf\n",
					&sample_dih_idx[jj][0], &sample_dih_idx[jj][1],
					&sample_dih_idx[jj][2], &sample_dih_idx[jj][3],
					&sample_dih_type[jj], &sample_c1[jj], &sample_c2[jj],
					&sample_c3[jj], &sample_c4[jj]);
			if (sample_dih_type[jj]==DIH_OPLS_COSIN)
			{
				sample_c1[jj] /= epsilon_base;
				sample_c2[jj] /= epsilon_base;
				sample_c3[jj] /= epsilon_base;
				sample_c4[jj] /= epsilon_base;
			}
			else if (sample_dih_type[jj]==DIH_CHARMM)
			{
				sample_c1[jj] /= epsilon_base;
			}
		}

		// Read in improper information
		sscanf(fgets(buffer, STRING_LENGTH, fpcfg), "%d", &sample_nimp_per_mole[ii]);
		nimp += sample_nimp_per_mole[ii]*nmole_per_specie[ii]; // Total number of impropers
		assert(nimp<=NIMP_MAX);
		sample_mole_first_imp_idx[ii] = iImp; // First improper of this sample molecule
		iImp += sample_nimp_per_mole[ii]; // Number of impropers of this sample molecule
		sample_mole_last_imp_idx[ii] = iImp; // Last impropers of this sample molecule
		// Read detailed improper dihedral information
		for (jj=sample_mole_first_imp_idx[ii]; jj<sample_mole_last_imp_idx[ii]; jj++)
		{
			fscanf(fpcfg, "%d %d %d %d %d %lf %lf\n", &sample_imp_idx[jj][0],
					&sample_imp_idx[jj][1], &sample_imp_idx[jj][2],
					&sample_imp_idx[jj][3], &sample_imp_type[jj],
					&sample_komega[jj], &sample_omega0[jj]);
			if (sample_imp_type[jj]==IMP_CHARMM)
			{
				sample_komega[jj] /= epsilon_base;
			}
		}

		// Read in non-bonded pair list
		sscanf(fgets(buffer, STRING_LENGTH, fpcfg), "%d", &sample_nnbp_per_mole[ii]);
		nnbp += sample_nnbp_per_mole[ii]*nmole_per_specie[ii]; // Total number of non-bonded pairs
		assert(nnbp<=NNBP_MAX);
		sample_mole_first_nbp_idx[ii] = iNbp; // First nbp of this sample molecule
		iNbp += sample_nnbp_per_mole[ii]; // Number of nbp of this sample molecule
		sample_mole_last_nbp_idx[ii] = iNbp; // Last nbp of this sample molecule
		// Read detailed nonbonded pair information
		for (jj=sample_mole_first_nbp_idx[ii]; jj<sample_mole_last_nbp_idx[ii]; jj++)
		{
			fscanf(fpcfg, "%d %d\n", &sample_nbp_idx[jj][0],
					&sample_nbp_idx[jj][1]);
		}
	}
	fclose(fpcfg);

	// Read in coordinates
	fprintf(stderr,"Reading initial coordinates of the system...\n");
	fprintf(fpouts, "Reading initial coordinates of the system...\n");
	fpcoords = fopen(COORDSIN,"r");
	fgets(buffer, STRING_LENGTH, fpcoords);
	fgets(buffer, STRING_LENGTH, fpcoords);
	for (ii=0; ii<natom; ii++)
	{
		fscanf(fpcoords, "%s %lf %lf %lf\n", buffer, &xx[ii], &yy[ii], &zz[ii]);
		// Reduce units
		xx[ii] /= sigma_base;
		yy[ii] /= sigma_base;
		zz[ii] /= sigma_base;
	}
	fclose(fpcoords);

	return 0;
}
