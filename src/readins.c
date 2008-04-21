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
	char buffer[LONG_STRING_LENGTH];

	fprintf(stderr,"Reading input file...\n");
	fprintf(fpouts, "Reading input file...\n");
	/* read input file */
	fpins = fopen(INPUT,"r");
	fgets(title, LONG_STRING_LENGTH, fpins);
	sscanf(fgets(buffer, LONG_STRING_LENGTH, fpins), "%d", &ij);
	sscanf(fgets(buffer, LONG_STRING_LENGTH, fpins), "%d", &jk);
	sscanf(fgets(buffer, LONG_STRING_LENGTH, fpins), "%lf", &sigma_base);
	sscanf(fgets(buffer, LONG_STRING_LENGTH, fpins), "%lf", &epsilon_base);
	sscanf(fgets(buffer, LONG_STRING_LENGTH, fpins), "%lf", &mass_base);
	// Calculate reduced unit for time step, unit is fs.
	time_base = sqrt(mass_base*1.0e-27/(epsilon_base*BOLTZMAN_CONSTANT))*sigma_base*1.0e5;
	// Calculate reduced unit for pressure, unit is bar.
	pressure_base = BOLTZMAN_CONSTANT*epsilon_base*1.0e25/(sigma_base*sigma_base*sigma_base);
	sscanf(fgets(buffer, LONG_STRING_LENGTH, fpins), "%lf", &treq);
	treq /= epsilon_base;
	sscanf(fgets(buffer, LONG_STRING_LENGTH, fpins), "%lf", &preq);
	preq /= pressure_base;
	sscanf(fgets(buffer, LONG_STRING_LENGTH, fpins), "%lf", &boxlx);
	sscanf(fgets(buffer, LONG_STRING_LENGTH, fpins), "%lf", &boxly);
	sscanf(fgets(buffer, LONG_STRING_LENGTH, fpins), "%lf", &boxlz);
	boxlx /= sigma_base;
	boxly /= sigma_base;
	boxlz /= sigma_base;
	sscanf(fgets(buffer, LONG_STRING_LENGTH, fpins), "%lf", &rcuton);
	sscanf(fgets(buffer, LONG_STRING_LENGTH, fpins), "%lf", &rcutoff);
	sscanf(fgets(buffer, LONG_STRING_LENGTH, fpins), "%lf", &rcutoffelec);
	rcuton /= sigma_base;
	rcutoff /= sigma_base;
	rcutoffelec /= sigma_base;
	sscanf(fgets(buffer, LONG_STRING_LENGTH, fpins), "%d", &isLJlrcOn);
	sscanf(fgets(buffer, LONG_STRING_LENGTH, fpins), "%d", &isLJswitchOn);
	sscanf(fgets(buffer, LONG_STRING_LENGTH, fpins), "%d", &iStart_option);
	sscanf(fgets(buffer, LONG_STRING_LENGTH, fpins), "%d", &nstep);
	sscanf(fgets(buffer, LONG_STRING_LENGTH, fpins), "%d", &nstep_eq);
	// Set the equilibrium run flag to true if number of equilibrium run is greater than 0
	bEquilibrium = (nstep_eq>0) ? true : false;
	sscanf(fgets(buffer, LONG_STRING_LENGTH, fpins), "%d", &nstep_ave);
	sscanf(fgets(buffer, LONG_STRING_LENGTH, fpins), "%d", &nstep_print);
	sscanf(fgets(buffer, LONG_STRING_LENGTH, fpins), "%d", &nstep_save);
	sscanf(fgets(buffer, LONG_STRING_LENGTH, fpins), "%d", &nstep_ss);
	sscanf(fgets(buffer, LONG_STRING_LENGTH, fpins), "%d", &nstep_trj);
	sscanf(fgets(buffer, LONG_STRING_LENGTH, fpins), "%lf", &delt);
	delt /= time_base;
	sscanf(fgets(buffer, LONG_STRING_LENGTH, fpins), "%d", &nstep_inner);
	sscanf(fgets(buffer, LONG_STRING_LENGTH, fpins), "%lf", &f0);
	sscanf(fgets(buffer, LONG_STRING_LENGTH, fpins), "%d", &what_simulation);
	sscanf(fgets(buffer, LONG_STRING_LENGTH, fpins), "%d", &what_ensemble);
	sscanf(fgets(buffer, LONG_STRING_LENGTH, fpins), "%d", &nconstraint);
	sscanf(fgets(buffer, LONG_STRING_LENGTH, fpins), "%d", &iInterMolePotType);
	sscanf(fgets(buffer, LONG_STRING_LENGTH, fpins), "%d", &iChargeType);
	sscanf(fgets(buffer, LONG_STRING_LENGTH, fpins), "%d", &iSF_type);
	sscanf(fgets(buffer, LONG_STRING_LENGTH, fpins), "%d", &iExternal_FF_type);
	fclose(fpins);
	
	PSAMPLE_MOLECULE pSampleMole;
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
	sscanf(fgets(buffer, LONG_STRING_LENGTH, fpcfg), "%d", &nspecie);
	assert(nspecie<NSPECIE_MAX);
	for (ii=0; ii<nspecie; ii++)
	{
	    pSampleMole = sample_mole + ii;
		fscanf(fpcfg, "%s\n", pSampleMole->mole_name);
		sscanf(fgets(buffer, LONG_STRING_LENGTH, fpcfg), "%d", &nmole_per_specie[ii]);
		
		// Read atom information
		sscanf(fgets(buffer, LONG_STRING_LENGTH, fpcfg), "%d", &pSampleMole->natom);
		natom_per_specie[ii] = pSampleMole->natom*nmole_per_specie[ii];
		natom += natom_per_specie[ii]; // Total number of atoms
		nmole += nmole_per_specie[ii]; // Total number of molecules
		assert(natom<NATOM_MAX);
		assert(nmole<NMOLE_MAX);
		// Set the weight of the molecule
		pSampleMole->mw = 0.0;
		// Read detailed atom information
		// Input parameters have units of Angstrom, K, Kg/mol
		for (jj=0;jj<pSampleMole->natom;jj++)
		{
			fscanf(fpcfg,
					"%d %s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d\n",
					&pSampleMole->type[jj], pSampleMole->atom_name[jj], &pSampleMole->xx[jj],
					&pSampleMole->yy[jj], &pSampleMole->zz[jj], &pSampleMole->ee[jj],
					&pSampleMole->ff[jj], &pSampleMole->gg[jj], &pSampleMole->aw[jj],
					&pSampleMole->epsilon[jj], &pSampleMole->sigma[jj], &pSampleMole->charge[jj],
					&pSampleMole->ghost_type[jj], &pSampleMole->tasos_type[jj]);
			// Reduce units
			pSampleMole->xx[jj] /=sigma_base;
			pSampleMole->yy[jj] /=sigma_base;
			pSampleMole->zz[jj] /=sigma_base;
			pSampleMole->ee[jj] /=sigma_base;
			pSampleMole->ff[jj] /=sigma_base;
			pSampleMole->gg[jj] /=sigma_base;
			pSampleMole->sigma[jj] /= sigma_base;
			pSampleMole->epsilon[jj] /= epsilon_base;
			// Calculate the reduced mass
			pSampleMole->aw[jj] = (pSampleMole->aw[jj]/AVOGADRO)/(mass_base*1.0e-27);
			pSampleMole->mw += pSampleMole->aw[jj];
		}
		
		// Read bond information
		sscanf(fgets(buffer, LONG_STRING_LENGTH, fpcfg), "%d", &pSampleMole->nbond);
		nbond += pSampleMole->nbond*nmole_per_specie[ii]; // Total number of bonds
		assert(nbond<=NBOND_MAX);
		// Read detailed bond information
		// Input parameters have units of Angstrom, K
		for (jj=0;jj<pSampleMole->nbond;jj++)
		{
			fscanf(fpcfg, "%d %d %d %lf %lf %lf\n", &pSampleMole->bnd_idx[jj][0],
					&pSampleMole->bnd_idx[jj][1], &pSampleMole->bnd_type[jj],
					&pSampleMole->Kb[jj], &pSampleMole->Req[jj], &pSampleMole->alpha[jj]);
			if (pSampleMole->bnd_type[jj]==BOND_MORSE)
			{ 
				pSampleMole->Kb[jj] /= epsilon_base;
			}
			else
			{
				pSampleMole->Kb[jj] = pSampleMole->Kb[jj]*sigma_base*sigma_base/epsilon_base;
			}
			pSampleMole->Req[jj] /= sigma_base;
			pSampleMole->alpha[jj] /= sigma_base;
		}

		// Read angle information
		sscanf(fgets(buffer, LONG_STRING_LENGTH, fpcfg), "%d", &pSampleMole->nangle);
		nangle += pSampleMole->nangle*nmole_per_specie[ii]; // Total number of angles
		assert(nangle<=NANGLE_MAX);
		// Read detailed angle information
		for (jj=0; jj<pSampleMole->nangle; jj++)
		{
			fscanf(fpcfg, "%d %d %d %d %lf %lf %lf %lf %lf\n",
					&pSampleMole->agl_idx[jj][0], &pSampleMole->agl_idx[jj][1],
					&pSampleMole->agl_idx[jj][2], &pSampleMole->agl_type[jj],
					&pSampleMole->Ktheta[jj], &pSampleMole->Thetaeq[jj],
					&pSampleMole->agl_para_3[jj], &pSampleMole->agl_para_4[jj],
					&pSampleMole->agl_para_5[jj]); // these 3 parameters only for TRwater
			if (pSampleMole->agl_type[jj]==ANGLE_TR_WATER)
			{
				pSampleMole->Ktheta[jj] = pSampleMole->Ktheta[jj]*sigma_base*sigma_base/epsilon_base;
				pSampleMole->Thetaeq[jj] = pSampleMole->Thetaeq[jj]*sigma_base*sigma_base/epsilon_base;
				pSampleMole->agl_para_3[jj] = pSampleMole->agl_para_3[jj]*sigma_base*sigma_base/epsilon_base;
				pSampleMole->agl_para_4[jj] /= sigma_base;
				pSampleMole->agl_para_5[jj] /= sigma_base;
			}
			else
			{
				pSampleMole->Ktheta[jj] /= epsilon_base;
			}
		}

		// Read in dihedral information
		sscanf(fgets(buffer, LONG_STRING_LENGTH, fpcfg), "%d", &pSampleMole->ndih);
		ndih += pSampleMole->ndih*nmole_per_specie[ii]; // Total number of dihedrals
		assert(ndih<=NDIH_MAX);
		// Readin detailed dihedral information
		for (jj=0; jj<pSampleMole->ndih; jj++)
		{
			fscanf(fpcfg, "%d %d %d %d %d %lf %lf %lf %lf\n",
					&pSampleMole->dih_idx[jj][0], &pSampleMole->dih_idx[jj][1],
					&pSampleMole->dih_idx[jj][2], &pSampleMole->dih_idx[jj][3],
					&pSampleMole->dih_type[jj], &pSampleMole->c1[jj], &pSampleMole->c2[jj],
					&pSampleMole->c3[jj], &pSampleMole->c4[jj]);
			if (pSampleMole->dih_type[jj]==DIH_OPLS_COSIN)
			{
				pSampleMole->c1[jj] /= epsilon_base;
				pSampleMole->c2[jj] /= epsilon_base;
				pSampleMole->c3[jj] /= epsilon_base;
				pSampleMole->c4[jj] /= epsilon_base;
			}
			else if (pSampleMole->dih_type[jj]==DIH_CHARMM)
			{
				pSampleMole->c1[jj] /= epsilon_base;
			}
		}

		// Read in improper information
		sscanf(fgets(buffer, LONG_STRING_LENGTH, fpcfg), "%d", &pSampleMole->nimp);
		nimp += pSampleMole->nimp*nmole_per_specie[ii]; // Total number of impropers
		assert(nimp<=NIMP_MAX);
		// Read detailed improper dihedral information
		for (jj=0; jj<pSampleMole->nimp; jj++)
		{
			fscanf(fpcfg, "%d %d %d %d %d %lf %lf\n", &pSampleMole->imp_idx[jj][0],
					&pSampleMole->imp_idx[jj][1], &pSampleMole->imp_idx[jj][2],
					&pSampleMole->imp_idx[jj][3], &pSampleMole->imp_type[jj],
					&pSampleMole->komega[jj], &pSampleMole->omega0[jj]);
			if (pSampleMole->imp_type[jj]==IMP_CHARMM)
			{
				pSampleMole->komega[jj] /= epsilon_base;
			}
		}

		// Read in non-bonded pair list
		sscanf(fgets(buffer, LONG_STRING_LENGTH, fpcfg), "%d", &pSampleMole->nnbp);
		nnbp += pSampleMole->nnbp*nmole_per_specie[ii]; // Total number of non-bonded pairs
		assert(nnbp<=NNBP_MAX);
		// Read detailed nonbonded pair information
		for (jj=0; jj<pSampleMole->nnbp; jj++)
		{
			fscanf(fpcfg, "%d %d\n", &pSampleMole->nbp_idx[jj][0],
					&pSampleMole->nbp_idx[jj][1]);
		}
	}
	fclose(fpcfg);

	// Read in coordinates
	fprintf(stderr,"Reading initial coordinates of the system...\n");
	fprintf(fpouts, "Reading initial coordinates of the system...\n");
	fpcoords = fopen(COORDSIN,"r");
	fgets(buffer, LONG_STRING_LENGTH, fpcoords);
	fgets(buffer, LONG_STRING_LENGTH, fpcoords);
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
