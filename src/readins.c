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
	const int datalen = 200;
	char buffer[200];
	int atomid, moleid;
	int isFirstAtom;
	int last_mole_id;
	int ispecie;
	int iAtom, iMole;
	int iBond, iAngle, iDih, iImp, iNbp;

	fprintf(stderr,"Reading input file...\n");
	fprintf(fpouts, "Reading input file...\n");
	/* read input file */
	fpins = fopen(INPUT,"r");
	sscanf(fgets(buffer, datalen, fpins), "%d %d", &ij, &jk);
	sscanf(fgets(buffer, datalen, fpins), "%lf", &treq);
	sscanf(fgets(buffer, datalen, fpins), "%lf", &preq);
	sscanf(fgets(buffer, datalen, fpins), "%lf %lf %lf", &boxlx, &boxly, &boxlz);
	sscanf(fgets(buffer, datalen, fpins), "%lf %lf %lf", &rcuton, &rcutoff,
			&rcutoffelec);
	// sscanf(fgets(buffer,datalen,fpins), "%s", coords_file);
	sscanf(fgets(buffer, datalen, fpins), "%d %d %d", &nstep, &fStart_option,
			&nstep_eq);
	// Set the equilibrium run flag to true if number of equilibrium run is greater than 0
	bEquilibrium = (nstep_eq>0)?true:false;
	sscanf(fgets(buffer, datalen, fpins), "%d %d %d %d %d", &nstep_ave,
			&nstep_print, &nstep_save, &nstep_ss, &nstep_trj);
	sscanf(fgets(buffer, datalen, fpins), "%lf %d", &delt, &nstep_inner);
	sscanf(fgets(buffer, datalen, fpins), "%lf", &f0);
	sscanf(fgets(buffer, datalen, fpins), "%d %d", &what_simulation,
			&what_ensemble);
	sscanf(fgets(buffer, datalen, fpins), "%d", &iInterMolePotType);
	sscanf(fgets(buffer, datalen, fpins), "%d %d", &isLJlrcOn, &isLJswitchOn);
	sscanf(fgets(buffer, datalen, fpins), "%d", &iChargeType);
	sscanf(fgets(buffer, datalen, fpins), "%d", &nconstraint);
	sscanf(fgets(buffer, datalen, fpins), "%d", &sf_type);
	sscanf(fgets(buffer, datalen, fpins), "%d", &fOtherFF);
	fclose(fpins);

	fprintf(stderr,"reading cfg file...\n");
	fprintf(fpouts, "reading cfg file...\n");
	/* read config file */
	fpcfg = fopen(CONFIG,"r");
	fscanf(fpcfg, "%s\n", sysname);
	sscanf(fgets(buffer, datalen, fpcfg), "%d", &nspecie);
	assert(nspecie<=NSPECIE_MAX);
	iAtom = 0;
	iMole = 0;
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
	for (ii=0; ii<nspecie; ii++)
	{
		fscanf(fpcfg, "%s\n", szSpecieName[ii]);
		sscanf(fgets(buffer, datalen, fpcfg), "%d", &nmole_per_specie[ii]);
		// readin atom information
		sscanf(fgets(buffer, datalen, fpcfg), "%d", &sample_natom_per_mole[ii]);
		natom_per_specie[ii] = sample_natom_per_mole[ii]*nmole_per_specie[ii];
		sample_mole_first_atom_idx[ii] = iAtom;
		natom += sample_natom_per_mole[ii]*nmole_per_specie[ii];
		assert(natom <= NATOM_MAX);
		nmole += nmole_per_specie[ii];
		assert(nmole <= NMOLE_MAX);
		iAtom += sample_natom_per_mole[ii];
		sample_mole_last_atom_idx[ii] = iAtom;
		sample_mw[ii] = 0.0;
		// readin detailed atom information
		for (jj=sample_mole_first_atom_idx[ii]; jj
				<sample_mole_last_atom_idx[ii]; jj++)
		{
			fscanf(fpcfg,
					"%d %s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d\n",
					&atomid, sample_atomname[jj], &sample_xx[jj],
					&sample_yy[jj], &sample_zz[jj], &sample_ee[jj],
					&sample_ff[jj], &sample_gg[jj], &sample_aw[jj],
					&sample_epsilon[jj], &sample_sigma[jj], &sample_charge[jj],
					&sample_ghost_type[jj], &sample_tasostype[jj]);
			sample_mw[ii] += sample_aw[jj];
		}

		// readin bond information
		sscanf(fgets(buffer, datalen, fpcfg), "%d", &sample_nbond_per_mole[ii]);
		sample_mole_first_bond_idx[ii] = iBond;
		nbond += sample_nbond_per_mole[ii]*nmole_per_specie[ii];
		assert(nbond<=NBOND_MAX);
		iBond += sample_nbond_per_mole[ii];
		sample_mole_last_bond_idx[ii] = iBond;
		// readin detailed bond information
		for (jj=sample_mole_first_bond_idx[ii]; jj
				<sample_mole_last_bond_idx[ii]; jj++)
		{
			fscanf(fpcfg, "%d %d %d %lf %lf %lf\n", &sample_bond_idx[jj][0],
					&sample_bond_idx[jj][1], &sample_bond_type[jj],
					&sample_Kb[jj], &sample_Req[jj], &sample_alpha[jj]);
		}

		// read in angle information
		sscanf(fgets(buffer, datalen, fpcfg), "%d", &sample_nangle_per_mole[ii]);
		sample_mole_first_angle_idx[ii] = iAngle;
		nangle += sample_nangle_per_mole[ii]*nmole_per_specie[ii];
		assert(nangle<=NANGLE_MAX);
		iAngle += sample_nangle_per_mole[ii];
		sample_mole_last_angle_idx[ii] = iAngle;
		// read in detailed angle information
		for (jj=sample_mole_first_angle_idx[ii]; jj
				<sample_mole_last_angle_idx[ii]; jj++)
		{
			fscanf(fpcfg, "%d %d %d %d %lf %lf %lf %lf %lf\n",
					&sample_angle_idx[jj][0], &sample_angle_idx[jj][1],
					&sample_angle_idx[jj][2], &sample_angle_type[jj],
					&sample_Ktheta[jj], &sample_Thetaeq[jj],
					&sample_agl_para_3[jj], &sample_agl_para_4[jj],
					&sample_agl_para_5[jj]); // these 3 parameters only for TRwater
		}

		// read in dihedral information
		sscanf(fgets(buffer, datalen, fpcfg), "%d", &sample_ndih_per_mole[ii]);
		sample_mole_first_dih_idx[ii] = iDih;
		ndih += sample_ndih_per_mole[ii]*nmole_per_specie[ii];
		assert(ndih<=NDIH_MAX);
		iDih += sample_ndih_per_mole[ii];
		sample_mole_last_dih_idx[ii] = iDih;
		// readin detailed dihedral information
		for (jj=sample_mole_first_dih_idx[ii]; jj<sample_mole_last_dih_idx[ii]; jj++)
		{
			fscanf(fpcfg, "%d %d %d %d %d %lf %lf %lf %lf\n",
					&sample_dih_idx[jj][0], &sample_dih_idx[jj][1],
					&sample_dih_idx[jj][2], &sample_dih_idx[jj][3],
					&sample_dih_type[jj], &sample_c1[jj], &sample_c2[jj],
					&sample_c3[jj], &sample_c4[jj]);
		}

		// read in improper information
		sscanf(fgets(buffer, datalen, fpcfg), "%d", &sample_nimp_per_mole[ii]);
		sample_mole_first_imp_idx[ii] = iImp;
		nimp += sample_nimp_per_mole[ii]*nmole_per_specie[ii];
		assert(nimp<=NIMP_MAX);
		iImp += sample_nimp_per_mole[ii];
		sample_mole_last_imp_idx[ii] = iImp;
		// readin detailed improper dihedral information
		for (jj=sample_mole_first_imp_idx[ii]; jj<sample_mole_last_imp_idx[ii]; jj++)
		{
			fscanf(fpcfg, "%d %d %d %d %d %lf %lf\n", &sample_imp_idx[jj][0],
					&sample_imp_idx[jj][1], &sample_imp_idx[jj][2],
					&sample_imp_idx[jj][3], &sample_imp_type[jj],
					&sample_komega[jj], &sample_omega0[jj]);
		}

		// read in nonbonded pair list
		sscanf(fgets(buffer, datalen, fpcfg), "%d", &sample_nnbp_per_mole[ii]);
		sample_mole_first_nbp_idx[ii] = iNbp;
		nnbp += sample_nnbp_per_mole[ii]*nmole_per_specie[ii];
		assert(nnbp<=NNBP_MAX);
		iNbp += sample_nnbp_per_mole[ii];
		sample_mole_last_nbp_idx[ii] = iNbp;
		// readin detailed nonbonded pair information
		for (jj=sample_mole_first_nbp_idx[ii]; jj<sample_mole_last_nbp_idx[ii]; jj++)
		{
			fscanf(fpcfg, "%d %d\n", &sample_nbp_idx[jj][0],
					&sample_nbp_idx[jj][1]);
		}
	}
	fclose(fpcfg);

	// read in coordinates
	fprintf(stderr,"reading initial coordinates of the system...\n");
	fprintf(fpouts, "reading initial coordinates of the system...\n");
	fpcoords = fopen(COORDSIN,"r");
	// fscanf(fpcoords, "%[^\n]", buffer);
	fgets(buffer, datalen, fpcoords);
	fgets(buffer, datalen, fpcoords);
	for (ii=0; ii<natom; ii++)
	{
		fscanf(fpcoords, "%s %lf %lf %lf\n", buffer, &xx[ii], &yy[ii], &zz[ii]);
	}
	fclose(fpcoords);

	return 0;
}
