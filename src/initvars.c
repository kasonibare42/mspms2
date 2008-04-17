#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <ctype.h>
#include <time.h>
#include "mspms2.h"

/**
 * Check if any dihedral, angle, bond share the same ending pairs.
 */
int CheckUniques()
{
	int ii, jj;
	int iFirst_Dih, iLast_Dih, iFirst_Agl, iLast_Agl, iFirst_Bnd, iLast_Bnd;

	// Check unique for dihedrals.
	// This is for possible ring structures where 1,4 atoms can form multiple dihedrals.
	// e.g. 1-2-3-4
	//       \5-6/
	// 1234 and 1564
	// The 14 pair should only be calculated once for energy/force.
	// That is what the unique check is for.
	iFirst_Dih = sample_mole_first_dih_idx[0]; // The very first dihedral
	iLast_Dih = sample_mole_last_dih_idx[nspecie-1]; // The very last dihedral
	iFirst_Agl = sample_mole_first_angle_idx[0]; // The very first angle
	iLast_Agl = sample_mole_last_angle_idx[nspecie-1]; // The very last angle
	iFirst_Bnd = sample_mole_first_bond_idx[0]; // The very first bond
	iLast_Bnd = sample_mole_last_bond_idx[nspecie-1]; // The very last bond
	for (ii=iFirst_Dih; ii<iLast_Dih; ii++)
	{
		sample_isDih_unique[ii] = true;
	}
	for (ii=iFirst_Dih; ii<iLast_Dih-1; ii++)
	{
		if (sample_isDih_unique[ii])
		{
			for (jj=ii+1; jj<iLast_Dih; jj++)
			{
				if (sample_isDih_unique[jj])
				{
					if ((sample_dih_idx[ii][0]==sample_dih_idx[jj][0] && sample_dih_idx[ii][3]==sample_dih_idx[jj][3])
							|| (sample_dih_idx[ii][0]==sample_dih_idx[jj][3] && sample_dih_idx[ii][3]==sample_dih_idx[jj][0]))
					{
						sample_isDih_unique[jj] = false;
						fprintf(
								fpouts,
								"Dihedral %d and dihedral %d have the same ending pairs.\n",
								ii, jj);
					} // Check uniqueness
				} // End if dihedral jj is unique
			} // End of dihedral jj
		} // End if dihedral ii is unique
	} // End of Dihedral ii

	// following codes make sure 14 and 13 do not share same ending pairs
	// this is also for ring kind structures
	for (ii=iFirst_Dih; ii<iLast_Dih; ii++)
	{
		if (sample_isDih_unique[ii])
		{
			for (jj=iFirst_Agl; jj<iLast_Agl; jj++)
			{
				if ((sample_dih_idx[ii][0]==sample_angle_idx[jj][0] && sample_dih_idx[ii][3]==sample_angle_idx[jj][2])
						|| (sample_dih_idx[ii][0]==sample_angle_idx[jj][2] && sample_dih_idx[ii][3]==sample_angle_idx[jj][0]))
				{
					sample_isDih_unique[jj] = false;
					fprintf(
							fpouts,
							"Dihedral %d and angle %d have the same ending pairs.\n",
							ii, jj);
				} // Check uniqueness
			} // End of Angle loop
		} // If dihedral is unique
	} // End of Dihedral loop

	// following codes make sure 14 and 12 do not share the same ending pairs
	for (ii=iFirst_Dih; ii<iLast_Dih; ii++)
	{
		if (sample_isDih_unique[ii])
		{
			for (jj=iFirst_Bnd; jj<iLast_Bnd; jj++)
			{
				if ((sample_dih_idx[ii][0]==sample_bond_idx[jj][0] && sample_dih_idx[ii][3]==sample_bond_idx[jj][1])
						|| (sample_dih_idx[ii][0]==sample_bond_idx[jj][1] && sample_dih_idx[ii][3]==sample_bond_idx[jj][0]))
				{
					sample_isDih_unique[jj] = false;
					fprintf(
							fpouts,
							"Dihedral %d and bond %d have the same ending pairs.\n",
							ii, jj);
				} // End of check of uniqueness
			} // End of Bond loop
		} // End if dihedral unique
	} // End of dihedral loop

	// check unique for angles
	// see above comments for dihedrals
	for (ii=iFirst_Agl; ii<iLast_Agl; ii++)
	{
		sample_isAngle_unique[ii] = true;
	}
	for (ii=iFirst_Agl; ii<iLast_Agl-1; ii++)
	{
		if (sample_isAngle_unique[ii])
		{
			for (jj=ii+1; jj<iLast_Agl; jj++)
			{
				if (sample_isAngle_unique[jj])
				{
					if ((sample_angle_idx[ii][0]==sample_angle_idx[jj][0] && sample_angle_idx[ii][2]==sample_angle_idx[jj][2])
							|| (sample_angle_idx[ii][0]==sample_angle_idx[jj][2] && sample_angle_idx[ii][2]==sample_angle_idx[jj][0]))
					{
						sample_isAngle_unique[jj] = false;
						fprintf(
								fpouts,
								"Angle %d and angle %d have the same ending pairs.\n",
								ii, jj);
					} // Check the angle uniqueness
				} // End angle jj uniqueness
			} // End of angle jj loop
		} // End angle ii uniqueness
	} // End of angle ii loop

	// following codes make sure 13 and 12 do not share the same ending pairs
	for (ii=iFirst_Agl; ii<iLast_Agl; ii++)
	{
		if (sample_isAngle_unique[ii])
		{
			for (jj=iFirst_Bnd; jj<iLast_Bnd; jj++)
			{
				if ((sample_angle_idx[ii][0]==sample_bond_idx[jj][0] && sample_angle_idx[ii][2]==sample_bond_idx[jj][1])
						|| (sample_angle_idx[ii][0]==sample_bond_idx[jj][1] && sample_angle_idx[ii][2]==sample_bond_idx[jj][0]))
				{
					sample_isAngle_unique[jj] = false;
					fprintf(
							fpouts,
							"Angle %d and bond %d have the same ending pairs.\n",
							ii, jj);
				} // End of angle uniqueness check
			} // End of bond loop
		} // End of Angle uniqueness
	} // End of Angle loop

	return 0;
}

/**
 * Initialize velocities
 */
int velinit()
{
	int ii;
	double px, py, pz;
	double stdvtmp, stdv;
	double totalmass;
	double scaling;
	int specie_id, rela_atom_id;

	totalmass = 0.0;
	px = py = pz = 0.0;
	// mv^2=kT, So, p = sqrt(kTm), v = sqrt(kT/m)
	// Use reduced units, v* = sqrt(T*/m*)
	stdvtmp = sqrt(treq);
	

	for (ii=0; ii<natom; ii++)
	{
		// From index ii, we calculate which specie this atom belongs to and
		// its position within a molecule
		get_specie_and_relative_atom_id(ii, &specie_id, &rela_atom_id);
		// stdv = stdvtmp/sqrt(sample_aw[]);
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
	char buffer[STRING_LENGTH];
	char keyword[100];

	fprintf(stderr,"Reading input data for Electrostatic interactions...\n");
	fprintf(fpouts, "Reading input data for Electrostatic interactions...\n");

	// re-open input file to read extra data section
	fpins = fopen(INPUT,"r");

	while (fgets(buffer, STRING_LENGTH, fpins)!=NULL)
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
				sscanf(fgets(buffer, STRING_LENGTH, fpins), "%lf", &kappa);
				sscanf(fgets(buffer, STRING_LENGTH, fpins), "%d %d %d %d",
						&KMAXX, &KMAXY, &KMAXZ, &KSQMAX);
				sscanf(fgets(buffer, STRING_LENGTH, fpins), "%d %d",
						&fEwald_BC, &fEwald_Dim);

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
				sscanf(fgets(buffer, STRING_LENGTH, fpins), "%lf", &kappa);

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
	// Set the proper end step to eq steps or data taking steps
	nstep_end = (nstep_eq>0) ? nstep_eq : nstep;

	// Zero the counts and accums
	rezero();

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
	shift4 = 4.0*shift1;

	// volume calculation
	boxv = boxlx*boxly*boxlz;

	// delt related
	deltby2 = delt/2.0;
	delts = delt/nstep_inner;
	deltsby2 = delts/2.0;
	dt_outer2 = deltby2;
	dt_outer4 = delt/4.0;
	dt_outer8 = delt/8.0;
	
	// Check if any bond, angle, dihedral share the same ending pairs.
	CheckUniques();

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
	fprintf(stderr, "initializing velocities...\n");
	fprintf(fpouts, "initializing velocities...\n");
	velinit();

	// If Solid-fluid interaction is required, initiliaze the related variables.
	if (iSF_type==SF_NANOTUBE_HYPERGEO)
	{
		init_sf_hypergeo();
	}
	else if (iSF_type==SF_NANOTUBE_ATOM_EXPLICIT)
	{
		init_sf_atom_explicit();
	}
	else if (iSF_type==SF_NANOTUBE_TASOS)
	{
		init_tasos_grid();
	}
	else if (iSF_type==SF_NANOTUBE_MY_INTERP)
	{
		init_my_interp();
	}

	// Read in electrostatic parametes and initialize if needed
	if (iChargeType != ELECTROSTATIC_NONE)
	{
		fnInitCharge();
	}

	// Initialize the common terms and calculate LJ lrc if it is requested.
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
	if (iStart_option!=NEW)
	{
		loadit();
	}

	return 0;
}
