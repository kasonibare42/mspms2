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
int InitCheckUniques()
{
	int mm, ii, jj;
	PSAMPLE_MOLECULE pSampleMole;
	int iFirst_Dih, iLast_Dih, iFirst_Agl, iLast_Agl, iFirst_Bnd, iLast_Bnd;

	// Check unique for dihedrals.
	// This is for possible ring structures where 1,4 atoms can form multiple dihedrals.
	// e.g. 1-2-3-4
	//       \5-6/
	// 1234 and 1564
	// The 14 pair should only be calculated once for energy/force.
	// That is what the unique check is for.
	for (mm=0;mm<nspecie;mm++)
	{
		pSampleMole = sample_mole + mm;
		for (ii=0;ii<pSampleMole->ndih;ii++)
		{
			pSampleMole->isDih_unique[ii] = true;
		}
		for (ii=0; ii<pSampleMole->ndih-1; ii++) // ii diehdral
		{
			if (pSampleMole->isDih_unique[ii])
			{
				for (jj=ii+1; jj<pSampleMole->ndih; jj++) // jj dihedral
				{
					if (pSampleMole->isDih_unique[jj])
					{
						if ((pSampleMole->dih_idx[ii][0]==pSampleMole->dih_idx[jj][0] && pSampleMole->dih_idx[ii][3]==pSampleMole->dih_idx[jj][3])
								|| (pSampleMole->dih_idx[ii][0]==pSampleMole->dih_idx[jj][3] && pSampleMole->dih_idx[ii][3]==pSampleMole->dih_idx[jj][0]))
						{
							pSampleMole->isDih_unique[jj] = false;
							fprintf(
									fpouts,
									"Sample molecule %d: Dihedral %d and dihedral %d have the same ending pairs.\n",
									mm, ii, jj);
						} // Check uniqueness
					} // End if dihedral jj is unique
				} // End of dihedral jj
			} // End if dihedral ii is unique
		} // End of Dihedral ii
		
		// following codes make sure 14 and 13 do not share same ending pairs
		// this is also for ring kind structures
		for (ii=0; ii<pSampleMole->ndih; ii++)
		{
			if (pSampleMole->isDih_unique[ii])
			{
				for (jj=0; jj<pSampleMole->nangle; jj++)
				{
					if ((pSampleMole->dih_idx[ii][0]==pSampleMole->agl_idx[jj][0] && pSampleMole->dih_idx[ii][3]==pSampleMole->agl_idx[jj][2])
							|| (pSampleMole->dih_idx[ii][0]==pSampleMole->agl_idx[jj][2] && pSampleMole->dih_idx[ii][3]==pSampleMole->agl_idx[jj][0]))
					{
						pSampleMole->isDih_unique[jj] = false;
						fprintf(
								fpouts,
								"Sample molecule %d: Dihedral %d and angle %d have the same ending pairs.\n",
								mm, ii, jj);
					} // Check uniqueness
				} // End of Angle loop
			} // If dihedral is unique
		} // End of Dihedral loop

		// following codes make sure 14 and 12 do not share the same ending pairs
		for (ii=0; ii<pSampleMole->ndih; ii++)
		{
			if (pSampleMole->isDih_unique[ii])
			{
				for (jj=0; jj<pSampleMole->nbond; jj++)
				{
					if ((pSampleMole->dih_idx[ii][0]==pSampleMole->bnd_idx[jj][0] && pSampleMole->dih_idx[ii][3]==pSampleMole->bnd_idx[jj][1])
							|| (pSampleMole->dih_idx[ii][0]==pSampleMole->bnd_idx[jj][1] && pSampleMole->dih_idx[ii][3]==pSampleMole->bnd_idx[jj][0]))
					{
						pSampleMole->isDih_unique[jj] = false;
						fprintf(
								fpouts,
								"Sample molecule %d: Dihedral %d and bond %d have the same ending pairs.\n",
								mm, ii, jj);
					} // End of check of uniqueness
				} // End of Bond loop
			} // End if dihedral unique
		} // End of dihedral loop

		// check unique for angles
		// see above comments for dihedrals
		for (ii=0; ii<pSampleMole->nangle; ii++)
		{
			pSampleMole->isAngle_unique[ii] = true;
		}
		for (ii=0; ii<pSampleMole->nangle-1; ii++)
		{
			if (pSampleMole->isAngle_unique[ii])
			{
				for (jj=ii+1; jj<pSampleMole->nangle; jj++)
				{
					if (pSampleMole->isAngle_unique[jj])
					{
						if ((pSampleMole->agl_idx[ii][0]==pSampleMole->agl_idx[jj][0] && pSampleMole->agl_idx[ii][2]==pSampleMole->agl_idx[jj][2])
								|| (pSampleMole->agl_idx[ii][0]==pSampleMole->agl_idx[jj][2] && pSampleMole->agl_idx[ii][2]==pSampleMole->agl_idx[jj][0]))
						{
							pSampleMole->isAngle_unique[jj] = false;
							fprintf(
									fpouts,
									"Sample molecule %d: Angle %d and angle %d have the same ending pairs.\n",
									mm, ii, jj);
						} // Check the angle uniqueness
					} // End angle jj uniqueness
				} // End of angle jj loop
			} // End angle ii uniqueness
		} // End of angle ii loop

		// following codes make sure 13 and 12 do not share the same ending pairs
		for (ii=0; ii<pSampleMole->nangle; ii++)
		{
			if (pSampleMole->isAngle_unique[ii])
			{
				for (jj=0; jj<pSampleMole->nbond; jj++)
				{
					if ((pSampleMole->agl_idx[ii][0]==pSampleMole->bnd_idx[jj][0] && pSampleMole->agl_idx[ii][2]==pSampleMole->bnd_idx[jj][1])
							|| (pSampleMole->agl_idx[ii][0]==pSampleMole->bnd_idx[jj][1] && pSampleMole->agl_idx[ii][2]==pSampleMole->bnd_idx[jj][0]))
					{
						pSampleMole->isAngle_unique[jj] = false;
						fprintf(
								fpouts,
								"Sample molecule %d: Angle %d and bond %d have the same ending pairs.\n",
								mm, ii, jj);
					} // End of angle uniqueness check
				} // End of bond loop
			} // End of Angle uniqueness
		} // End of Angle loop
	}

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
	int specie_id, sample_atom_id;

	px = py = pz = 0.0;
	// mv^2=kT, So, p = sqrt(kTm), v = sqrt(kT/m)
	// Use reduced units, v* = sqrt(T*/m*)
	stdvtmp = sqrt(treq);
	
	for (ii=0; ii<natom; ii++)
	{
		// From index ii, we calculate which specie this atom belongs to and
		// its position within a molecule
		get_specie_and_relative_atom_id(ii, &specie_id, &sample_atom_id);
		stdv = stdvtmp/sqrt(sample_mole[specie_id].aw[sample_atom_id]);
		vx[ii] = stdv*gaussran();
		vy[ii] = stdv*gaussran();
		vz[ii] = stdv*gaussran();
		px += vx[ii]*sample_mole[specie_id].aw[sample_atom_id];
		py += vy[ii]*sample_mole[specie_id].aw[sample_atom_id];
		pz += vz[ii]*sample_mole[specie_id].aw[sample_atom_id];
	}
	// zero the momentum
	px /= system_mass;
	py /= system_mass;
	pz /= system_mass;
	for (ii=0; ii<natom; ii++)
	{
		get_specie_and_relative_atom_id(ii, &specie_id, &sample_atom_id);
		vx[ii] -= px;
		vy[ii] -= py;
		vz[ii] -= pz;
	}
	
	px = 0.0;
	py = 0.0;
	pz = 0.0;
	for (ii=0; ii<natom; ii++)
	{
		// From index ii, we calculate which specie this atom belongs to and
		// its position within a molecule
		get_specie_and_relative_atom_id(ii, &specie_id, &sample_atom_id);
		px += vx[ii]*sample_mole[specie_id].aw[sample_atom_id];
		py += vy[ii]*sample_mole[specie_id].aw[sample_atom_id];
		pz += vz[ii]*sample_mole[specie_id].aw[sample_atom_id];
	}
	
	// Calculate the energy and instantaneous temperature
	ukin = 0.0;
	for (ii=0; ii<natom; ii++)
	{
		get_specie_and_relative_atom_id(ii, &specie_id, &sample_atom_id);
		ukin += sample_mole[specie_id].aw[sample_atom_id]*(vx[ii]*vx[ii]+vy[ii]*vy[ii]+vz[ii]*vz[ii]);
	}
	ukin = 0.5*ukin;
	tinst = 2.0*ukin/nfree;
	
	// Rescale velocity for required temperature
	double scaling;
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
		get_specie_and_relative_atom_id(ii, &specie_id, &sample_atom_id);
		ukin += sample_mole[specie_id].aw[sample_atom_id]*(vx[ii]*vx[ii]+vy[ii]*vy[ii]+vz[ii]*vz[ii]);
	}
	ukin = 0.5*ukin;
	tinst = 2.0*ukin/nfree;
	
	return 0;
}

// Read and initialize the necessary variables for HMC
int init_hmc()
{
	int ii;
	char buffer[LONG_STRING_LENGTH];
	char keyword[100];
	int position_counter;

	fprintf(stderr,"Reading input data for HMC simulation...\n");
	fprintf(fpouts, "Reading input data for HMC simulation...\n");

	// re-open input file to read extra data section
	fpins = fopen(INPUT,"r");

	while (fgets(buffer, LONG_STRING_LENGTH, fpins)!=NULL)
	{
		sscanf(buffer, "%s", keyword);
		for (ii=0; ii<strlen(keyword); ii++)
		{
			keyword[ii] = toupper(keyword[ii]);
		}
		if (!strcmp(keyword, "HMC"))
		{
			fprintf(stderr,"Data section for HMC simulation found...\n");
			fprintf(fpouts, "Data section for HMC simulation found...\n");

			sscanf(fgets(buffer, LONG_STRING_LENGTH, fpins), "%d", &nstep_md_per_hmc);
			sscanf(fgets(buffer, LONG_STRING_LENGTH, fpins), "%lf", &pdisp);
			sscanf(fgets(buffer, LONG_STRING_LENGTH, fpins), "%lf", &pvolm);
			sscanf(fgets(buffer, LONG_STRING_LENGTH, fpins), "%lf", &pmake);
			sscanf(fgets(buffer, LONG_STRING_LENGTH, fpins), "%lf", &pkill);
			sscanf(fgets(buffer, LONG_STRING_LENGTH, fpins), "%lf", &rreq_disp);
			sscanf(fgets(buffer, LONG_STRING_LENGTH, fpins), "%lf", &rreq_volm);
			sscanf(fgets(buffer, LONG_STRING_LENGTH, fpins), "%lf", &delv);
			sscanf(fgets(buffer, LONG_STRING_LENGTH, fpins), "%d", &nstep_delt_adj_cycle);
			sscanf(fgets(buffer, LONG_STRING_LENGTH, fpins), "%d", &nstep_delv_adj_cycle);

			// Read insertion/deletion input data if required
			if (pmake+pkill > 0.0)
			{
				for (ii=0; ii<nspecie; ii++)
				{
					sscanf(fgets(buffer, LONG_STRING_LENGTH, fpins), "%lf %lf", &cpt[ii], &pcomp[ii]);
					// initialize zact
					zact[ii] = exp(cpt[ii]/treq)*boxv;
				}
			}

			// Allocate memory for saving positions
			_safealloc(xx_old,natom,sizeof(double)) ;
			_safealloc(yy_old,natom,sizeof(double)) ;
			_safealloc(zz_old,natom,sizeof(double)) ;

			fclose(fpins);
			return 0;
		} // if keyword found
	} // read through the lines
	fprintf(stderr,"Error: Data for HMC simulation not found.\n");
	fprintf(fpouts, "Error: Data for HMC simulation not found.\n");
	fclose(fpins);
	exit(1);
}

/// Read parameters from the input file for MDNVT and initialize the variables needed for NVT simulation
int init_nvt()
{
	int ii;
	char buffer[LONG_STRING_LENGTH];
	char keyword[100];

	fprintf(stderr,"Reading input data for MD NVT simulation...\n");
	fprintf(fpouts, "Reading input data for MD NVT simulation...\n");

	// re-open input file to read extra data section
	fpins = fopen(INPUT,"r");

	while (fgets(buffer, LONG_STRING_LENGTH, fpins)!=NULL)
	{
		sscanf(buffer, "%s", keyword);
		for (ii=0; ii<strlen(keyword); ii++)
		{
			keyword[ii] = toupper(keyword[ii]);
		}
		if (!strcmp(keyword, "MDNVT"))
		{
			fprintf(stderr,"Data section for MD NVT simulation found...\n");
			fprintf(fpouts, "Data section for MD NVT simulation found...\n");
			sscanf(fgets(buffer, LONG_STRING_LENGTH, fpins), "%lf %lf", &qq, &qqs);

			/**
			 * Qts = RGAS*treq*nfree/Omega
			 * where Omega is a parameter related to the mass of the thermostat
			 * for this program, we read in the Qts directly.
			 */
			// variables for nose-hoover NVT, see frenkel and smit
			delt_sqby2 = delt*delt/2.0;
			delts_sqby2 = delts*delts/2.0;
			unhts = 0.0;
			gg = nfree; // need double check
			ss = 0.0;
			ps = 0.0;
			ggs = nfree;
			sss = 0.0;
			pss = 0.0;
			unhtss = 0.0;

			fclose(fpins);
			return 0;
		} // if keyword found
	} // read through lines
	fprintf(stderr,"Error: data for MD NVT not found.\n");
	fprintf(fpouts, "Error: data for MD NVT not found.\n");
	fclose(fpins);
	exit(1);
}

// Read and initialize NPT related variables
int init_npt_respa()
{
	int ii;
	char buffer[LONG_STRING_LENGTH];
	char keyword[100];

	fprintf(stderr,"Reading input data for MD NPT simulation...\n");
	fprintf(fpouts, "Reading input data for MD NPT simulation...\n");

	// re-open input file to read extra data section
	fpins = fopen(INPUT,"r");

	while (fgets(buffer, LONG_STRING_LENGTH, fpins)!=NULL)
	{
		sscanf(buffer, "%s", keyword);
		for (ii=0; ii<strlen(keyword); ii++)
		{
			keyword[ii] = toupper(keyword[ii]);
		}
		if (!strcmp(keyword, "MDNPT"))
		{
			fprintf(stderr,"Data section for MD NPT simulation found...\n");
			fprintf(fpouts, "Data section for MD NPT simulation found...\n");

			sscanf(fgets(buffer, LONG_STRING_LENGTH, fpins), "%lf %lf", &Qts, &Qbs);

			// thermo/barostat
			utsbs = 0.0;
			vts = sqrt((nfree+1.0)/Qts);
			vbs = sqrt((nfree+1.0)/Qbs);
			rts = 0.0;

			utsbs = 0.5*Qbs*vbs*vbs + 0.5*Qts*vts*vts + (nfree+1)*treq*rts
					+ preq*boxv;

			fclose(fpins);
			return 0;
		} // if keyword found
	} // read through lines
	fprintf(stderr,"Error: data for MD NPT not found.\n");
	fprintf(fpouts, "Error: data for MD NPT not found.\n");
	fclose(fpins);
	exit(1);
}

/// Read in electrostatic parametes and initialize related variables
int init_charge()
{
	int ii;
	char buffer[LONG_STRING_LENGTH];
	char keyword[100];

	fprintf(stderr,"Reading input data for Electrostatic interactions...\n");
	fprintf(fpouts, "Reading input data for Electrostatic interactions...\n");

	// re-open input file to read extra data section
	fpins = fopen(INPUT,"r");

	while (fgets(buffer, LONG_STRING_LENGTH, fpins)!=NULL)
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

			if (iChargeType == ELECTROSTATIC_EWALD)
			{
				sscanf(fgets(buffer, LONG_STRING_LENGTH, fpins), "%lf", &kappa);
				sscanf(fgets(buffer, LONG_STRING_LENGTH, fpins), "%d", &KMAXX);
				sscanf(fgets(buffer, LONG_STRING_LENGTH, fpins), "%d", &KMAXY);
				sscanf(fgets(buffer, LONG_STRING_LENGTH, fpins), "%d", &KMAXZ);
				sscanf(fgets(buffer, LONG_STRING_LENGTH, fpins), "%d", &KSQMAX);
				sscanf(fgets(buffer, LONG_STRING_LENGTH, fpins), "%d", &fEwald_BC);
				sscanf(fgets(buffer, LONG_STRING_LENGTH, fpins), "%d", &fEwald_Dim);

				// Initialization
				// set ewald parameters
				kappa = kappa*sigma_base; // Reduce the kappa
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
				sscanf(fgets(buffer, LONG_STRING_LENGTH, fpins), "%lf", &kappa);
				
				// set wolf parameters
				kappa = kappa*sigma_base; // Reduce the kappa
				wolfvcon1 = -erfc(kappa*rcutoffelec)/rcutoffelec;
				wolfvcon2 = erfc(kappa*rcutoffelec)/rcutoffelecsq + 2.0*kappa
						*exp(-(kappa *rcutoffelec)*(kappa*rcutoffelec))
						/(sqrt(pi) *rcutoffelec);
				wolfvcon3 = kappa/sqrt(pi)+erfc(kappa*rcutoffelec)/(2.0*rcutoffelec);
				wolffcon1 = 2.0*kappa/sqrt(pi);
				wolffcon2 = -wolfvcon2;
			}
			else if (iChargeType == ELECTROSTATIC_SIMPLE_COULOMB)
			{
			}
			
			// calculate the reduced coulomb constant
			coulomb_prefactor = COULOMB_PREFACTOR/epsilon_base/sigma_base;

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
			for (ii=0; ii<sample_mole[mm].natom; ii++)
			{
				// use the first molecule of one specie to do the calculation
				atomid_1 = ii;
				// loop through all the atoms in one molecule of specie nn
				for (jj=0; jj<sample_mole[nn].natom; jj++)
				{
					atomid_2 = jj;
					sigmaij = 0.5*(sample_mole[mm].sigma[atomid_1]+sample_mole[nn].sigma[atomid_2]);
					epsilonij = sqrt(sample_mole[mm].epsilon[atomid_1]*sample_mole[nn].epsilon[atomid_2]);
					temp1 = pow(sigmaij, 9.0)*uljlrc_term1 + pow(sigmaij, 3.0)*uljlrc_term2;
					temp2 = epsilonij*pow(sigmaij, 3.0);
					temp3 = pow(sigmaij, 9.0)*pljlrc_term1 + pow(sigmaij, 3.0)*pljlrc_term2;
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

	fprintf(stderr,"Initializing variables...\n");
	fprintf(fpouts, "Initializing variables...\n");
	
	for (ii=0;ii<natom;ii++)
	{
		fxl[ii] = fyl[ii] = fzl[ii] = 0.0;
		fxs[ii] = fys[ii] = fzs[ii] = 0.0;
	}
	uvdw = 0.0;
	ushift = 0.0;
	uelec = 0.0;
	ureal = 0.0;
	uexcl = 0.0;
	ucoulomb = 0.0;
	virial_inter = 0.0;
	uinter = 0.0;
	ubond = 0.0;
	uangle = 0.0;
	udih = 0.0;
	uimp = 0.0;
	virial_intra = 0.0;
	uintra = 0.0;

	/// Initiate file variables, LOG, TRAJECTORY\n
	/// Output file is initialized already at the very beginning of the run.
	fplog = fopen(LOG,"w");
	fptrj = fopen(TRAJECTORY,"wb");

	/// Calculate molecule weight for the real list.
	system_mass = 0.0;
	for (ii=0; ii<nspecie; ii++)
	{
		system_mass += sample_mole[ii].mw*nmole_per_specie[ii];
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
	InitCheckUniques();
	
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

	// initialize velocities for needed simulations
	fprintf(stderr, "initializing velocities...\n");
	fprintf(fpouts, "initializing velocities...\n");
	velinit();
	
	// initialize thermostat/baron stat input data
	if (what_ensemble == NVT)
	{
		init_nvt();
	}
	else if (what_ensemble == NPT)
	{
		init_npt_respa();
	}

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
		init_charge();
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
