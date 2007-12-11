/*
 * Functions related to HMC operation go here
 *
 * Written by Yang Wang 12-07-2007
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <ctype.h>
#include "vars.h"

extern int vver();
extern int erfrc();
extern int rafrc();
extern int averages();
extern int loadit();
extern int saveit();
extern int printit();
extern int snapshot();
extern int trajectory();
extern int velinit();
extern int echo();

int init_hmc()
{
	int ii, jj;
	const int datalen = 200;
	char buffer[200];
	double auc, buc, cuc, nanotuberadius;
	int na, nb, nc;
	int nspecies_yang;
	int tnuatoms;
	char szgrid[12][200]; // assuming only 12 types of grids
	char keyword[100];

	fprintf(stderr,"Reading input data for HMC simulation...\n");
	fprintf(fpouts, "Reading input data for HMC simulation...\n");

	// re-open input file to read extra data section
	fpins = fopen(INPUT,"r");

	while (fgets(buffer, datalen, fpins)!=NULL)
	{
		sscanf(buffer, "%s", keyword);
		for (ii=0; ii<strlen(keyword); ii++)
			keyword[ii] = toupper(keyword[ii]);
		if (!strcmp(keyword, "HMC"))
		{
			fprintf(stderr,"Data section for HMC simulation found...\n");
			fprintf(fpouts, "Data section for HMC simulation found...\n");

			sscanf(fgets(buffer, datalen, fpins), "%d", &nstep_md_per_hmc);
			sscanf(fgets(buffer, datalen, fpins), "%lf %lf %lf", &prob_cm,
					&prob_vc, &prob_id);
			sscanf(fgets(buffer, datalen, fpins), "%lf %lf", &ratio_cm_req,
					&ratio_vc_req);
			sscanf(fgets(buffer, datalen, fpins), "%d %d",
					&nstep_delt_adj_cycle, &nstep_delv_adj_cycle);
			sscanf(fgets(buffer, datalen, fpins), "%lf", &delv);

			if (prob_cm+prob_vc+prob_id>1.0)
			{
				fprintf( stderr, "Warning: prob_cm (%lf) + prob_vc (%lf) + prob_id (%lf) > 1.0 \n", prob_cm, prob_vc, prob_id);
				fprintf(
						fpouts,
						"Warning: prob_cm (%lf) + prob_vc (%lf) + prob_id (%lf) > 1.0 \n",
						prob_cm, prob_vc, prob_id);

			}
			fclose(fpins);
			return 0;
		} // if keyword found
	} // read through the lines
	fprintf(stderr,"Error: data for HMC simulation not found.\n");
	fprintf(fpouts, "Error: data for HMC simulation not found.\n");
	fclose(fpins);
	exit(1);
}

int hmc()
{
	int ii;
	double prob_vc_upper;
	double *xx_old, *yy_old, *zz_old;
	double upot_old, utot_old;
	double upot_new, utot_new;
	double ukin_old, tinst_old, uinter_old, uintra_old, uvdw_old, ubond_old,
			uangle_old, udih_old, uimp_old, uewald_old, usflj_old, unhts_old,
			unhtss_old, virial_old, virial_inter_old, virial_intra_old,
			utsbs_old, pinst_old, boxlx_old, boxly_old, boxlz_old, boxv_old;
	double dH;
	int isAccept;
	double rndnum[3];
	double ratio;

	// calculate the upper limit of volume changes
	prob_vc_upper = prob_cm + prob_vc;

	// allocate memory for saving positions
	_safealloc(xx_old,natom,sizeof(double))
	;
	_safealloc(yy_old,natom,sizeof(double))
	;
	_safealloc(zz_old,natom,sizeof(double))
	;

	// start of the HMC simulation
	// initialize velocities
	fprintf(fpouts, "initializing velocities...\n");
	velinit();

	// if not new run, load from old file
	if (fStart_option!=new_run)
		loadit();

	// calculate total energies
	erfrc();
	rafrc();

	// the energies
	upot_old = uinter + uintra;
	utot_old = upot_old + ukin;
	// add energy of thermostat, if nose hoover is not used, they will just be zero
	utot_old = utot_old + unhts + unhtss + utsbs;
	// add long range corrections into total energy
	utot_old = utot_old + uljlrc;

	// save old positions and energies
	for (ii=0; ii<natom; ii++)
	{
		xx_old[ii] = xx[ii];
		yy_old[ii] = yy[ii];
		zz_old[ii] = zz[ii];
	}
	ukin_old = ukin;
	tinst_old = tinst;
	uinter_old = uinter;
	uintra_old = uintra;
	uvdw_old = uvdw;
	ubond_old = ubond;
	uangle_old = uangle;
	udih_old = udih;
	uimp_old = uimp;
	uewald_old = uewald;
	usflj_old = usflj;
	unhts_old = unhts;
	unhtss_old = unhtss_old;
	virial_old = virial;
	virial_inter_old = virial_inter;
	virial_intra_old = virial_intra;
	utsbs_old = utsbs;
	pinst_old = pinst;
	boxlx_old = boxlx;
	boxly_old = boxly;
	boxlz_old = boxlz;
	boxv_old = boxv;

	// print out initial values
	echo();
	// print initial properties
	printit();
	// make snapshots & movies
	trajectory();

	// simulation loop
	for (istep=nstep_start; istep<=nstep; istep++) // NOTE: start from 1 and <=
	{
		ranmar(rndnum, 1);

		if (rndnum[0] <= prob_cm) // canonical moves
		{
			// vel init and energy calcualtions have been done for the 1st step outside the loop
			if (istep!=nstep_start)
			{
				velinit();
				erfrc();
				rafrc();
				// the energies
				upot_old = uinter + uintra;
				utot_old = upot_old + ukin;
				// add energy of thermostat, if nose hoover is not used, they will just be zero
				utot_old = utot_old + unhts + unhtss + utsbs;
				// add long range corrections into total energy
				utot_old = utot_old + uljlrc;

				// save old positions and energies
				for (ii=0; ii<natom; ii++)
				{
					xx_old[ii] = xx[ii];
					yy_old[ii] = yy[ii];
					zz_old[ii] = zz[ii];
				}
				ukin_old = ukin;
				tinst_old = tinst;
				uinter_old = uinter;
				uintra_old = uintra;
				uvdw_old = uvdw;
				ubond_old = ubond;
				uangle_old = uangle;
				udih_old = udih;
				uimp_old = uimp;
				uewald_old = uewald;
				usflj_old = usflj;
				unhts_old = unhts;
				unhtss_old = unhtss_old;
				virial_old = virial;
				virial_inter_old = virial_inter;
				virial_intra_old = virial_intra;
				utsbs_old = utsbs;
				pinst_old = pinst;
				boxlx_old = boxlx;
				boxly_old = boxly;
				boxlz_old = boxlz;
				boxv_old = boxv;
			}
			// MD moves
			for (ii=0; ii<nstep_md_per_hmc; ii++)
			{
				vver();
			}
			upot_new = uinter + uintra;
			utot_new = upot_new + ukin;
			// add energy of thermostat, if nose hoover is not used, they will just be zero
			utot_new = utot_new + unhts + unhtss + utsbs;
			// add long range corrections into total energy
			utot_new = utot_new + uljlrc;

			// Hamotonial difference
			dH = (utot_new - utot_old)*rRgas/treq;
			isAccept = 0;
			if (dH <= 0.0)
			{
				isAccept = 1;
			}
			else
			{
				ranmar(rndnum, 1);
				if (rndnum[0] < exp(-dH))
					isAccept = 1;
			}

			icounter[20]++; // canonical moves
			if (isAccept == 1)
			{
				icounter[21]++; // accepted canonical moves
			}
			else
			{
				// restore old values
				for (ii=0; ii<natom; ii++)
				{
					xx[ii] = xx_old[ii];
					yy[ii] = yy_old[ii];
					zz[ii] = zz_old[ii];
				}
				ukin = ukin_old;
				tinst = tinst_old;
				uinter = uinter_old;
				uintra = uintra_old;
				uvdw = uvdw_old;
				ubond = ubond_old;
				uangle = uangle_old;
				udih = udih_old;
				uimp = uimp_old;
				uewald = uewald_old;
				usflj = usflj_old;
				unhts = unhts_old;
				unhtss = unhtss_old;
				virial = virial_old;
				virial_inter = virial_inter_old;
				virial_intra = virial_intra_old;
				utsbs = utsbs_old;
				pinst = pinst_old;
				boxlx = boxlx_old;
				boxly = boxly_old;
				boxlz = boxlz_old;
				boxv = boxv_old;
			}
		}
		else if (rndnum[0]<=prob_vc_upper) // volume change moves
		{
			// volchg(pBox);
		}
		else // insertions or deletions
		{
			// insdel(pBox);
		}

		// print out, snapshot, trajectory, save
		if (istep%nstep_print == 0)
		{
			printit();
		}
		if (nstep_ss && istep%nstep_ss == 0)
		{
			snapshot();
		}
		if (nstep_trj && istep%nstep_trj==0)
		{
			trajectory();
		}
		if (istep%nstep_save==0)
		{
			saveit();
		}

		icounter[11]--;
		// if still in equilibrium run
		// adjust the delt, deltv
		// do not do averages
		if (icounter[11]>=0)
		{
			if (icounter[20]%nstep_delt_adj_cycle==0) // delt adjustment
			{
				ratio = icounter[21]*1.0/nstep_delt_adj_cycle;
				delt *= (1.0 - ratio_cm_req + ratio);
				icounter[20] = 0;
				icounter[21] = 0;
				// delt update
				deltby2 = delt/2.0;
				delts = delt/nstep_inner;
				deltsby2 = delts/2.0;
			}
			if (icounter[22]%nstep_delv_adj_cycle==0) // delv adjustment
			{
				ratio = icounter[23]*1.0/nstep_delv_adj_cycle;
				delv = delv*(1.0 - ratio_vc_req + ratio);
				icounter[22] = 0;
				icounter[23] = 0;
				// length update
			}

			continue;
		}

		// accumulators
		do_accumu();

		// do averages
		if ((istep-nstep_eq)%nstep_ave==0)
		{
			averages();
		}
	}

	// release the dynamically allocated memory for saving old positions
	free(xx_old);
	free(yy_old);
	free(zz_old);

	return 0;
}
