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
#include "funcs.h"

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
		{
			keyword[ii] = toupper(keyword[ii]);
		}
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

			// allocate memory for saving positions
			_safealloc(xx_old,natom,sizeof(double))
			;
			_safealloc(yy_old,natom,sizeof(double))
			;
			_safealloc(zz_old,natom,sizeof(double))
			;

			// calculate the upper limit of volume changes
			prob_vc_upper = prob_cm + prob_vc;

			fclose(fpins);
			return 0;
		} // if keyword found
	} // read through the lines
	fprintf(stderr,"Error: data for HMC simulation not found.\n");
	fprintf(fpouts, "Error: data for HMC simulation not found.\n");
	fclose(fpins);
	exit(1);
}

/**
 * \brief Hybrid Monte Carlo simulation
 */
int hmc()
{
	double ratio;

	// start of the HMC simulation
	// initialize velocities
	fprintf(fpouts, "initializing velocities...\n");
	velinit();

	// if not new run, load from old file
	if (fStart_option!=new_run)
	{
		loadit();
	}

	// calculate total energies
	erfrc();
	rafrc();
	
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
			fnMDmove();
		}
		else if (rndnum[0]<=prob_vc_upper) // volume change moves
		{
			fnVolumeChange();
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

	return 0;
}
