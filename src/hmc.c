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
#include "mspms2.h"

int init_hmc()
{
	int ii;
	const int datalen = 200;
	char buffer[200];
	char keyword[100];
	int position_counter;

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
			
			// printf("prob_id = %lf\n",prob_id);
			// printf("delv = %lf\n",delv);

			// readin insertion/deletion input data if required
			if (prob_id > 0.0)
			{
				fgets(buffer, datalen, fpins);
				position_counter = 0;
				for (ii=0; ii<nspecie; ii++)
				{
					sscanf(&buffer[position_counter], "%lf %lf %lf%n",
							&probability_to_be_selected[ii],
							&probability_to_insert[ii], &fugacity_required[ii],
							&position_counter);
					// initialize zact
					zact[ii] = fugacity_required[ii]/KB_OVER_1E30/treq;
					// printf("fugacity[%d] = %lf\n",ii, fugacity_required[ii]);
					// printf("zact[%d] = %lf\n",ii,zact[ii]);
					// initialize vacancy
					specie_first_vacancy_idx[ii] = -1;
				}
			}

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
	fprintf(stderr,"Error: Data for HMC simulation not found.\n");
	fprintf(fpouts, "Error: Data for HMC simulation not found.\n");
	fclose(fpins);
	exit(1);
}

/**
 * \brief Hybrid Monte Carlo simulation
 */
int hmc()
{
	double ratio;
	void (*pfnRezero)();
	int (*pfnMDtype)();
	
	if (what_ensemble == NVT)
	{ 
		pfnRezero = &rezero_nvt_ts;
		pfnMDtype = &vver_nh_3;
	}
	else if (what_ensemble == NVE)
	{
		pfnRezero = NULL;
		pfnMDtype = &vver;
	}
	else if (what_ensemble == NPT)
	{
		pfnRezero = &rezero_npt_ts;
		pfnMDtype = &npt_respa;
	}
	
	// if not new run, load from old file
	if (fStart_option!=NEW)
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
	if (nstep_trj)
	{ 	
		trajectory();
	}

	// above counts as the first step
	counts[11]--;

	// simulation loop
	for (istep=nstep_start; istep<=nstep; istep++) // NOTE: start from 1 and <=
	{
		ranmar(rndnum, 1);

		if (rndnum[0] <= prob_cm) // canonical moves
		{
			// velocity init and energy calcualtions have been done for the 1st step outside the loop
			if (istep!=nstep_start)
			{
				velinit();
				erfrc();
				rafrc(); // we may not need rafrc here??
			}
			fnMDmove(nstep_md_per_hmc, pfnRezero, pfnMDtype);
		}
		else if (rndnum[0]<=prob_vc_upper) // volume change moves
		{
			fnVolumeChange();
		}
		else // insertions or deletions
		{
			fnInsDelMole();
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

		counts[11]--;
		// if still in equilibrium run
		// adjust the delt, deltv
		// do not do averages
		if (counts[11]>=0)
		{
			if (counts[20]==nstep_delt_adj_cycle) // delt adjustment
			{
				ratio = counts[21]*1.0/nstep_delt_adj_cycle;
				delt = delt*(1.0 - ratio_cm_req + ratio);
				counts[20] = 0;
				counts[21] = 0;
				// delt update
				deltby2 = delt/2.0;
				delts = delt/nstep_inner;
				deltsby2 = delts/2.0;
			}
			if (counts[23]==nstep_delv_adj_cycle) // delv adjustment
			{
				ratio = counts[24]*1.0/nstep_delv_adj_cycle;
				delv = delv*(1.0 - ratio_vc_req + ratio);
				counts[23] = 0;
				counts[24] = 0;
			}
			continue;
		}

		// accumulators
		collect_aves();

		// do averages
		if ((istep-nstep_eq)%nstep_ave==0)
		{
			averages();
		}
	}

	return 0;
}