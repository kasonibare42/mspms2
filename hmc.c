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

int hmc()
{
	int ii;
	double prob_vc_upper;
	double *xx_old, *yy_old, *zz_old;
	double upot_old, utot_old;
	double upot_new, utot_new;
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
	velinit();

	// if not new run, load from old file
	if (fStart_option!=new_run)
		loadit();

	// calculate total energies
	erfrc();
	rafrc();

	// save old positions and energies
	for (ii=0; ii<natom; ii++)
	{
		xx_old[ii] = xx[ii];
		yy_old[ii] = yy[ii];
		zz_old[ii] = zz[ii];
	}
	upot_old = uinter + uintra;
	utot_old = upot_old + ukin;

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
				// save old positions and energies
				for (ii=0; ii<natom; ii++)
				{
					xx_old[ii] = xx[ii];
					yy_old[ii] = yy[ii];
					zz_old[ii] = zz[ii];
				}
				upot_old = uinter + uintra;
				utot_old = upot_old + ukin;
			}
			for (ii=0; ii<nstep_md_per_hmc; ii++)
			{
				vver();
			}
			upot_new = uinter + uintra;
			utot_new = upot_new + ukin;

			dH = (utot_new - utot_old)*rRgas/treq;
			isAccept = 0;
			if (dH <= 0.0)
				isAccept = 1;
			else
			{
				ranmar(rndnum, 1);
				if (rndnum[0] < exp(-dH))
					isAccept = 1;
			}

			if (isAccept == 1)
			{
				icounter[20]++; // accepted canonical moves
			}
			else
			{
				xx[ii] = xx_old[ii];
				yy[ii] = yy_old[ii];
				zz[ii] = zz_old[ii];
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

		if ((istep-nstep_eq)%nstep_delt_adj_cycle==0) // delt adjustment
		{
			ratio = icounter[20]*1.0/nstep_delt_adj_cycle;
			delt *= (1.0 - ratio_cm_req + ratio);
			icounter[20] = 0;
			// delt update
			deltby2 = delt/2.0;
			delts = delt/nstep_inner;
			deltsby2 = delts/2.0;
		}
		if ((istep-nstep_eq)%nstep_delv_adj_cycle==0) // delv adjustment
		{
			ratio = icounter[21]*1.0/nstep_delv_adj_cycle;
			delv = delv*(1.0 - ratio_vc_req + ratio);
			icounter[21] = 0;
			// length update
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
		// do not do averages
		if (icounter[11]>=0)
		{
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
