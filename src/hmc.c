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

/**
 * \brief Hybrid Monte Carlo simulation
 */
int hmc()
{
	static double ratio;
	static void (*pfnRezero)();
	static int (*pfnMDtype)();

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

	ranmar(rndnum, 1);

	if (rndnum[0] <= pdisp) // canonical moves
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
	else if (rndnum[0]<=pvolm_upper) // volume change moves
	{
		fnVolumeChange();
	}
	else // insertions or deletions
	{
		fnInsDelMole();
	}

	if (bEquilibrium==true)
	{
		if (counts[20]==nstep_delt_adj_cycle) // delt adjustment
		{
			ratio = counts[21]*1.0/nstep_delt_adj_cycle;
			delt = delt*(1.0 - rreq_disp + ratio);
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
			delv = delv*(1.0 - rreq_volm + ratio);
			counts[23] = 0;
			counts[24] = 0;
		}
	}

	return 0;
}
