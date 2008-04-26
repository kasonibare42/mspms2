/**
 * Project: mspms2
 * File: hmc.c
 * 
 * Copyright (C) 2008    Yang Wang <ywangd@gmail.com>
 * Created @ Dec 07, 2007
 * 
 * Description:
 *   Functions related to HMC operation go here.
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
	ranmar(rndnum, 1);

	if (rndnum[0] <= pdisp) // canonical moves
	{
		// velocity init and energy calcualtions have been done for the 1st step outside the loop
		if (istep!=nstep_start)
		{
			velinit();
			frclong();
			frcshort();
		}
		md_move(nstep_md_per_hmc);
	}
	else if (rndnum[0]<=pvolm_upper) // volume change moves
	{
		volume_change();
	}
	else // insertions or deletions
	{
		fnInsDelMole();
	}

	if (bEquilibrium==true) // Only adjust delt and delv while it is still in equilibrium
	{
		if (counts[20]==nstep_delt_adj_cycle) // delt adjustment
		{
			rinst_disp = counts[21]*1.0/nstep_delt_adj_cycle;
			delt = delt*(1.0 - rreq_disp + rinst_disp);
			counts[20] = 0;
			counts[21] = 0;
			// delt update
			deltby2 = delt/2.0;
			delts = delt/nstep_inner;
			deltsby2 = delts/2.0;
		}
		if (counts[23]==nstep_delv_adj_cycle) // delv adjustment
		{
			rinst_volm = counts[24]*1.0/nstep_delv_adj_cycle;
			delv = delv*(1.0 - rreq_volm + rinst_volm);
			counts[23] = 0;
			counts[24] = 0;
		}
	} // End of equilibrium check

	return 0;
}
