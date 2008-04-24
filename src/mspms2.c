/**
 * Project: mspms2
 * File: mspms2.c
 * 
 * Copyright (C) 2008    Yang Wang <ywangd@gmail.com>
 * Created @ 2007
 * Modified @ Apr 24, 2008
 * 
 * Description:
 *   Maintainable Simplex Purpose Molecular Simulator 2
 *   Rewritten with standard C language.
 *   The goal is simpler, quicker, easier, better.
 *   No class and other complex data structures.
 *   Use common names for input files. So the every job must
 *   have its own working directory.
 * 
 * NOTE:
 *      112607: Test move from CVS to SVN
 *      112707: Successfully imported the SVN repository to Google Code host
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

int printit()
{
	fprintf(stdout,"%10d %10.4le %10.4le %10.4le %10.4le %10.4le %10.4le\n",
	istep,utot*epsilon_base*RGAS,upot*epsilon_base*RGAS,ukin*epsilon_base*RGAS,
	tinst*epsilon_base,pinst*pressure_base*1.0e5,boxv*sigma_base*sigma_base*sigma_base);

	fprintf(
			fplog,
			"%10d %10.4le %10.4le %10.4le %10.4le %10.4le %10.4le %10.4le %10.4le %10.4le ",
			istep, utot, upot, ukin, tinst, uinter, uintra, uvdw, ubond, uangle);

	fprintf(
			fplog,
			"%10.4le %10.4le %10.4le %10.4le %10.4le %10.4le %10.4le %10.4le %10.4le ",
			udih, uimp, uewald, usflj, unhts, unhtss, virial, virial_inter,
			virial_intra);

	fprintf(fplog, "%10.4le %10.4le %10.4le\n", utsbs, pinst, boxv);

	return 0;
}


/**
 * \brief The entrance and main loop of the program.
 */
int main(int argc, char *argv[])
{
	bool run;
	int (*pfnSimulation)(); // the function pointer to different type of simulations
	
	/// Open output file and print text graph.
	opening();
	/// Read in.mspms, cfg.mspms, and xyz.mspms.
	readins();
	/// Validate the input data.
	fnValidateInput();
	/// Initialize variables according to the input data.
	init_vars();
	/// Validate the initialized variables.
	fnValidateInit();
	// Calculate energy and forces.
	frclong();
	frcshort();
	// Echo the input
	echo();
	// Print out the results from the zero step.
	printit();
	
	/*
	int ii;
	for (ii=0;ii<natom;ii++)
	{
		printf("%d: fxl=%lf, fyl=%lf, fzl=%lf | fxs=%lf, fys=%lf, fzs=%lf\n",
				ii, fxl[ii]*force_base, fyl[ii]*force_base, fzl[ii]*force_base, 
				fxs[ii]*force_base, fys[ii]*force_base, fzs[ii]*force_base);
	}
	for (ii=0;ii<natom;ii++)
	{
		printf("%d, vx=%lf, vy=%lf, vz=%lf\n", 
				ii, vx[ii]*velocity_base, vy[ii]*velocity_base, vz[ii]*velocity_base);
	}
	*/
	
	if (nstep_trj!=0)
	{
		trajectory();
	}
	
	pfnSimulation = NULL;
	/// Decide what type of simulation to run. 
	if (what_simulation == MOLECULAR_DYNAMICS)
	{
		if (what_ensemble == NVT)
		{
			pfnSimulation = &vver_nh_3;
		}
		else if (what_ensemble == NPT)
		{
			pfnSimulation = &npt_respa;
		}
		else if (what_ensemble == NVE)
		{
			pfnSimulation = &vver;
		}
	}
	else if (what_simulation == HYBRID_MONTE_CARLO)
	{
		pfnSimulation = &hmc;
	} // End of simulation functtion assignment 
	
	if (pfnSimulation == NULL)
	{
		fprintf(stderr,"Error: Unknown simulation type.\n");
		exit(1);
	}
	
	run = true;
	while (run)
	{
		for (istep=nstep_start;istep<=nstep_end;istep++)
		{
			pfnSimulation(); // one simulation step
			
			collect_aves();

			// cm_positions();
			
			// do averages
			if (istep%nstep_ave==0)
			{
				averages();
			}
			// print out, snapshot, trajectory, save
			if (istep%nstep_print==0)
			{
				printit();
			}
			// snapshot
			if (nstep_ss>0 && istep%nstep_ss==0)
			{
				snapshot();
			}
			// trajectory
			if (nstep_trj>0 && istep%nstep_trj==0)
			{
				trajectory();
			}
			// save
			if (istep%nstep_save==0)
			{
				saveit();
			}
		} // End of simulation step loop
		if (bEquilibrium==true) // Reset for data taking run
		{
			nstep_start = 1;
			nstep_end = nstep;
			calres();
			ending();
			rezero();
			bEquilibrium = false;
		}
		else // Stop the simulation after data taking run
		{
			run = false;
		}
	} // Simulation block

	/// Take the last snapshot of the system.
	snapshot();
	/// Calculate the final results.
	calres();
	/// Write out results and clean up.
	ending();

	return 0;
}

