/*
 * Maintainable Simplex Purpose Molecular Simulator 2
 * Rewritten with standard C language
 * The goal is simpler, quicker, easier, better.
 * No class and other complex data structures.
 * Use common names for input files. So the every job must
 * have its own working directory.
 * 
 * NOTE:
 *      112607: Test move from CVS to SVN
 *      112707: Successfully imported the SVN repository to Google Code host
 * 
 * 
 * Written by Yang Wang 2007
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
	upot = uinter + uintra;
	utot = upot + ukin;
	virial = virial_inter + virial_intra;
	// add energy of thermostat, if nose hoover is not used, they will just be zero
	utot = utot + unhts + unhtss + utsbs;
	// calculate ideal pressure part, rho*K*T = rho(*)*T(*)
	pideal=natom/(boxlx*boxly*boxlz)*tinst;
	// do not need to recalculate lrc here, it should be calculated
	// elsewhere when variables changed
	pinst = pideal + (virial_inter+virial_intra)*VIRIAL_TO_PRESSURE/(boxlx
			*boxly*boxlz);
	// add long range corrections into total energy and pressure if needed
	if (isLJlrcOn)
	{
		utot += uljlrc;
		pinst += pljlrc;
	}
	fprintf(stdout,"%10d %10.4le %10.4le %10.4le %10.4le %10.4le %10.4le\n",
	istep,utot,upot,ukin,tinst,pinst,boxv);

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
	}
	else if (what_simulation == SIMULATED_ANNEALING)
	{
		pfnSimulation = &siman;
		pfnSimulation();
		snapshot();
		calres();
		ending();
		exit(1);
	}
	
	if (pfnSimulation == NULL)
	{
		fprintf(stderr,"Error: Cannot determine the simulation type.");
		exit(1);
	}
	exit(1);
	
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
			// adjust delt?
		}
		if (bEquilibrium==true)
		{
			nstep_start = 1;
			nstep_end = nstep;
			calres();
			ending();
			rezero();
			bEquilibrium = false;
		}
		else
		{
			run = false;
		}
	}

	/// Take the last snapshot of the system.
	snapshot();
	/// Calculate the final results.
	calres();
	/// Write out results and clean up.
	ending();

	return 0;
}

