#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "vars.h"
#include "funcs.h"

/* set up parameters for this simulated annealing run */

/* how many points do we try before stepping (move to the next temperature) */
double n_tries;
/* how many iterations for each T? (how many MD steps before use the acceptance critirea) */
int iters_fixed_t;
/* initial temperature */
double t_initial;
/* damping factor for temperature */
double mu_t, t_min;

int init_siman()
{
	int ii;
	const int datalen = 200;
	char buffer[200];
	char keyword[100];
	int position_counter;

	fprintf(stderr, "Reading input data for Simulated Annealing...\n");
	fprintf(fpouts, "Reading input data for Simulated Annealing...\n");

	// re-open input file to read extra data section
	fpins = fopen(INPUT,"r");

	while (fgets(buffer, datalen, fpins)!=NULL)
	{
		sscanf(buffer, "%s", keyword);
		for (ii=0; ii<strlen(keyword); ii++)
		{
			keyword[ii] = toupper(keyword[ii]);
		}
		if (!strcmp(keyword, "SIMAN"))
		{
			fprintf(stderr,"Data section for Simulated Annealing found...\n");
			fprintf(fpouts, "Data section for Simulated Annealing found...\n");

			sscanf(fgets(buffer, datalen, fpins), "%d %d", &n_tries,
					&iters_fixed_t);
			sscanf(fgets(buffer, datalen, fpins), "%lf %lf %lf", &t_initial,
					&t_min, &mu_t);

			// set the required temperature for velocity initialization
			treq = t_initial;

			// allocate memory for saving positions
			_safealloc(xx_old,natom,sizeof(double))
			;
			_safealloc(yy_old,natom,sizeof(double))
			;
			_safealloc(zz_old,natom,sizeof(double))
			;

			fclose(fpins);
			return 0;
		} // if keyword found
	} // read through the lines
	fprintf(stderr,"Error: Data for Simulated Annealing not found.\n");
	fprintf(fpouts, "Error: Data for Simulated Annealing not found.\n");
	fclose(fpins);
	exit(1);
}

int siman()
{
	int ii, done;
	double T;
	int n_accepts, n_rejects;

	istep = 0;

	// if not new run, load from old file
	// FIXME: NOT tested for simulated annealing
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

	n_accepts = 0;
	n_rejects = 0;

	T = t_initial;
	done = 0;

	istep = nstep_start;

	while (!done)
	{
		for (ii = 0; ii < n_tries; ii++)
		{
			if (istep != nstep_start)
			{
				treq = T;
				velinit();
				erfrc();
				rafrc(); // we may not need rafrc here??
				// rezero nvt related variables
				unhts = 0.0;
				ss = 0.0;
				ps = 0.0;
				sss = 0.0;
				pss = 0.0;
				unhtss = 0.0;
			}
			if (fnMDmove(iters_fixed_t, &vver_nh_3) == 1) // if accepted
			{
				++n_accepts;
			}
			else
			{
				++n_rejects;
			}

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
		
		istep++;
		if (T < t_min)
		{
			done = 1;
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
		
		/* apply the cooling schedule to the temperature */
		T /= mu_t;
	}

	fprintf(stderr, "accepted = %d     rejected = %d\n", n_accepts, n_rejects);

	return 0;
}
