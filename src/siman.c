#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "mspms2.h"

int init_siman()
{
	int ii;
	char buffer[LONG_STRING_LENGTH];
	char keyword[100];
	int position_counter;

	fprintf(stderr, "Reading input data for Simulated Annealing...\n");
	fprintf(fpouts, "Reading input data for Simulated Annealing...\n");

	// re-open input file to read extra data section
	fpins = fopen(INPUT,"r");

	while (fgets(buffer, LONG_STRING_LENGTH, fpins)!=NULL)
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

			sscanf(fgets(buffer, LONG_STRING_LENGTH, fpins), "%d %d", &n_tries,
					&iters_fixed_t);
			sscanf(fgets(buffer, LONG_STRING_LENGTH, fpins), "%lf %lf %lf",
					&t_initial, &t_min, &mu_t);

			// set the required temperature for velocity initialization
			t_initial /= epsilon_base;
			treq = t_initial;

			t_min /= epsilon_base;

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

	// if not new run, load from old file
	// FIXME: NOT tested for simulated annealing
	// if (iStart_option!=NEW)
	// {
	// 	loadit();
	// }

	// calculate total energies
	// frclong();
	// frcshort();

	// print out initial values
	// echo();
	// print initial properties
	// printit();
	// make snapshots & movies
	// if (nstep_trj)
	// { 
	// 	trajectory();
	// }

	n_accepts = 0;
	n_rejects = 0;

	T = t_initial;
	done = 0;

	istep = nstep_start;

	while (!done)
	{
		treq = T;
		velinit();
		frclong();
		frcshort();

		for (ii=0; ii<n_tries; ii++)
		{
			vver_nh_3();
			upot = uinter + uintra;
			utot = upot + ukin;
			virial = virial_inter + virial_intra;
			utot = utot + unhts + unhtss;

			if (ii%iters_fixed_t==0)
			{
				fprintf(stdout, "    %10d  %8d %10.4le %10.4le %10.4le %10.4le\n",
				istep, ii, utot*epsilon_base,
				upot*epsilon_base, ukin*epsilon_base, tinst*epsilon_base);
				fprintf(fplog, "    %10d  %8d %10.4le %10.4le %10.4le %10.4le\n",
				istep, ii, utot*epsilon_base,
				upot*epsilon_base, ukin*epsilon_base, tinst*epsilon_base);
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

		counts[11]--;
		// if still in equilibrium run
		// do not do averages
		if (counts[11]>=0)
		{
			continue;
		}

		// accumulators
		collect_aves();

		// do averages
		if ((istep-nstep_eq)%nstep_ave==0)
		{
			averages();
		}

		// apply the cooling schedule to the temperature 
		T /= mu_t;
		unhts = 0.0;
		gg = nfree;
		ss = 0.0;
		ps = 0.0;
		ggs = nfree;
		sss = 0.0;
		pss = 0.0;
		unhtss = 0.0;
	}
	
	snapshot();
	calres();
	ending();

	return 0;
}
