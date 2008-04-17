#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <ctype.h>
#include "mspms2.h"

/**
 * \brief Validate the input data.
 */
int fnValidateInput()
{
	bool isError;
	double fMinBoxLength;

	isError = false;

	if (ij<0 || ij>31328 || jk<0 || jk>30081)
	{
		fprintf(stderr, "Error: The first random number seed must have a value between 0 and 31328. ij = %d\n",ij);
		fprintf(
		stderr,
		"Error: The second seed must have a value between 0 and 30081. jk = %d\n",
		jk);
		fprintf(
				fpouts,
				"Error: The first random number seed must have a value between 0 and 31328. ij = %d\n",
				ij);
		fprintf(
				fpouts,
				"Error: The second seed must have a value between 0 and 30081. jk = %d\n",
				jk);
		isError = true;
	}
	
	if (what_simulation<MOLECULAR_DYNAMICS || what_simulation>SIMULATED_ANNEALING)
	{
		fprintf(stderr,"Error: unknown simulation type.\n");
		fprintf(fpouts, "Error: unknown simulation type.\n");
		isError = true;
	}

	if (treq < 0.0)
	{
		fprintf(stderr,"Error: Invalid required temperature. treq = %lf\n",treq);
		fprintf(fpouts, "Error: Invalid required temperature. treq = %lf\n",
				treq);
		isError = true;
	}

	if (boxlx<=0.0 || boxly<=0.0 || boxlz<=0.0)
	{
		fprintf(stderr,"Error: Invalid box size, boxlx = %lf, boxly = %lf, boxlz = %lf\n",boxlx, boxly, boxlz);
		fprintf(
				fpouts,
				"Error: Invalid box size, boxlx = %lf, boxly = %lf, boxlz = %lf\n",
				boxlx, boxly, boxlz);
		isError = true;
	}

	fMinBoxLength = boxlx;
	fMinBoxLength = fMinBoxLength<boxly ? fMinBoxLength : boxly;
	fMinBoxLength = fMinBoxLength<boxlz ? fMinBoxLength : boxlz;
	fMinBoxLength /= 2.0;
	if (rcuton>fMinBoxLength || rcutoff>fMinBoxLength || rcutoffelec
			>fMinBoxLength)
	{
		fprintf(stderr,"Error: Invalid cutoff. rcuton = %lf, rcutoff = %lf, rcutoffelec = %lf\n", rcuton, rcutoff, rcutoffelec);
		fprintf(
				fpouts,
				"Error: Invalid cutoff. rcuton = %lf, rcutoff = %lf, rcutoffelec = %lf\n",
				rcuton, rcutoff, rcutoffelec);
		isError = true;
	}

	if (iSF_type < SF_NONE || iSF_type > SF_NANOTUBE_MY_INTERP)
	{
		fprintf(stderr,"Error: Invalid solid-fluid interaction type. iSF_type = %d\n",iSF_type);
		fprintf(fpouts,
				"Error: Invalid solid-fluid interaction type. iSF_type = %d\n",
				iSF_type);
		isError = true;
	}

	if (iChargeType < ELECTROSTATIC_NONE || iChargeType > ELECTROSTATIC_SIMPLE_COULOMB)
	{
		fprintf(stderr,"Error: Invalid electrostatic interaction type. iChargeType = %d\n", iChargeType);
		fprintf(
				fpouts,
				"Error: Invalid electrostatic interaction type. iChargeType = %d\n",
				iChargeType);
		isError = true;
	}

	if (isError == true)
	{
		exit(1);
	}

	return 0;
}

/**
 * \brief Validate the initialized data.
 */
int fnValidateInit()
{
	if (pdisp+pvolm+pmake+pkill>1.0)
	{
		fprintf( stderr, "Warning: pdisp (%lf) + pvolm (%lf) + pmake (%lf) + pkill (%lf) > 1.0 \n", pdisp, pvolm, pmake, pkill);
		fprintf( fpouts, "Warning: pdisp (%lf) + pvolm (%lf) + pmake (%lf) + pkill (%lf) > 1.0 \n", pdisp, pvolm, pmake, pkill);
	}
	
	return 0;
}

