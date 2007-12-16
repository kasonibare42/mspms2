#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <ctype.h>
#include "vars.h"

/**
 * \brief Validate the input data.
 */
int fnValidateInput()
{
	int isError;
	double fMinBoxLength;

	isError = 0;

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
		isError = 1;
	}

	if (treq < 0.0)
	{
		fprintf(stderr,"Error: Invalid required temperature. treq = %lf\n",treq);
		fprintf(fpouts, "Error: Invalid required temperature. treq = %lf\n",
				treq);
		isError = 1;
	}

	if (boxlx<=0.0 || boxly<=0.0 || boxlz<=0.0)
	{
		fprintf(stderr,"Error: Invalid box size, boxlx = %lf, boxly = %lf, boxlz = %lf\n",boxlx, boxly, boxlz);
		fprintf(
				fpouts,
				"Error: Invalid box size, boxlx = %lf, boxly = %lf, boxlz = %lf\n",
				boxlx, boxly, boxlz);
		isError = 1;
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
		isError = 1;
	}

	if (sf_type < _NO_SF_POTENTIAL || sf_type > nanotube_my_interp)
	{
		fprintf(stderr,"Error: Invalid solid-fluid interaction type. sf_type = %d\n",sf_type);
		fprintf(fpouts,
				"Error: Invalid solid-fluid interaction type. sf_type = %d\n",
				sf_type);
		isError = 1;
	}

	if (iChargeType < _NO_ELECTROSTATIC_INTERACTION || iChargeType
			> elec_simple_coulomb)
	{
		fprintf(stderr,"Error: Invalid electrostatic interaction type. iChargeType = %d\n", iChargeType);
		fprintf(
				fpouts,
				"Error: Invalid electrostatic interaction type. iChargeType = %d\n",
				iChargeType);
		isError = 1;
	}

	if (isError == 1)
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
	if (prob_cm+prob_vc+prob_id>1.0)
	{
		fprintf( stderr, "Warning: prob_cm (%lf) + prob_vc (%lf) + prob_id (%lf) > 1.0 \n", prob_cm, prob_vc, prob_id);
		fprintf(
				fpouts,
				"Warning: prob_cm (%lf) + prob_vc (%lf) + prob_id (%lf) > 1.0 \n",
				prob_cm, prob_vc, prob_id);
	}
	
	return 0;
}

