#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_siman.h>
#include "vars.h"
#include "funcs.h"

/* set up parameters for this simulated annealing run */

/* how many points do we try before stepping */
double n_tries;
/* how many iterations for each T? */
int iters_fixed_t;
/* initial temperature */
double t_initial;
/* damping factor for temperature */
double mu_t, t_min;

// The x array used for simulated annealing
double *siman_x;

int init_siman()
{
	int iSize;
	iSize = natom*sizeof(double);

	siman_x = (double *) malloc(3*iSize); // factor of 3 for x,y,z
	// copy the x, y, z coordinate values
	memcpy(siman_x, xx, iSize);
	memcpy(siman_x+natom, yy, iSize);
	memcpy(siman_x+2*natom, zz, iSize);

	return 0;
}

// Add a initialize function for initialize the *xp needed for the calculations

/* now some functions to test in one dimension */
double target_energy_func(void *xp)
{
	return utot;
}

// This function is not used in the current version of gsl_siman
double distance_func(void *xp, void *yp)
{

	return 0.0;
}

// This function is used to generate the new value set for x
void take_step_func(const gsl_rng * r, void *xp, double step_size)
{
	// We do not need step_size 
	// We use MD code to generate the new value set for x (i.e. x, y, z)
	// Use the initial values of xp
	int ii;
	int iSize;
	double *ptr_x, *ptr_y, *ptr_z;
	double target_energy;

	// set the initial x,y,z values based on the passed parameter xp
	iSize = natom*sizeof(double);

	ptr_x = (double *)xp;
	ptr_y = ptr_x + natom;
	ptr_z = ptr_y + natom;

	memcpy(xx, ptr_x, iSize);
	memcpy(yy, ptr_y, iSize);
	memcpy(zz, ptr_z, iSize);

	// Now we can start the MD moves
	vver();

	// potential energy
	upot = uinter + uintra;
	// total energy
	utot = upot + ukin;
	// add energy of thermostat, if nose hoover is not used, they will just be zero
	utot = utot + unhts + unhtss + utsbs;
	// add long range corrections into total energy and pressure if needed
	if (isLJlrcOn)
	{
		utot += uljlrc;
	}

	// printit();

	// copy the values back to the parameter xp as the new value set
	memcpy(ptr_x, xx, iSize);
	memcpy(ptr_y, yy, iSize);
	memcpy(ptr_z, zz, iSize);
}

void print_func(void *xp)
{
	printit();
}

int siman()
{
	gsl_siman_params_t params;

	init_siman();

	n_tries = 20;
	iters_fixed_t = 1000;
	t_initial = 1000;
	mu_t = 1.01;
	t_min = 50;

	params.n_tries = n_tries;
	params.iters_fixed_T = iters_fixed_t;
	params.k = Rgas;
	params.t_initial = t_initial;
	params.mu_t = mu_t;
	params.t_min = t_min;
	/* max step size in random walk for simulated annealing */
	// It is not used for the 
	params.step_size = 0.0;

	const gsl_rng_type * T;
	gsl_rng * r;

	gsl_rng_env_setup();

	T = gsl_rng_default;
	r = gsl_rng_alloc(T);

	gsl_siman_solve(r, siman_x, target_energy_func, take_step_func,
			distance_func, print_func, NULL, NULL, NULL,
			3*natom*sizeof(double), params);

	gsl_rng_free(r);

	return 0;
}
