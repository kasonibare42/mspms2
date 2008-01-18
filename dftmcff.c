/**
 * \brief Calculate the energy of metal clusters (DFT fit)
 * 
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_deriv.h>
#include "vars.h"

#define X_AXIS	0
#define Y_AXIS	1
#define Z_AXIS	2

extern void ffieldcu_(int* nm, int* ndata, double* x, double* y, double* z,
		double* ff);

typedef struct _dftmcffparam
{
	int iWhichAxis;
	int index;
} DFTMCFFPARAM;

double deriv_ffieldcu(double pos, void* params)
{
	int ndata;
	double posold;
	double energy;

	DFTMCFFPARAM *p = (DFTMCFFPARAM *)params;
	if (p->iWhichAxis == X_AXIS)
	{
		posold = xx[p->index];
		xx[p->index] = pos;
		ffieldcu_(&natom, &ndata, xx, yy, zz, &energy);
		xx[p->index] = posold;
	}
	else if (p->iWhichAxis == Y_AXIS)
	{
		posold = yy[p->index];
		yy[p->index] = pos;
		ffieldcu_(&natom, &ndata, xx, yy, zz, &energy);
		yy[p->index] = posold;
	}
	else // Z_AXIS
	{
		posold = zz[p->index];
		zz[p->index] = pos;
		ffieldcu_(&natom, &ndata, xx, yy, zz, &energy);
		zz[p->index] = posold;
	}

	return energy*EV_TO_J_PER_MOLE;
}

int fnMetalClusterFF()
{
	int ii;
	int ndata; // no use at all, just for consistency with the FORTRAN code
	double energy;

	// energy calculation
	ffieldcu_(&natom, &ndata, xx, yy, zz, &energy);

	gUMetalClusterSession = energy*EV_TO_J_PER_MOLE; // conver to J/mol

	// numerical forces
	DFTMCFFPARAM param;
	gsl_function FF;
	FF.function = &deriv_ffieldcu;
	FF.params = &param;
	double value, abserr;

	for (ii=0; ii<natom; ii++)
	{
		param.index = ii;

		param.iWhichAxis = X_AXIS;
		gsl_deriv_central(&FF, xx[ii], 1.0e-8, &value, &abserr);
		fxl[ii] += value;
		// printf("%lf   ", value);

		param.iWhichAxis = Y_AXIS;
		gsl_deriv_central(&FF, yy[ii], 1.0e-8, &value, &abserr);
		fyl[ii] += value;
		// printf("%lf   ", value);

		param.iWhichAxis = Z_AXIS;
		gsl_deriv_central(&FF, zz[ii], 1.0e-8, &value, &abserr);
		fzl[ii] += value;
		// printf("%lf   \n", value);
	}

	return 0;
}
