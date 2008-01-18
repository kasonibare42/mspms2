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

extern void ffieldcu_single__(int* nm, int* ndata, double* x, double* y, double* z,
		double* ff);


typedef struct _dftmcffparam
{
	int iWhichAxis;
	int index;
} DFTMCFFPARAM;

double deriv_ffieldcu(double pos, void* params)
{
	int ndata, indexF;
	double posold;
	double energy;

	DFTMCFFPARAM *p = (DFTMCFFPARAM *)params;
	
	ndata = p->index;
	indexF = ndata + 1; // Fortran index
	
	if (p->iWhichAxis == X_AXIS)
	{
		posold = xx[ndata];
		xx[ndata] = pos;
		ffieldcu_single__(&natom, &indexF, xx, yy, zz, &energy);
		xx[ndata] = posold;
	}
	else if (p->iWhichAxis == Y_AXIS)
	{
		posold = yy[ndata];
		yy[ndata] = pos;
		ffieldcu_single__(&natom, &indexF, xx, yy, zz, &energy);
		yy[ndata] = posold;
	}
	else // Z_AXIS
	{
		posold = zz[ndata];
		zz[ndata] = pos;
		ffieldcu_single__(&natom, &indexF, xx, yy, zz, &energy);
		zz[ndata] = posold;
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

	/*
	printf("energy = %lf\n",energy);
	double tmp = 0.0;
	int index;
	for (ii=0;ii<natom;ii++)
	{
		index = ii+1;
		ffieldcu_single__(&natom, &index, xx, yy, zz, &energy);
		tmp += energy;
	}
	tmp = tmp/natom/2.0;
	printf("tmp = %lf\n",tmp);
	exit(1);
	*/
	
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
		fxl[ii] -= value;
		// printf("%lf   ", value);

		param.iWhichAxis = Y_AXIS;
		gsl_deriv_central(&FF, yy[ii], 1.0e-8, &value, &abserr);
		fyl[ii] -= value;
		// printf("%lf   ", value);

		param.iWhichAxis = Z_AXIS;
		gsl_deriv_central(&FF, zz[ii], 1.0e-8, &value, &abserr);
		fzl[ii] -= value;
		// printf("%lf   \n", value);
	}

	return 0;
}
