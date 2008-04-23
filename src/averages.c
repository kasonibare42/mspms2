/**
 * Project: mspms2
 * File: averages.c
 * 
 * Copyright (C) 2008    Yang Wang <ywangd@gmail.com>
 * Created @ 15/04/2008
 * 
 * Description:
 * 		Calculate the accums and averages.  
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <ctype.h>
#include "mspms2.h"


int collect_aves()
{
	// Calculate energies and pressures
	upot = uinter + uintra;
	utot = upot + ukin;
	virial = virial_inter + virial_intra;
	// Add energy of thermostat. 
	// If nose hoover is not used, they will just be zero
	utot = utot + unhts + unhtss + utsbs;
	// calculate ideal pressure part, rho*K*T = rho(*)*T(*)
	pideal=natom/(boxlx*boxly*boxlz)*tinst;
	// Do not need to recalculate lrc here, it should be calculated
	// elsewhere when variables changed.
	pinst = pideal + (virial_inter+virial_intra)/(boxlx*boxly*boxlz);
	// Add long range corrections into total energy and pressure if needed.
	if (isLJlrcOn)
	{
		utot += uljlrc;
		pinst += pljlrc;
	}

	// accumulators
	accum[0][0] += utot;
	accum[0][1] += utot*utot;
	accum[1][0] += upot;
	accum[1][1] += upot*upot;
	accum[2][0] += ukin;
	accum[2][1] += ukin*upot;
	accum[3][0] += uinter;
	accum[3][1] += uinter*uinter;
	accum[4][0] += uintra;
	accum[4][1] += uintra*uintra;
	accum[5][0] += uvdw;
	accum[5][1] += uvdw*uvdw;
	accum[6][0] += ubond;
	accum[6][1] += ubond*ubond;
	accum[7][0] += uangle;
	accum[7][1] += uangle*uangle;
	accum[8][0] += udih;
	accum[8][1] += udih*udih;
	accum[9][0] += uimp;
	accum[9][1] += uimp*uimp;
	accum[10][0] += uewald;
	accum[10][1] += uewald*uewald;
	accum[11][0] += ureal;
	accum[11][1] += ureal*ureal;
	accum[12][0] += ufourier;
	accum[12][1] += ufourier*ufourier;
	accum[13][0] += uself;
	accum[13][1] += uself*uself;
	accum[14][0] += usflj;
	accum[14][1] += usflj*usflj;
	accum[15][0] += tinst;
	accum[15][1] += tinst*tinst;
	accum[16][0] += uvacuum;
	accum[16][1] += uvacuum*uvacuum;
	accum[17][0] += uwolf;
	accum[17][1] += uwolf*uwolf;
	accum[18][0] += pinst;
	accum[18][1] += pinst*pinst;
	accum[19][0] += boxv;
	accum[19][1] += boxv*boxv;
	accum[20][0] += pideal;
	accum[20][1] += pideal*pideal;

	return 0;
}


int averages()
{
	int ii;
	double temp1;

	for (ii=0; ii<NCOUNTS_MAX; ii++)
	{
		temp1 = accum[ii][0]/nstep_ave;
		accum[ii][2] += temp1;
		accum[ii][3] += accum[ii][1]/nstep_ave;
		accum[ii][4] += temp1*temp1;
		// rezero
		accum[ii][0] = 0.0;
		accum[ii][1] = 0.0;
	}
	counts[10]++; // number of average cycles

	return 0;
}

int calres()
{
	int ii;
	double ave_of_square, ave_of_ave_square;
	double ave, err, fluc;
	for (ii=0; ii<NCOUNTS_MAX; ii++)
	{
		ave = accum[ii][5] = accum[ii][2]/counts[10]; // ave
		ave_of_square = accum[ii][3]/counts[10];
		ave_of_ave_square = accum[ii][4]/counts[10];
		err = accum[ii][6] = sqrt(fabs(ave_of_ave_square-ave*ave)); // err
		fluc = accum[ii][7] = sqrt(fabs(ave_of_square-ave*ave)); // fluc
	}

	return 0;
}

int rezero()
{
	int ii, jj;
	for (ii=0;ii<NCOUNTS_MAX;ii++)
	{
		counts[ii] = 0;
		for (jj=0;jj<5;jj++)
		{
			accum[ii][jj] = 0.0;
		}
	}
	
	return 0;
}

