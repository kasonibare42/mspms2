#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <ctype.h>
#include "vars.h"

int fnVolumeChange()
{
	// zero the delta energies
	fDeltaUljlrc = 0.0;
	fDeltaPljlrc = 0.0;
	fDeltaU = 0.0;

	// save the old states
	uinter_old = uinter;
	uvdw_old = uvdw;
	uewald_old = uewald;
	uwolf_old = uwolf;
	ucoulomb_old = ucoulomb;
	usflj_old = usflj;
	virial_inter_old = virial_inter;
	boxlx_old = boxlx;
	boxly_old = boxly;
	boxlz_old = boxlz;
	boxv_old = boxv;
	uljlrc_old = uljlrc;
	pljlrc_old = pljlrc;

	// Random volume change
	ranmar(rndnum, 1);
	// calculate the new box size
	fVolumeNew = boxv + 2.0*delv*(rndnum[0]-0.5);
	fRatioNewV2OldV = fVolumeNew/boxv;
	fRatioNewL2OldL = pow(fRatioNewV2OldV, 0.333333333);
	fLengthNew = boxlx*fRatioNewL2OldL;
	fWidthNew = boxly*fRatioNewL2OldL;
	fHeightNew = boxlz*fRatioNewL2OldL;

	/** 
	 * Check whether the the cutoff is still valid. If not, exit.
	 */
	fMinHalf = fLengthNew;
	if (fMinHalf > fWidthNew)
	{
		fMinHalf = fWidthNew;
	}
	if (fMinHalf > fHeightNew)
	{
		fMinHalf = fHeightNew;
	}
	fMinHalf *= 0.5;
	if (rcuton>fMinHalf || rcutoff>fMinHalf || rcutoffelec>fMinHalf)
	{
		fprintf(
		stderr,
		"Error: Cutoff is larger than the half of the minimal box size. fMinHalf = %lf\n",
		fMinHalf);
		fprintf(
				fpouts,
				"Error: Cutoff is larger than the half of the minimal box size. fMinHalf = %lf\n",
				fMinHalf);
		exit(1);
	}

	/// Calculate new LJ lrc is LJlrc is enabled
	if (isLJlrc)
	{
		// delv_ljlrc(pBox, del_uljlrc, del_pljlrc, delv);
		fDeltaUljlrc = uljlrc - uljlrc_old;
		fDeltaU += fDeltaUljlrc;
	}

	// Set the box size to the new data for energy calculation
	boxv = fVolumeNew;
	boxlx = fLengthNew;
	boxly = fWidthNew;
	boxlz = fHeightNew;

	// calculate the total energy, only inter-molecular energy
	erfrc();

	// calculate the energy difference
	fDeltaU += (uinter - uinter_old);

	// calculate Hamotonial difference
	dH = fDeltaU*rRgas/treq + preq*delv/kb_1e30/treq;
	dH = dH - nmole*log(fRatioNewV2OldV);

	// attempted vc moves
	icounter[23]++;
	
	// check if the move is accepted
	isAccept = 0;
	if (dH <= 0.0)
	{
		isAccept = 1;
	}
	else
	{
		ranmar(rndnum, 1);
		if (rndnum[0] < exp(-dH))
		{
			isAccept = 1;
		}
	}
	if (isAccept == 1)
	{
		// accepted vc moves
		icounter[24]++;
		
		// re-calculate box size related variables
		Vfactor_ewald = 2.0*pi/(boxlx*boxly*boxlz);
		TWOPI_LX = 2.0*pi/boxlx;
		TWOPI_LY = 2.0*pi/boxly;
		TWOPI_LZ = 2.0*pi/boxlz;
		// 1D ewald constant
		twopi_over_3v = 2.0*pi/3.0/boxlx/boxly/boxlz;
	}
	else
	{
		uinter = uinter_old;
		uvdw = uvdw_old;
		uewald = uewald_old;
		uwolf = uwolf_old;
		ucoulomb = ucoulomb_old;
		usflj = usflj_old;
		virial_inter = virial_inter_old;
		boxlx = boxlx_old;
		boxly = boxly_old;
		boxlz = boxlz_old;
		boxv = boxv_old;
		uljlrc = uljlrc_old;
		pljlrc = pljlrc_old;
	}

	return 0;
}
