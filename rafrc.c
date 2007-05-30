/* Calculate intra energies and forces
   include bond, angle, dihedral, impropers
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "vars.h"

int bndfrc()
{
    int ii;
    int ii1, ii2;
    float rxij, ryij, rzij;
    float rijsq, rij;
    float delta_r;
    float ubond_temp;
    float fij, fxij, fyij, fzij;
    float De;
    float exp_term, one_minus_exp_term;
    float rmax2;
    // bond stretching calculations
    for (ii=0;ii<nbond;ii++)
    {
	ii1 = bond_idx[ii][0];
	ii2 = bond_idx[ii][1];
	rxij = xx[ii1] - xx[ii2];
	ryij = yy[ii1] - yy[ii2];
	rzij = zz[ii1] - zz[ii2];
	rijsq = rxij*rxij+ryij*ryij+rzij*rzij;
	rij = sqrt(rijsq);

	switch (bond_type[ii])
	{
	    case bond_none: // 0
		break;
	    case bond_harmonic: // 1
		delta_r = rij - Req[ii];
		ubond_temp = 0.5*Kb[ii]*delta_r*delta_r;
		ubond += ubond_temp;
		fij = -Kb[ii]*delta_r/rij;
		fxij = fij*rxij;
		fyij = fij*ryij;
		fzij = fij*rzij;
		// force on atom 1
		fxs[ii1] += fxij;
		fys[ii1] += fyij;
		fzs[ii1] += fzij;
		// force on atom 2
		fxs[ii2] -= fxij;
		fys[ii2] -= fyij;
		fzs[ii2] -= fzij;
		break;
	    case bond_morse: // 2
		delta_r = rij - Req[ii];
		De = Kb[ii];
		exp_term = exp(-alpha[ii]*delta_r);
		one_minus_exp_term = 1.0 - exp_term;
		ubond_temp = De*one_minus_exp_term*one_minus_exp_term;
		ubond += ubond_temp;
		fij = -2.0*De*alpha[ii]*one_minus_exp_term*exp_term/rij;
		fxij = fij*rxij;
		fyij = fij*ryij;
		fzij = fij*rzij;
		// force on atom 1
		fxs[ii1] += fxij;
		fys[ii1] += fyij;
		fzs[ii1] += fzij;
		// force on atom 2
		fxs[ii2] -= fxij;
		fys[ii2] -= fyij;
		fzs[ii2] -= fzij;
		break;
	    case bond_fene: // 3
		rmax2 = Req[ii]*Req[ii];
		ubond_temp = Kb[ii]*rmax2*log(1.0-rijsq/rmax2);
		ubond += ubond_temp;
		fij = -2.0*Kb[ii]*rmax2/(rmax2-rijsq);
		fxij = fij*rxij;
		fyij = fij*ryij;
		fzij = fij*rzij;
		// force on atom 1
		fxs[ii1] += fxij;
		fys[ii1] += fyij;
		fzs[ii1] += fzij;
		// force on atom 2
		fxs[ii2] -= fxij;
		fys[ii2] -= fyij;
		fzs[ii2] -= fzij;
		break;
	    default:
		printf("unknown bond type.\n");
		exit(1);
	} // switch for different bond type
    } // end of bond stretching term calculations
}

int aglfrc()
{
    int ii;
    int iia, iib, iic;
    float xdab, ydab, zdab, xdbc, ydbc, zdbc;
    float rab, rrab, xab, yab, zab;
    float rbc, rrbc, xbc, ybc, zbc;
    float cost, theta, sint;
    float delta_theta;
    float uangle_temp;

    // angle bending calculations
    for (ii=0;ii<nangle;ii++)
    {
	iia = angle_idx[ii][0];
	iib = angle_idx[ii][1];
	iic = angle_idx[ii][2];
	switch (angle_type[ii])
	{
	    case angle_none:
		break;
	    case angle_harmonic:
		// first bond vector



		break;
	    case angle_TRwater:
		break;
	    default:
		printf("unknown angle type.\n");
		exit(1);
	} // ending of switch for different angle types
    } // end of angle bending calculations

}

int rafrc()
{
    int ii;

    ubond = 0.0;
    uangle = 0.0;
    udih = 0.0;
    uimp = 0.0;

    for (ii=0;ii<natom;ii++)
	fxs[ii] = fys[ii] = fzs[ii] = 0.0;

    if (nbond > 0)
	bndfrc();

    if (nangle > 0)
	aglfrc();


}

