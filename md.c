/*
 * Functions related to MD operation go here
 *
 * Written by Yang Wang 11-29-2007
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <ctype.h>
#include "vars.h"

extern int erfrc();
extern int rafrc();
extern int vver_nh_3();
extern int npt_respa();
extern int averages();
extern int loadit();
extern int saveit();
extern int printit();
extern int snapshot();
extern int trajectory();
extern int velinit();
extern int echo();

int vver() // velocity verlet
{
    int ii, ll;

    for (ii=0;ii<natom;ii++)
    {
		// the factor of 1.0e-5 is based on Angstrom (from the force)
		// and femto second (from delt)
		vx[ii] += (deltby2*fxl[ii]*1.0e-5/aw[ii]);
		vy[ii] += (deltby2*fyl[ii]*1.0e-5/aw[ii]);
		vz[ii] += (deltby2*fzl[ii]*1.0e-5/aw[ii]);
    }

    for (ll=0;ll<nstep_inner;ll++)
    {
		for (ii=0;ii<natom;ii++)
		{
		    vx[ii] += (deltsby2*fxs[ii]*1.0e-5/aw[ii]);
		    vy[ii] += (deltsby2*fys[ii]*1.0e-5/aw[ii]);
		    vz[ii] += (deltsby2*fzs[ii]*1.0e-5/aw[ii]);
		    xx[ii] = xx[ii] + delts*vx[ii]*1.0e-5;
		    yy[ii] = yy[ii] + delts*vy[ii]*1.0e-5;
		    zz[ii] = zz[ii] + delts*vz[ii]*1.0e-5;
		}
	
		// intra forces, short ranged
		rafrc();
	
		// compute the pseudo velocity at delts
		for (ii=0;ii<natom;ii++)
		{
		    vx[ii] += (deltsby2*fxs[ii]*1.0e-5/aw[ii]);
		    vy[ii] += (deltsby2*fys[ii]*1.0e-5/aw[ii]);
		    vz[ii] += (deltsby2*fzs[ii]*1.0e-5/aw[ii]);
		}
    }

    // inter forces, long ranged
    erfrc();

    // use the new forces to calculate the new velocities at t+delt
    ukin = 0.0;
    for (ii=0;ii<natom;ii++)
    {
		vx[ii] += (deltby2*fxl[ii]*1.0e-5/aw[ii]);
		vy[ii] += (deltby2*fyl[ii]*1.0e-5/aw[ii]);
		vz[ii] += (deltby2*fzl[ii]*1.0e-5/aw[ii]);
		ukin += aw[ii]*(vx[ii]*vx[ii]+vy[ii]*vy[ii]+vz[ii]*vz[ii]);
    }
    ukin = ukin/2.0;

    // calculate instant temperature
    tinst = 2.0*ukin*rRgas/nfree;
    
    return 0;
}

int do_accumu()
{
    // accumulators
    accumulator[0][0] += utot;
    accumulator[0][1] += utot*utot;
    accumulator[1][0] += upot;
    accumulator[1][1] += upot*upot;
    accumulator[2][0] += ukin;
    accumulator[2][1] += ukin*upot;
    accumulator[3][0] += uinter;
    accumulator[3][1] += uinter*uinter;
    accumulator[4][0] += uintra;
    accumulator[4][1] += uintra*uintra;
    accumulator[5][0] += uvdw;
    accumulator[5][1] += uvdw*uvdw;
    accumulator[6][0] += ubond;
    accumulator[6][1] += ubond*ubond;
    accumulator[7][0] += uangle;
    accumulator[7][1] += uangle*uangle;
    accumulator[8][0] += udih;
    accumulator[8][1] += udih*udih;
    accumulator[9][0] += uimp;
    accumulator[9][1] += uimp*uimp;
    accumulator[10][0] += uewald;
    accumulator[10][1] += uewald*uewald;
    accumulator[11][0] += ureal;
    accumulator[11][1] += ureal*ureal;
    accumulator[12][0] += ufourier;
    accumulator[12][1] += ufourier*ufourier;
    accumulator[13][0] += uself;
    accumulator[13][1] += uself*uself;
    accumulator[14][0] += usflj;
    accumulator[14][1] += usflj*usflj;
    accumulator[15][0] += tinst;
    accumulator[15][1] += tinst*tinst;
    accumulator[16][0] += uvacuum;
    accumulator[16][1] += uvacuum*uvacuum;
    accumulator[17][0] += uwolf;
    accumulator[17][1] += uwolf*uwolf;
    accumulator[18][0] += pinst;
    accumulator[18][1] += pinst*pinst;
    accumulator[19][0] += boxv;
    accumulator[19][1] += boxv*boxv;
    accumulator[20][0] += pideal;
    accumulator[20][1] += pideal*pideal;

    return 0;
}

int md ()
{
    // start of the MD simulation
    // initialize velocities
    velinit(); 

    // if not new run, load from old file
    if (fStart_option!=new_run) loadit(); 

    // calculate total energies
    erfrc(); 
    rafrc();
    // print out initial values
    echo();
    // print initial properties
    printit();
    // make snapshots & movies
    trajectory();

    if (what_ensemble == nvt_run) // velocity verlet with nose hoover for NVT MD
    {
		for (istep=nstep_start;istep<=nstep;istep++) // NOTE: start from 1 and <=
		{	
		    vver_nh_3();
		    // print out, snapshot, trajectory, save
		    if (istep%nstep_print == 0) printit();
		    if (nstep_ss && istep%nstep_ss == 0) snapshot();
		    if (nstep_trj && istep%nstep_trj==0) trajectory();
		    if (istep%nstep_save==0) saveit();
	
		    icounter[11]--;
		    // if still in equilibrium run
		    // do not do averages
		    if (icounter[11]>=0) continue;
		    // accumulators
		    do_accumu();
		    if ((istep-nstep_eq)%nstep_ave==0) averages();
		}
    }
    else if (what_ensemble == npt_run) // baronstat for NPT MD
    {
		for (istep=nstep_start;istep<=nstep;istep++) // NOTE: start from 1 and <=
		{	
		    npt_respa();
		    // print out, snapshot, trajectory, save
		    if (istep%nstep_print == 0) printit();
		    if (nstep_ss && istep%nstep_ss == 0) snapshot();
		    if (nstep_trj && istep%nstep_trj==0) trajectory();
		    if (istep%nstep_save==0) saveit();
	
		    icounter[11]--;
		    // if still in equilibrium run
		    // do not do averages
		    if (icounter[11]>=0) continue;
		    // accumulators
		    do_accumu();
		    if ((istep-nstep_eq)%nstep_ave==0) averages();
		}
    }
    else // velocity verlet for normal NVE MD
    {
		for (istep=nstep_start;istep<=nstep;istep++) // NOTE: start from 1 and <=
		{	
		    vver(); 
		    // print out, snapshot, trajectory, save
		    if (istep%nstep_print == 0) printit();
		    if (nstep_ss && istep%nstep_ss == 0) snapshot();
		    if (nstep_trj && istep%nstep_trj==0) trajectory();
		    if (istep%nstep_save==0) saveit();
	
		    icounter[11]--;
		    // if still in equilibrium run
		    // do not do averages
		    if (icounter[11]>=0) continue;
		    // accumulators
		    do_accumu();
		    if ((istep-nstep_eq)%nstep_ave==0) averages();
		}
    }

    return 0;
}

