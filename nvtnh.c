#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include "vars.h"

extern int erfrc();
extern int rafrc();

// init the variables needed for NVT simulation
// Also read parameters from the input file
int init_nvt()
{
    int ii;
    const int datalen = 200;
    char buffer[200];
    char keyword[100];

    fprintf(stderr,"Reading input data for MD NVT simulation...\n");
    fprintf(fpouts,"Reading input data for MD NVT simulation...\n");

    // re-open input file to read extra data section
    fpins = fopen(INPUT,"r");

    while (fgets(buffer,datalen,fpins)!=NULL)
    {
	sscanf(buffer,"%s",keyword);
	for (ii=0;ii<strlen(keyword);ii++)
	    keyword[ii] = toupper(keyword[ii]);
	if (!strcmp(keyword,"MDNVT"))
	{
	    fprintf(stderr,"Data section for MD NVT simulation found...\n");
	    fprintf(fpouts,"Data section for MD NVT simulation found...\n");
	    sscanf(fgets(buffer,datalen,fpins), "%lf %lf", &qq, &qqs);

	    // nose hoover
	    // following for Dr. Maginn's nose hoover
	    // Gts = 0.0;
	    // vts = 0.0;
	    // rts = 0.0;
	    /*
	     * Qts = Rgas*treq*nfree/Omega
	     * where Omega is a parameter related to the mass of the thermostat
	     * for this program, we read in the Qts directly.
	     */

	    // variables for nose-hoover NVT, see frenkel and smit
	    delt_sqby2 = delt*delt/2.0;
	    delts_sqby2 = delts*delts/2.0;
	    unhts = 0.0;
	    gg = nfree; // need double check
	    ss = 0.0;
	    ps = 0.0;
	    ggs = nfree;
	    sss = 0.0;
	    pss = 0.0;
	    unhtss = 0.0;

	    fclose(fpins);
	    return 0;
	} // if keyword found
    } // read through lines
    fprintf(stderr,"Error: data for MD NVT not found.\n");
    fprintf(fpouts,"Error: data for MD NVT not found.\n");
    fclose(fpins);
    exit(1);
}


// velocity verlet with nose hoover thermostat
// based on Frenkel and Smit's codes
// use one thermostat for long range forces
// ignore the inner step
// calculate ss, ps after inner step is done
int vver_nh_2()
{
    int ii, ll;

    double delps;
    double sumv2;
    double err = 1.0e-10;
    double psn, pso;
    int	ready;
    int iter;
    double ri;
    double di; 
    int ipart;

    for (ii=0;ii<natom;ii++)
    {
	// the factor of 1.0e-5 is based on Angstrom (from the force)
	// and femto second (from delt)
	vx[ii] = vx[ii] + deltby2*(fxl[ii]/aw[ii]*1.0e-5 - ps*vx[ii]);
	vy[ii] = vy[ii] + deltby2*(fyl[ii]/aw[ii]*1.0e-5 - ps*vy[ii]);
	vz[ii] = vz[ii] + deltby2*(fzl[ii]/aw[ii]*1.0e-5 - ps*vz[ii]);
    }

    for (ll=0;ll<nstep_inner;ll++)
    {
	for (ii=0;ii<natom;ii++)
	{
	    vx[ii] = vx[ii] + deltsby2*fxs[ii]/aw[ii]*1.0e-5;
	    vy[ii] = vy[ii] + deltsby2*fys[ii]/aw[ii]*1.0e-5;
	    vz[ii] = vz[ii] + deltsby2*fzs[ii]/aw[ii]*1.0e-5;
	    xx[ii] = xx[ii] + delts*vx[ii]*1.0e-5;
	    yy[ii] = yy[ii] + delts*vy[ii]*1.0e-5;
	    zz[ii] = zz[ii] + delts*vz[ii]*1.0e-5;
	}

	// intra forces, short ranged
	rafrc();

	sumv2 = 0.0;
	for (ii=0;ii<natom;ii++)
	{
	    vx[ii] = vx[ii] + deltsby2*fxs[ii]*1.0e-5/aw[ii];
	    vy[ii] = vy[ii] + deltsby2*fys[ii]*1.0e-5/aw[ii];
	    vz[ii] = vz[ii] + deltsby2*fzs[ii]*1.0e-5/aw[ii];
	    sumv2 = sumv2 + aw[ii]*(vx[ii]*vx[ii]+vy[ii]*vy[ii]+vz[ii]*vz[ii]);
	}
    }
    // ps has unit of  1/(femto second) =  1/fs
    // qq has unit of kg*m^2/mol*1.0*e-30
    // ss has no unit
    ss = ss + ps*delt + (sumv2-gg*treq*Rgas)*delt_sqby2/qq;
    ps = ps + (sumv2-gg*treq*Rgas)*deltby2/qq;

    // inter forces, long ranged
    erfrc();

    sumv2 = 0.0;
    for (ii=0;ii<natom;ii++)
    {
	vxn[ii] = vx[ii];
	vyn[ii] = vy[ii];
	vzn[ii] = vz[ii];
	sumv2 = sumv2 + aw[ii]*(vxn[ii]*vxn[ii]+vyn[ii]*vyn[ii]+vzn[ii]*vzn[ii]);
    }
    psn = ps;
    ready = false;
    iter = 0;
    while (ready==false && iter<100)
    {
	iter++;
	pso = psn;
	delps = 0.0;
	for (ii=0;ii<natom;ii++)
	{
	    vxo[ii] = vxn[ii];
	    vyo[ii] = vyn[ii];
	    vzo[ii] = vzn[ii];

	    bx[ii] = -deltby2*(fxl[ii]/aw[ii]*1.0e-5 - pso*vxo[ii]) - (vx[ii]-vxo[ii]);
	    ri = aw[ii]*vxo[ii]*delt/qq;
	    delps = delps + ri*bx[ii];

	    by[ii] = -deltby2*(fyl[ii]/aw[ii]*1.0e-5 - pso*vyo[ii]) - (vy[ii]-vyo[ii]);
	    ri = aw[ii]*vyo[ii]*delt/qq;
	    delps = delps + ri*by[ii];

	    bz[ii] = -deltby2*(fzl[ii]/aw[ii]*1.0e-5 - pso*vzo[ii]) - (vz[ii]-vzo[ii]);
	    ri = aw[ii]*vzo[ii]*delt/qq;
	    delps = delps + ri*bz[ii];
	}
	di = -(pso*deltby2 + 1.0);
	delps = delps - di*((-sumv2+gg*treq*Rgas)*deltby2/qq - (ps-pso));
	delps = delps/(-delt*deltby2*sumv2/qq + di);

	sumv2 = 0.0;
	for (ii=0;ii<natom;ii++)
	{
	    vxn[ii] = vxn[ii] + (bx[ii] + deltby2*vxo[ii]*delps)/di;
	    vyn[ii] = vyn[ii] + (by[ii] + deltby2*vyo[ii]*delps)/di;
	    vzn[ii] = vzn[ii] + (bz[ii] + deltby2*vzo[ii]*delps)/di;
	    sumv2 = sumv2 + aw[ii]*(vxn[ii]*vxn[ii]+vyn[ii]*vyn[ii]+vzn[ii]*vzn[ii]);
	}
	psn = pso + delps;
	// test for convergence
	ready = true;
	ipart = 0;
	while (ipart<=natom && ready==true) // NOTE: <=
	{
	    ipart++;
	    if (ipart<=natom)
	    {
		if (fabs((vxn[ii]-vxo[ii])/vxn[ii]) > err)
		    ready = false;
		if (fabs((vyn[ii]-vyo[ii])/vyn[ii]) > err)
		    ready = false;
		if (fabs((vzn[ii]-vzo[ii])/vzn[ii]) > err)
		    ready = false;
	    }
	    else if (fabs((psn-pso)/psn) > err)
		ready = false;
	}
    } // end of first while

    for (ii=0;ii<natom;ii++)
    {
	vx[ii] = vxn[ii];
	vy[ii] = vyn[ii];
	vz[ii] = vzn[ii];
    }
    ps = psn;
    ukin = sumv2/2.0;
    // H = ukin + upot + (ps*ps*qq)/2 + gg*treq*Rgas*ss;

    // energy of thermostat
    unhts = (ps*ps*qq)/2.0 + gg*treq*Rgas*ss;

    // calculate instant temperature
    tinst = 2.0*ukin*rRgas/nfree;
    
    return 0;
}






// velocity verlet with nose hoover thermostat
// based on Frenkel and Smit's codes
// use one thermostat for both outter and inner steps
int vver_nh_1()
{
    int ii, ll;

    double delps;
    double sumv2;
    double err = 1.0e-10;
    double psn, pso;
    int	ready;
    int iter;
    double ri;
    double di; 
    int ipart;

    for (ii=0;ii<natom;ii++)
    {
	// the factor of 1.0e-5 is based on Angstrom (from the force)
	// and femto second (from delt)
	vx[ii] = vx[ii] + deltby2*(fxl[ii]/aw[ii]*1.0e-5 - ps*vx[ii]);
	vy[ii] = vy[ii] + deltby2*(fyl[ii]/aw[ii]*1.0e-5 - ps*vy[ii]);
	vz[ii] = vz[ii] + deltby2*(fzl[ii]/aw[ii]*1.0e-5 - ps*vz[ii]);
    }

    for (ll=0;ll<nstep_inner;ll++)
    {
	sumv2 = 0.0;
	for (ii=0;ii<natom;ii++)
	{
	    sumv2 = sumv2 + aw[ii]*(vx[ii]*vx[ii]+vy[ii]*vy[ii]+vz[ii]*vz[ii]);

	    vx[ii] = vx[ii] + deltsby2*fxs[ii]/aw[ii]*1.0e-5;
	    vy[ii] = vy[ii] + deltsby2*fys[ii]/aw[ii]*1.0e-5;
	    vz[ii] = vz[ii] + deltsby2*fzs[ii]/aw[ii]*1.0e-5;

	    xx[ii] = xx[ii] + delts*vx[ii]*1.0e-5;
	    yy[ii] = yy[ii] + delts*vy[ii]*1.0e-5;
	    zz[ii] = zz[ii] + delts*vz[ii]*1.0e-5;

	}
	ss = ss + ps*delts + (sumv2-gg*treq*Rgas)*delts_sqby2/qq;
	ps = ps + (sumv2-gg*treq*Rgas)*deltsby2/qq;

	// intra forces, short ranged
	rafrc();

	sumv2 = 0.0;
	for (ii=0;ii<natom;ii++)
	{
	    vxn[ii] = vx[ii];
	    vyn[ii] = vy[ii];
	    vzn[ii] = vz[ii];
	    sumv2 = sumv2 + aw[ii]*(vxn[ii]*vxn[ii]+vyn[ii]*vyn[ii]+vzn[ii]*vzn[ii]);
	}
	psn = ps;
	ready = false;
	iter = 0;

	while (ready==false && iter<100)
	{
	    iter++;
	    pso = psn;
	    delps = 0.0;
	    for (ii=0;ii<natom;ii++)
	    {
		vxo[ii] = vxn[ii];
		vyo[ii] = vyn[ii];
		vzo[ii] = vzn[ii];

		bx[ii] = -deltsby2*(fxs[ii]/aw[ii]*1.0e-5 - pso*vxo[ii]) - (vx[ii]-vxo[ii]);
		ri = aw[ii]*vxo[ii]*delts/qq;
		delps = delps + ri*bx[ii];

		by[ii] = -deltsby2*(fys[ii]/aw[ii]*1.0e-5 - pso*vyo[ii]) - (vy[ii]-vyo[ii]);
		ri = aw[ii]*vyo[ii]*delts/qq;
		delps = delps + ri*by[ii];

		bz[ii] = -deltsby2*(fzs[ii]/aw[ii]*1.0e-5 - pso*vzo[ii]) - (vz[ii]-vzo[ii]);
		ri = aw[ii]*vzo[ii]*delts/qq;
		delps = delps + ri*bz[ii];
	    }
	    di = -(pso*deltsby2 + 1.0);
	    delps = delps - di*((-sumv2+gg*treq*Rgas)*deltsby2/qq - (ps-pso));
	    delps = delps/(-delts*deltsby2*sumv2/qq + di);

	    sumv2 = 0.0;
	    for (ii=0;ii<natom;ii++)
	    {
		vxn[ii] = vxn[ii] + (bx[ii] + deltsby2*vxo[ii]*delps)/di;
		vyn[ii] = vyn[ii] + (by[ii] + deltsby2*vyo[ii]*delps)/di;
		vzn[ii] = vzn[ii] + (bz[ii] + deltsby2*vzo[ii]*delps)/di;
		sumv2 = sumv2 + aw[ii]*(vxn[ii]*vxn[ii]+vyn[ii]*vyn[ii]+vzn[ii]*vzn[ii]);
	    }
	    psn = pso + delps;
	    // test for convergence
	    ready = true;
	    ipart = 0;
	    while (ipart<=natom && ready==true) // NOTE: <=
	    {
		ipart++;
		if (ipart<=natom)
		{
		    if (fabs((vxn[ii]-vxo[ii])/vxn[ii]) > err)
			ready = false;
		    if (fabs((vyn[ii]-vyo[ii])/vyn[ii]) > err)
			ready = false;
		    if (fabs((vzn[ii]-vzo[ii])/vzn[ii]) > err)
			ready = false;
		}
		else if (fabs((psn-pso)/psn) > err)
		    ready = false;
	    }
	} // end of while

	for (ii=0;ii<natom;ii++)
	{
	    vx[ii] = vxn[ii];
	    vy[ii] = vyn[ii];
	    vz[ii] = vzn[ii];
	}
	ps = psn;
    }
    // ps has unit of  1/(femto second) =  1/fs
    // qq has unit of kg*m^2/mol*1.0*e-30
    // ss has no unit

    // inter forces, long ranged
    erfrc();

    sumv2 = 0.0;
    for (ii=0;ii<natom;ii++)
    {
	vxn[ii] = vx[ii];
	vyn[ii] = vy[ii];
	vzn[ii] = vz[ii];
	sumv2 = sumv2 + aw[ii]*(vxn[ii]*vxn[ii]+vyn[ii]*vyn[ii]+vzn[ii]*vzn[ii]);
    }
    psn = ps;
    ready = false;
    iter = 0;
    while (ready==false && iter<100)
    {
	iter++;
	pso = psn;
	delps = 0.0;
	for (ii=0;ii<natom;ii++)
	{
	    vxo[ii] = vxn[ii];
	    vyo[ii] = vyn[ii];
	    vzo[ii] = vzn[ii];

	    bx[ii] = -deltby2*(fxl[ii]/aw[ii]*1.0e-5 - pso*vxo[ii]) - (vx[ii]-vxo[ii]);
	    ri = aw[ii]*vxo[ii]*delt/qq;
	    delps = delps + ri*bx[ii];

	    by[ii] = -deltby2*(fyl[ii]/aw[ii]*1.0e-5 - pso*vyo[ii]) - (vy[ii]-vyo[ii]);
	    ri = aw[ii]*vyo[ii]*delt/qq;
	    delps = delps + ri*by[ii];

	    bz[ii] = -deltby2*(fzl[ii]/aw[ii]*1.0e-5 - pso*vzo[ii]) - (vz[ii]-vzo[ii]);
	    ri = aw[ii]*vzo[ii]*delt/qq;
	    delps = delps + ri*bz[ii];
	}
	di = -(pso*deltby2 + 1.0);
	delps = delps - di*((-sumv2+gg*treq*Rgas)*deltby2/qq - (ps-pso));
	delps = delps/(-delt*deltby2*sumv2/qq + di);

	sumv2 = 0.0;
	for (ii=0;ii<natom;ii++)
	{
	    vxn[ii] = vxn[ii] + (bx[ii] + deltby2*vxo[ii]*delps)/di;
	    vyn[ii] = vyn[ii] + (by[ii] + deltby2*vyo[ii]*delps)/di;
	    vzn[ii] = vzn[ii] + (bz[ii] + deltby2*vzo[ii]*delps)/di;
	    sumv2 = sumv2 + aw[ii]*(vxn[ii]*vxn[ii]+vyn[ii]*vyn[ii]+vzn[ii]*vzn[ii]);
	}
	psn = pso + delps;
	// test for convergence
	ready = true;
	ipart = 0;
	while (ipart<=natom && ready==true) // NOTE: <=
	{
	    ipart++;
	    if (ipart<=natom)
	    {
		if (fabs((vxn[ii]-vxo[ii])/vxn[ii]) > err)
		    ready = false;
		if (fabs((vyn[ii]-vyo[ii])/vyn[ii]) > err)
		    ready = false;
		if (fabs((vzn[ii]-vzo[ii])/vzn[ii]) > err)
		    ready = false;
	    }
	    else if (fabs((psn-pso)/psn) > err)
		ready = false;
	}
    } // end of first while

    for (ii=0;ii<natom;ii++)
    {
	vx[ii] = vxn[ii];
	vy[ii] = vyn[ii];
	vz[ii] = vzn[ii];
    }
    ps = psn;
    ukin = sumv2/2.0;
    // H = ukin + upot + (ps*ps*qq)/2 + gg*treq*Rgas*ss;

    // energy of thermostat
    unhts = (ps*ps*qq)/2.0 + gg*treq*Rgas*ss;

    // calculate instant temperature
    tinst = 2.0*ukin*rRgas/nfree;
    
    return 0;
}








// velocity verlet with nose hoover thermostat
// based on Frenkel and Smit's codes
// use two thermostats for outter and inner steps
int vver_nh_3()
{
    int ii, ll;

    double delps;
    double delpss;
    double sumv2;
    double err = 1.0e-10;
    double psn, pso;
    double pssn, psso;
    int	ready;
    int iter;
    double ri;
    double di; 
    int ipart;

    for (ii=0;ii<natom;ii++)
    {
	// the factor of 1.0e-5 is based on Angstrom (from the force)
	// and femto second (from delt)
	vx[ii] = vx[ii] + deltby2*(fxl[ii]/aw[ii]*1.0e-5 - ps*vx[ii]);
	vy[ii] = vy[ii] + deltby2*(fyl[ii]/aw[ii]*1.0e-5 - ps*vy[ii]);
	vz[ii] = vz[ii] + deltby2*(fzl[ii]/aw[ii]*1.0e-5 - ps*vz[ii]);
    }

    for (ll=0;ll<nstep_inner;ll++)
    {
	sumv2 = 0.0;
	for (ii=0;ii<natom;ii++)
	{
	    sumv2 = sumv2 + aw[ii]*(vx[ii]*vx[ii]+vy[ii]*vy[ii]+vz[ii]*vz[ii]);

	    vx[ii] = vx[ii] + deltsby2*(fxs[ii]/aw[ii]*1.0e-5 - pss*vx[ii]);
	    vy[ii] = vy[ii] + deltsby2*(fys[ii]/aw[ii]*1.0e-5 - pss*vy[ii]);
	    vz[ii] = vz[ii] + deltsby2*(fzs[ii]/aw[ii]*1.0e-5 - pss*vz[ii]);

	    xx[ii] = xx[ii] + delts*vx[ii]*1.0e-5;
	    yy[ii] = yy[ii] + delts*vy[ii]*1.0e-5;
	    zz[ii] = zz[ii] + delts*vz[ii]*1.0e-5;
	}
	sss = sss + pss*delts + (sumv2-ggs*treq*Rgas)*delts_sqby2/qqs;
	pss = pss + (sumv2-ggs*treq*Rgas)*deltsby2/qqs;

	// intra forces, short ranged
	rafrc();

	sumv2 = 0.0;
	for (ii=0;ii<natom;ii++)
	{
	    vxn[ii] = vx[ii];
	    vyn[ii] = vy[ii];
	    vzn[ii] = vz[ii];
	    sumv2 = sumv2 + aw[ii]*(vxn[ii]*vxn[ii]+vyn[ii]*vyn[ii]+vzn[ii]*vzn[ii]);
	}
	pssn = pss;
	ready = false;
	iter = 0;

	while (ready==false && iter<100)
	{
	    iter++;
	    psso = pssn;
	    delpss = 0.0;
	    for (ii=0;ii<natom;ii++)
	    {
		vxo[ii] = vxn[ii];
		vyo[ii] = vyn[ii];
		vzo[ii] = vzn[ii];

		bx[ii] = -deltsby2*(fxs[ii]/aw[ii]*1.0e-5 - psso*vxo[ii]) - (vx[ii]-vxo[ii]);
		ri = aw[ii]*vxo[ii]*delts/qqs;
		delpss = delpss + ri*bx[ii];

		by[ii] = -deltsby2*(fys[ii]/aw[ii]*1.0e-5 - psso*vyo[ii]) - (vy[ii]-vyo[ii]);
		ri = aw[ii]*vyo[ii]*delts/qqs;
		delpss = delpss + ri*by[ii];

		bz[ii] = -deltsby2*(fzs[ii]/aw[ii]*1.0e-5 - psso*vzo[ii]) - (vz[ii]-vzo[ii]);
		ri = aw[ii]*vzo[ii]*delts/qqs;
		delpss = delpss + ri*bz[ii];
	    }
	    di = -(psso*deltsby2 + 1.0);
	    delpss = delpss - di*((-sumv2+ggs*treq*Rgas)*deltsby2/qqs - (pss-psso));
	    delpss = delpss/(-delts*deltsby2*sumv2/qqs + di);

	    sumv2 = 0.0;
	    for (ii=0;ii<natom;ii++)
	    {
		vxn[ii] = vxn[ii] + (bx[ii] + deltsby2*vxo[ii]*delpss)/di;
		vyn[ii] = vyn[ii] + (by[ii] + deltsby2*vyo[ii]*delpss)/di;
		vzn[ii] = vzn[ii] + (bz[ii] + deltsby2*vzo[ii]*delpss)/di;
		sumv2 = sumv2 + aw[ii]*(vxn[ii]*vxn[ii]+vyn[ii]*vyn[ii]+vzn[ii]*vzn[ii]);
	    }
	    pssn = psso + delpss;
	    // test for convergence
	    ready = true;
	    ipart = 0;
	    while (ipart<=natom && ready==true) // NOTE: <=
	    {
		ipart++;
		if (ipart<=natom)
		{
		    if (fabs((vxn[ii]-vxo[ii])/vxn[ii]) > err)
			ready = false;
		    if (fabs((vyn[ii]-vyo[ii])/vyn[ii]) > err)
			ready = false;
		    if (fabs((vzn[ii]-vzo[ii])/vzn[ii]) > err)
			ready = false;
		}
		else if (fabs((pssn-psso)/pssn) > err)
		    ready = false;
	    }
	} // end of while

	for (ii=0;ii<natom;ii++)
	{
	    vx[ii] = vxn[ii];
	    vy[ii] = vyn[ii];
	    vz[ii] = vzn[ii];
	}
	pss = pssn;

	// energy of inner thermostat
	unhtss = (pss*pss*qqs)/2.0 + ggs*treq*Rgas*sss;
    }
    // ps has unit of  1/(femto second) =  1/fs
    // qq has unit of kg*m^2/mol*1.0*e-30
    // ss has no unit
    ss = ss + ps*delt + (sumv2-gg*treq*Rgas)*delt_sqby2/qq;
    ps = ps + (sumv2-gg*treq*Rgas)*deltby2/qq;

    // inter forces, long ranged
    erfrc();

    sumv2 = 0.0;
    for (ii=0;ii<natom;ii++)
    {
	vxn[ii] = vx[ii];
	vyn[ii] = vy[ii];
	vzn[ii] = vz[ii];
	sumv2 = sumv2 + aw[ii]*(vxn[ii]*vxn[ii]+vyn[ii]*vyn[ii]+vzn[ii]*vzn[ii]);
    }
    psn = ps;
    ready = false;
    iter = 0;
    while (ready==false && iter<100)
    {
	iter++;
	pso = psn;
	delps = 0.0;
	for (ii=0;ii<natom;ii++)
	{
	    vxo[ii] = vxn[ii];
	    vyo[ii] = vyn[ii];
	    vzo[ii] = vzn[ii];

	    bx[ii] = -deltby2*(fxl[ii]/aw[ii]*1.0e-5 - pso*vxo[ii]) - (vx[ii]-vxo[ii]);
	    ri = aw[ii]*vxo[ii]*delt/qq;
	    delps = delps + ri*bx[ii];

	    by[ii] = -deltby2*(fyl[ii]/aw[ii]*1.0e-5 - pso*vyo[ii]) - (vy[ii]-vyo[ii]);
	    ri = aw[ii]*vyo[ii]*delt/qq;
	    delps = delps + ri*by[ii];

	    bz[ii] = -deltby2*(fzl[ii]/aw[ii]*1.0e-5 - pso*vzo[ii]) - (vz[ii]-vzo[ii]);
	    ri = aw[ii]*vzo[ii]*delt/qq;
	    delps = delps + ri*bz[ii];
	}
	di = -(pso*deltby2 + 1.0);
	delps = delps - di*((-sumv2+gg*treq*Rgas)*deltby2/qq - (ps-pso));
	delps = delps/(-delt*deltby2*sumv2/qq + di);

	sumv2 = 0.0;
	for (ii=0;ii<natom;ii++)
	{
	    vxn[ii] = vxn[ii] + (bx[ii] + deltby2*vxo[ii]*delps)/di;
	    vyn[ii] = vyn[ii] + (by[ii] + deltby2*vyo[ii]*delps)/di;
	    vzn[ii] = vzn[ii] + (bz[ii] + deltby2*vzo[ii]*delps)/di;
	    sumv2 = sumv2 + aw[ii]*(vxn[ii]*vxn[ii]+vyn[ii]*vyn[ii]+vzn[ii]*vzn[ii]);
	}
	psn = pso + delps;
	// test for convergence
	ready = true;
	ipart = 0;
	while (ipart<=natom && ready==true) // NOTE: <=
	{
	    ipart++;
	    if (ipart<=natom)
	    {
		if (fabs((vxn[ii]-vxo[ii])/vxn[ii]) > err)
		    ready = false;
		if (fabs((vyn[ii]-vyo[ii])/vyn[ii]) > err)
		    ready = false;
		if (fabs((vzn[ii]-vzo[ii])/vzn[ii]) > err)
		    ready = false;
	    }
	    else if (fabs((psn-pso)/psn) > err)
		ready = false;
	}
    } // end of first while

    for (ii=0;ii<natom;ii++)
    {
	vx[ii] = vxn[ii];
	vy[ii] = vyn[ii];
	vz[ii] = vzn[ii];
    }
    ps = psn;
    ukin = sumv2/2.0;
    // H = ukin + upot + (ps*ps*qq)/2 + gg*treq*Rgas*ss;

    // energy of thermostat
    unhts = (ps*ps*qq)/2.0 + gg*treq*Rgas*ss;

    // calculate instant temperature
    tinst = 2.0*ukin*rRgas/nfree;
    
    return 0;
}

// Following codes are modified from Dr. Maginn's sample codes.
int nvt_nh_operator()
{
    int ii;
    double mvsq;
    double AA;

    mvsq = 0.0;
    for (ii=0;ii<natom;ii++)
	mvsq += aw[ii]*(vx[ii]*vx[ii]+vy[ii]*vy[ii]+vz[ii]*vz[ii]);

    // compute the driving force for the thermostat
    Gts = (mvsq - nfree*Rgas*treq)/Qts;

    // advance the thermostat velocity 1/4 time step
    vts = vts + dt_outer4*Gts;

    // advance the thermostat position 1/2 timestep
    rts = rts + dt_outer2*vts;

    // advance velocities of atoms @ 1/2 time step
    AA = exp(-dt_outer2*vts);

    for (ii=0;ii<natom;ii++)
    {
	vx[ii] *= AA;
	vy[ii] *= AA;
	vz[ii] *= AA;
    }

    // compute new driving force for the thermostat and advance velocities @ 1/4 time step
    mvsq = mvsq*AA*AA;

    Gts = (mvsq - nfree*Rgas*treq)/Qts;

    vts = vts + dt_outer4*Gts;

    // extra energy for the thermostat for conserving energy
    unhts = 0.5*Gts*vts*vts + treq*Rgas*rts*nfree;

    // calculate kinetic energy and temperature
    ukin = 0.5*mvsq;
    tinst = 2.0*ukin*rRgas/nfree;
    
    return 0;
}

int nvt_respa() // velocity verlet
{
    int ii, ll;

    nvt_nh_operator();

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
    for (ii=0;ii<natom;ii++)
    {
	vx[ii] += (deltby2*fxl[ii]*1.0e-5/aw[ii]);
	vy[ii] += (deltby2*fyl[ii]*1.0e-5/aw[ii]);
	vz[ii] += (deltby2*fzl[ii]*1.0e-5/aw[ii]);
    }

    nvt_nh_operator();
    
    return 0;
}
