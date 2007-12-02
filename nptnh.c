// functions for NPT-RESPA
// modified from Dr. Maginn's code
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include "vars.h"

extern int erfrc();
extern int rafrc();
extern int calculate_ljlrc();

int init_npt_respa()
{
    int ii;
    const int datalen = 200;
    char buffer[200];
    char keyword[100];

    fprintf(stderr,"Reading input data for MD NPT simulation...\n");
    fprintf(fpouts,"Reading input data for MD NPT simulation...\n");

    // re-open input file to read extra data section
    fpins = fopen(INPUT,"r");

    while (fgets(buffer,datalen,fpins)!=NULL)
    {
	sscanf(buffer,"%s",keyword);
	for (ii=0;ii<strlen(keyword);ii++)
	    keyword[ii] = toupper(keyword[ii]);
	if (!strcmp(keyword,"MDNPT"))
	{
	    fprintf(stderr,"Data section for MD NPT simulation found...\n");
	    fprintf(fpouts,"Data section for MD NPT simulation found...\n");
	    sscanf(fgets(buffer,datalen,fpins), "%lf", &preq);
	    sscanf(fgets(buffer,datalen,fpins), "%lf %lf", &Qts, &Qbs);

	    // thermo/barostat
	    utsbs = 0.0;
	    vts = sqrt((nfree+1.0)*Rgas/Qts);
	    vbs = sqrt((nfree+1.0)*Rgas/Qbs);
	    rts = 0.0;

	    // extra enery from the thermo/barostat for conserve energy
	    // utsbs = 0.5*Qbs*vbs*vbs + 0.5*Qts*vts*vts + (nfree+1)*Rgas*treq*rts*1.0e-10 + preq*boxv*PascalA3_to_J_mol;

	    utsbs = 0.5*Qbs*vbs*vbs + 0.5*Qts*vts*vts + (nfree+1)*Rgas*treq*rts + preq*boxv*PascalA3_to_J_mol;
 
	    // printf("initiate   vts=%lf  vbs=%lf  utsbs=%lf\n",vts,vbs,utsbs);

	    fclose(fpins);
	    return 0;
	} // if keyword found
    } // read through lines
    fprintf(stderr,"Error: data for MD NPT not found.\n");
    fprintf(fpouts,"Error: data for MD NPT not found.\n");
    fclose(fpins);
    exit(1);
}

int npt_nh_operator()
{
    int ii;
    double mvsq, AA, BB, pdiff;
    double N_plus1RT = (nfree+1)*Rgas*treq;
    double one_3N = 1.0+3.0/nfree;
    // G - force
    // Q - Mass
    // v - velocity
    // r - position

    mvsq = 0.0;
    for (ii=0;ii<natom;ii++)
	mvsq += aw[ii]*(vx[ii]*vx[ii]+vy[ii]*vy[ii]+vz[ii]*vz[ii]);

    // the virial is 3.0*real_virial, unit is J/mol
    // preq here should be just the external pressure?

    // calculate the long range corrections
    calculate_ljlrc();

    // calculate the difference 3*V*(P_internal - P_external)
    pdiff = virial_inter + virial_intra + pljlrc*boxv*3.0*PascalA3_to_J_mol // turn pascal to J/mol
	- preq*boxv*3.0*PascalA3_to_J_mol; //6.0221415e-7 is Na*1e-30 turn preq*boxv to J/mol

    // printf("0   pdiff=%lf Gts=%lf vts=%lf vbs=%lf\n",pdiff,Gts,vts,vbs);

    // thermostat driving force
    // mvsq is (kg/mol)*(m/s)^2 = J/mol
    // Qts is kg/mol
    // Gts should have the unit of m/s^2
    // Qts is (kg/mol)
    // So the equation below should have a hidden 1/m, i.e. divided by unit meter
    // the unit is (m^2/s^2)/m
    Gts = (mvsq + Qbs*vbs*vbs - N_plus1RT)/Qts; 

    // thermostat velocity
    // vts should have unit of m/s
    // dt_outer4 is fs = 1.0e-15
    // the 2nd term should be multiplied by 1.0e-5 to get the unit right
    // vts = vts + dt_outer4*Gts*1.0e-15;
    vts = vts + dt_outer4*Gts;

    // barostat velocity
    // vts is m/s, dt_outer8 is fs
    // BB = exp(-vts*dt_outer8*1.0e-15);
    BB = exp(-vts*dt_outer8);
    vbs = vbs*BB;

    // printf("1   pdiff=%lf Gts=%lf vts=%lf BB=%lf vbs=%lf\n",pdiff,Gts,vts,BB,vbs);

    // barostat driving force
    // Gbs is m/s^2 and Gbs should be (kg/mol)*m
    Gbs = (one_3N*mvsq + pdiff)/Qbs;
    // vbs = vbs + dt_outer4*Gbs*1.0e-15;
    vbs = vbs + dt_outer4*Gbs;

    vbs = vbs*BB;

    // printf("2   pdiff=%lf Gts=%lf Gbs=%lf vts=%lf BB=%lf vbs=%lf\n",pdiff,Gts,Gbs,vts,BB,vbs);

    // thermostat position
    // rts should have unit of angstrom
    // deltby2 is fs, vts is m/s
    // rts = rts + deltby2*vts*1.0e-5; // deltby2 is dt_outer2
    rts = rts + deltby2*vts; // deltby2 is dt_outer2


    // dt_outer is fs, vts, vbs is m/s
    // AA = exp(-dt_outer2*(vts+one_3N*vbs)*1.0e-15);
    AA = exp(-dt_outer2*(vts+one_3N*vbs));

    // printf("rts=%lf  AA=%lf\n",rts,AA);

    for (ii=0;ii<natom;ii++)
    {
	vx[ii] = vx[ii]*AA;
	vy[ii] = vy[ii]*AA;
	vz[ii] = vz[ii]*AA;
    }

    vbs = vbs*BB;

    // new mv^2
    mvsq = mvsq*AA*AA;

    // recompute barostat force
    Gbs = (one_3N*mvsq + pdiff)/Qbs;

    // new barostat velocity
    // vbs = vbs + dt_outer4*Gbs*1.0e-15;
    vbs = vbs + dt_outer4*Gbs;

    vbs = vbs*BB;

    // new thermostat force
    Gts = (mvsq + Qbs*vbs*vbs - (N_plus1RT) ) / Qts;

    // thermostat velocity
    // vts  = vts + dt_outer4*Gts*1.0e-15;
    vts  = vts + dt_outer4*Gts;

    // extra enery from the thermo/barostat for conserve energy
    // utsbs = 0.5*Qbs*vbs*vbs + 0.5*Qts*vts*vts + (nfree+1)*Rgas*treq*rts*1.0e-10 + preq*boxv*PascalA3_to_J_mol;

    utsbs = 0.5*Qbs*vbs*vbs + 0.5*Qts*vts*vts + (nfree+1)*Rgas*treq*rts + preq*boxv*PascalA3_to_J_mol;

    // printf("3   utsbs=%lf Gts=%lf Gbs=%lf vts=%lf BB=%lf vbs=%lf\n",utsbs,Gts,Gbs,vts,BB,vbs);

    // calcualte kinetic energy and temperature
    ukin = mvsq/2.0;
    tinst = 2.0*ukin*rRgas/nfree;

    // calculate the total pressure
    pideal=natom/boxv*tinst*kb_1e30;
    // pljlrc is calculated already at the beginning of this function
    // and did not change during above calculations
    pinst = pideal
	+ (virial_inter+virial_intra)*virial_to_pressure/boxv
	+ pljlrc;

    // printf("tinst=%lf  vts=%lf  vbs=%lf\n",tinst,vts,vbs);

    return 0;
}

int npt_respa() 
{
    int ii, ll;
    double AA;
    double expfactor;

    npt_nh_operator();

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
	// AA = exp(deltsby2*vbs*1.0e-15);  // deltsby2 = dt_inner2
	AA = exp(deltsby2*vbs);  // deltsby2 = dt_inner2

	for (ii=0;ii<natom;ii++)
	{
	    vx[ii] += (deltsby2*fxs[ii]*1.0e-5/aw[ii]);
	    vy[ii] += (deltsby2*fys[ii]*1.0e-5/aw[ii]);
	    vz[ii] += (deltsby2*fzs[ii]*1.0e-5/aw[ii]);

	    xx[ii] = xx[ii]*AA;
	    yy[ii] = yy[ii]*AA;
	    zz[ii] = zz[ii]*AA;

	    xx[ii] = xx[ii] + delts*vx[ii]*1.0e-5;
	    yy[ii] = yy[ii] + delts*vy[ii]*1.0e-5;
	    zz[ii] = zz[ii] + delts*vz[ii]*1.0e-5;

	    xx[ii] = xx[ii]*AA;
	    yy[ii] = yy[ii]*AA;
	    zz[ii] = zz[ii]*AA;
	}

	// adjust box volume
	// expfactor = exp(delts*vbs*1.0e-15);
	expfactor = exp(delts*vbs);
	boxlx = boxlx*expfactor;
	boxly = boxly*expfactor;
	boxlz = boxlz*expfactor;
	boxv = boxlx*boxly*boxlz;

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
    // ukin = 0.0;
    for (ii=0;ii<natom;ii++)
    {
	vx[ii] += (deltby2*fxl[ii]*1.0e-5/aw[ii]);
	vy[ii] += (deltby2*fyl[ii]*1.0e-5/aw[ii]);
	vz[ii] += (deltby2*fzl[ii]*1.0e-5/aw[ii]);
	// ukin += aw[ii]*(vx[ii]*vx[ii]+vy[ii]*vy[ii]+vz[ii]*vz[ii]);
    }
    // ukin = ukin/2.0;

    // calculate instant temperature
    // tinst = 2.0*ukin*rRgas/nfree;

    npt_nh_operator();

    return 0;
}
