// functions for NPT-RESPA
// modified from Dr. Maginn's code
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include "mspms2.h"

int npt_nh_operator()
{
	int ii, iSpecie, iAtom;
	double mvsq, AA, BB, pdiff;
	double N_plus1RT = (nfree+1)*treq;
	double one_3N = 1.0+3.0/nfree;
	double r_mass;
	// G - force
	// Q - Mass
	// v - velocity
	// r - position

	mvsq = 0.0;
	for (ii=0; ii<natom; ii++)
	{
		// Get specie ID and relative atom ID
		get_specie_and_relative_atom_id(ii, &iSpecie, &iAtom);
		r_mass = 1.0/sample_mole[iSpecie].aw[iAtom];
		mvsq += (vx[ii]*vx[ii]+vy[ii]*vy[ii]+vz[ii]*vz[ii])/r_mass;
	}

	// the virial is 3.0*real_virial
	// preq here should be just the external pressure?

	// calculate the difference 3*V*(P_internal - P_external)
	pdiff = virial_inter + virial_intra + pljlrc*boxv*3.0 - preq*boxv*3.0; 

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

	for (ii=0; ii<natom; ii++)
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
	vts = vts + dt_outer4*Gts;

	// extra enery from the thermo/barostat for conserve energy

	utsbs = 0.5*Qbs*vbs*vbs + 0.5*Qts*vts*vts + (nfree+1)*treq*rts + preq*boxv;

	// printf("3   utsbs=%lf Gts=%lf Gbs=%lf vts=%lf BB=%lf vbs=%lf\n",utsbs,Gts,Gbs,vts,BB,vbs);

	// calcualte kinetic energy and temperature
	ukin = mvsq/2.0;
	tinst = 2.0*ukin/nfree;

	// calculate the total pressure
	pideal = natom/boxv*tinst;
	// pljlrc is calculated already after the box size changes
	// and does not change during above calculations
	pinst = pideal + (virial_inter+virial_intra)/3.0/boxv + pljlrc;

	// printf("tinst=%lf  vts=%lf  vbs=%lf\n",tinst,vts,vbs);

	return 0;
}

int npt_respa()
{
	int ii, ll, iSpecie, iAtom;
	double AA;
	double expfactor;
	double r_mass;

	npt_nh_operator();

	for (ii=0; ii<natom; ii++)
	{
		// Get specie ID and relative atom ID
		get_specie_and_relative_atom_id(ii, &iSpecie, &iAtom);
		r_mass = 1.0/sample_mole[iSpecie].aw[iAtom];
		vx[ii] += (deltby2*fxl[ii]*r_mass);
		vy[ii] += (deltby2*fyl[ii]*r_mass);
		vz[ii] += (deltby2*fzl[ii]*r_mass);
	}

	for (ll=0; ll<nstep_inner; ll++)
	{
		// AA = exp(deltsby2*vbs*1.0e-15);  // deltsby2 = dt_inner2
		AA = exp(deltsby2*vbs); // deltsby2 = dt_inner2

		for (ii=0; ii<natom; ii++)
		{
			// Get specie ID and relative atom ID
			get_specie_and_relative_atom_id(ii, &iSpecie, &iAtom);
			r_mass = 1.0/sample_mole[iSpecie].aw[iAtom];
			vx[ii] += (deltsby2*fxs[ii]*r_mass);
			vy[ii] += (deltsby2*fys[ii]*r_mass);
			vz[ii] += (deltsby2*fzs[ii]*r_mass);

			xx[ii] = xx[ii]*AA;
			yy[ii] = yy[ii]*AA;
			zz[ii] = zz[ii]*AA;

			xx[ii] = xx[ii] + delts*vx[ii];
			yy[ii] = yy[ii] + delts*vy[ii];
			zz[ii] = zz[ii] + delts*vz[ii];

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
 
		// Re-calculate the long range corrections due to box size changes
		calculate_ljlrc();

		// Re-calculate box size related variables for ewald summation
		if (iChargeType == ELECTROSTATIC_EWALD)
		{
			Vfactor_ewald = 2.0*pi/boxv;
			TWOPI_LX = 2.0*pi/boxlx;
			TWOPI_LY = 2.0*pi/boxly;
			TWOPI_LZ = 2.0*pi/boxlz;
			// 1D ewald constant
			twopi_over_3v = 2.0*pi/3.0/boxv;
		}
		
		// intra forces, short ranged
		frcshort();

		// compute the pseudo velocity at delts
		for (ii=0; ii<natom; ii++)
		{
			// Get specie ID and relative atom ID
			get_specie_and_relative_atom_id(ii, &iSpecie, &iAtom);
			r_mass = 1.0/sample_mole[iSpecie].aw[iAtom];
			vx[ii] += (deltsby2*fxs[ii]*r_mass);
			vy[ii] += (deltsby2*fys[ii]*r_mass);
			vz[ii] += (deltsby2*fzs[ii]*r_mass);
		}
	}

	// inter forces, long ranged
	frclong();

	// use the new forces to calculate the new velocities at t+delt
	for (ii=0; ii<natom; ii++)
	{
		// Get specie ID and relative atom ID
		get_specie_and_relative_atom_id(ii, &iSpecie, &iAtom);
		r_mass = 1.0/sample_mole[iSpecie].aw[iAtom];
		vx[ii] += (deltby2*fxl[ii]*r_mass);
		vy[ii] += (deltby2*fyl[ii]*r_mass);
		vz[ii] += (deltby2*fzl[ii]*r_mass);
	}

	npt_nh_operator();

	return 0;
}
