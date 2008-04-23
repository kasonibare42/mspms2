#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include "mspms2.h"


void rezero_nvt_ts()
{
	// rezero nvt related variables
	unhts = 0.0;
	ss = 0.0;
	ps = 0.0;
	sss = 0.0;
	pss = 0.0;
	unhtss = 0.0;
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
	int ready;
	int iter;
	double ri;
	double di;
	int ipart;

	for (ii=0; ii<natom; ii++)
	{
		// the factor of 1.0e-5 is based on Angstrom (from the force)
		// and femto second (from delt)
		vx[ii] = vx[ii] + deltby2*(fxl[ii]/aw[ii]*1.0e-5 - ps*vx[ii]);
		vy[ii] = vy[ii] + deltby2*(fyl[ii]/aw[ii]*1.0e-5 - ps*vy[ii]);
		vz[ii] = vz[ii] + deltby2*(fzl[ii]/aw[ii]*1.0e-5 - ps*vz[ii]);
	}

	for (ll=0; ll<nstep_inner; ll++)
	{
		for (ii=0; ii<natom; ii++)
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
		for (ii=0; ii<natom; ii++)
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
	ss = ss + ps*delt + (sumv2-gg*treq*RGAS)*delt_sqby2/qq;
	ps = ps + (sumv2-gg*treq*RGAS)*deltby2/qq;

	// inter forces, long ranged
	erfrc();

	sumv2 = 0.0;
	for (ii=0; ii<natom; ii++)
	{
		vxn[ii] = vx[ii];
		vyn[ii] = vy[ii];
		vzn[ii] = vz[ii];
		sumv2 = sumv2 + aw[ii]
				*(vxn[ii]*vxn[ii]+vyn[ii]*vyn[ii]+vzn[ii]*vzn[ii]);
	}
	psn = ps;
	ready = false;
	iter = 0;
	while (ready==false && iter<100)
	{
		iter++;
		pso = psn;
		delps = 0.0;
		for (ii=0; ii<natom; ii++)
		{
			vxo[ii] = vxn[ii];
			vyo[ii] = vyn[ii];
			vzo[ii] = vzn[ii];

			bx[ii] = -deltby2*(fxl[ii]/aw[ii]*1.0e-5 - pso*vxo[ii]) - (vx[ii]
					-vxo[ii]);
			ri = aw[ii]*vxo[ii]*delt/qq;
			delps = delps + ri*bx[ii];

			by[ii] = -deltby2*(fyl[ii]/aw[ii]*1.0e-5 - pso*vyo[ii]) - (vy[ii]
					-vyo[ii]);
			ri = aw[ii]*vyo[ii]*delt/qq;
			delps = delps + ri*by[ii];

			bz[ii] = -deltby2*(fzl[ii]/aw[ii]*1.0e-5 - pso*vzo[ii]) - (vz[ii]
					-vzo[ii]);
			ri = aw[ii]*vzo[ii]*delt/qq;
			delps = delps + ri*bz[ii];
		}
		di = -(pso*deltby2 + 1.0);
		delps = delps - di*((-sumv2+gg*treq*RGAS)*deltby2/qq - (ps-pso));
		delps = delps/(-delt*deltby2*sumv2/qq + di);

		sumv2 = 0.0;
		for (ii=0; ii<natom; ii++)
		{
			vxn[ii] = vxn[ii] + (bx[ii] + deltby2*vxo[ii]*delps)/di;
			vyn[ii] = vyn[ii] + (by[ii] + deltby2*vyo[ii]*delps)/di;
			vzn[ii] = vzn[ii] + (bz[ii] + deltby2*vzo[ii]*delps)/di;
			sumv2 = sumv2 + aw[ii]*(vxn[ii]*vxn[ii]+vyn[ii]*vyn[ii]+vzn[ii]
					*vzn[ii]);
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

	for (ii=0; ii<natom; ii++)
	{
		vx[ii] = vxn[ii];
		vy[ii] = vyn[ii];
		vz[ii] = vzn[ii];
	}
	ps = psn;
	ukin = sumv2/2.0;
	// H = ukin + upot + (ps*ps*qq)/2 + gg*treq*RGAS*ss;

	// energy of thermostat
	unhts = (ps*ps*qq)/2.0 + gg*treq*RGAS*ss;

	// calculate instant temperature
	tinst = 2.0*ukin*R_RGAS/nfree;

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
	int ready;
	int iter;
	double ri;
	double di;
	int ipart;

	for (ii=0; ii<natom; ii++)
	{
		// the factor of 1.0e-5 is based on Angstrom (from the force)
		// and femto second (from delt)
		vx[ii] = vx[ii] + deltby2*(fxl[ii]/aw[ii]*1.0e-5 - ps*vx[ii]);
		vy[ii] = vy[ii] + deltby2*(fyl[ii]/aw[ii]*1.0e-5 - ps*vy[ii]);
		vz[ii] = vz[ii] + deltby2*(fzl[ii]/aw[ii]*1.0e-5 - ps*vz[ii]);
	}

	for (ll=0; ll<nstep_inner; ll++)
	{
		sumv2 = 0.0;
		for (ii=0; ii<natom; ii++)
		{
			sumv2 = sumv2 + aw[ii]*(vx[ii]*vx[ii]+vy[ii]*vy[ii]+vz[ii]*vz[ii]);

			vx[ii] = vx[ii] + deltsby2*fxs[ii]/aw[ii]*1.0e-5;
			vy[ii] = vy[ii] + deltsby2*fys[ii]/aw[ii]*1.0e-5;
			vz[ii] = vz[ii] + deltsby2*fzs[ii]/aw[ii]*1.0e-5;

			xx[ii] = xx[ii] + delts*vx[ii]*1.0e-5;
			yy[ii] = yy[ii] + delts*vy[ii]*1.0e-5;
			zz[ii] = zz[ii] + delts*vz[ii]*1.0e-5;

		}
		ss = ss + ps*delts + (sumv2-gg*treq*RGAS)*delts_sqby2/qq;
		ps = ps + (sumv2-gg*treq*RGAS)*deltsby2/qq;

		// intra forces, short ranged
		rafrc();

		sumv2 = 0.0;
		for (ii=0; ii<natom; ii++)
		{
			vxn[ii] = vx[ii];
			vyn[ii] = vy[ii];
			vzn[ii] = vz[ii];
			sumv2 = sumv2 + aw[ii]*(vxn[ii]*vxn[ii]+vyn[ii]*vyn[ii]+vzn[ii]
					*vzn[ii]);
		}
		psn = ps;
		ready = false;
		iter = 0;

		while (ready==false && iter<100)
		{
			iter++;
			pso = psn;
			delps = 0.0;
			for (ii=0; ii<natom; ii++)
			{
				vxo[ii] = vxn[ii];
				vyo[ii] = vyn[ii];
				vzo[ii] = vzn[ii];

				bx[ii] = -deltsby2*(fxs[ii]/aw[ii]*1.0e-5 - pso*vxo[ii])
						- (vx[ii]-vxo[ii]);
				ri = aw[ii]*vxo[ii]*delts/qq;
				delps = delps + ri*bx[ii];

				by[ii] = -deltsby2*(fys[ii]/aw[ii]*1.0e-5 - pso*vyo[ii])
						- (vy[ii]-vyo[ii]);
				ri = aw[ii]*vyo[ii]*delts/qq;
				delps = delps + ri*by[ii];

				bz[ii] = -deltsby2*(fzs[ii]/aw[ii]*1.0e-5 - pso*vzo[ii])
						- (vz[ii]-vzo[ii]);
				ri = aw[ii]*vzo[ii]*delts/qq;
				delps = delps + ri*bz[ii];
			}
			di = -(pso*deltsby2 + 1.0);
			delps = delps - di*((-sumv2+gg*treq*RGAS)*deltsby2/qq - (ps-pso));
			delps = delps/(-delts*deltsby2*sumv2/qq + di);

			sumv2 = 0.0;
			for (ii=0; ii<natom; ii++)
			{
				vxn[ii] = vxn[ii] + (bx[ii] + deltsby2*vxo[ii]*delps)/di;
				vyn[ii] = vyn[ii] + (by[ii] + deltsby2*vyo[ii]*delps)/di;
				vzn[ii] = vzn[ii] + (bz[ii] + deltsby2*vzo[ii]*delps)/di;
				sumv2 = sumv2 + aw[ii]*(vxn[ii]*vxn[ii]+vyn[ii]*vyn[ii]+vzn[ii]
						*vzn[ii]);
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

		for (ii=0; ii<natom; ii++)
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
	for (ii=0; ii<natom; ii++)
	{
		vxn[ii] = vx[ii];
		vyn[ii] = vy[ii];
		vzn[ii] = vz[ii];
		sumv2 = sumv2 + aw[ii]
				*(vxn[ii]*vxn[ii]+vyn[ii]*vyn[ii]+vzn[ii]*vzn[ii]);
	}
	psn = ps;
	ready = false;
	iter = 0;
	while (ready==false && iter<100)
	{
		iter++;
		pso = psn;
		delps = 0.0;
		for (ii=0; ii<natom; ii++)
		{
			vxo[ii] = vxn[ii];
			vyo[ii] = vyn[ii];
			vzo[ii] = vzn[ii];

			bx[ii] = -deltby2*(fxl[ii]/aw[ii]*1.0e-5 - pso*vxo[ii]) - (vx[ii]
					-vxo[ii]);
			ri = aw[ii]*vxo[ii]*delt/qq;
			delps = delps + ri*bx[ii];

			by[ii] = -deltby2*(fyl[ii]/aw[ii]*1.0e-5 - pso*vyo[ii]) - (vy[ii]
					-vyo[ii]);
			ri = aw[ii]*vyo[ii]*delt/qq;
			delps = delps + ri*by[ii];

			bz[ii] = -deltby2*(fzl[ii]/aw[ii]*1.0e-5 - pso*vzo[ii]) - (vz[ii]
					-vzo[ii]);
			ri = aw[ii]*vzo[ii]*delt/qq;
			delps = delps + ri*bz[ii];
		}
		di = -(pso*deltby2 + 1.0);
		delps = delps - di*((-sumv2+gg*treq*RGAS)*deltby2/qq - (ps-pso));
		delps = delps/(-delt*deltby2*sumv2/qq + di);

		sumv2 = 0.0;
		for (ii=0; ii<natom; ii++)
		{
			vxn[ii] = vxn[ii] + (bx[ii] + deltby2*vxo[ii]*delps)/di;
			vyn[ii] = vyn[ii] + (by[ii] + deltby2*vyo[ii]*delps)/di;
			vzn[ii] = vzn[ii] + (bz[ii] + deltby2*vzo[ii]*delps)/di;
			sumv2 = sumv2 + aw[ii]*(vxn[ii]*vxn[ii]+vyn[ii]*vyn[ii]+vzn[ii]
					*vzn[ii]);
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

	for (ii=0; ii<natom; ii++)
	{
		vx[ii] = vxn[ii];
		vy[ii] = vyn[ii];
		vz[ii] = vzn[ii];
	}
	ps = psn;
	ukin = sumv2/2.0;
	// H = ukin + upot + (ps*ps*qq)/2 + gg*treq*RGAS*ss;

	// energy of thermostat
	unhts = (ps*ps*qq)/2.0 + gg*treq*RGAS*ss;

	// calculate instant temperature
	tinst = 2.0*ukin*R_RGAS/nfree;

	return 0;
}

// velocity verlet with nose hoover thermostat
// based on Frenkel and Smit's codes
// use two thermostats for outter and inner steps
int vver_nh_3()
{
	int ii, ll;
	int iSpecie, iAtom;
	double r_mass;
	double delps;
	double delpss;
	double sumv2;
	double err = 1.0e-10;
	double psn, pso;
	double pssn, psso;
	int ready;
	int iter;
	double ri;
	double di;
	int ipart;

	for (ii=0; ii<natom; ii++)
	{
		// Get specie ID and relative atom ID
		get_specie_and_relative_atom_id(ii, &iSpecie, &iAtom);
		r_mass = 1.0/sample_mole[iSpecie].aw[iAtom];
		vx[ii] = vx[ii] + deltby2*(fxl[ii]*r_mass - ps*vx[ii]);
		vy[ii] = vy[ii] + deltby2*(fyl[ii]*r_mass - ps*vy[ii]);
		vz[ii] = vz[ii] + deltby2*(fzl[ii]*r_mass - ps*vz[ii]);
	}

	for (ll=0; ll<nstep_inner; ll++)
	{
		sumv2 = 0.0;
		for (ii=0; ii<natom; ii++)
		{
			// Get specie ID and relative atom ID
			get_specie_and_relative_atom_id(ii, &iSpecie, &iAtom);
			r_mass = 1.0/sample_mole[iSpecie].aw[iAtom];
			sumv2 = sumv2 + (vx[ii]*vx[ii]+vy[ii]*vy[ii]+vz[ii]*vz[ii])/r_mass;
			vx[ii] = vx[ii] + deltsby2*(fxs[ii]*r_mass - pss*vx[ii]);
			vy[ii] = vy[ii] + deltsby2*(fys[ii]*r_mass - pss*vy[ii]);
			vz[ii] = vz[ii] + deltsby2*(fzs[ii]*r_mass - pss*vz[ii]);
			xx[ii] = xx[ii] + delts*vx[ii];
			yy[ii] = yy[ii] + delts*vy[ii];
			zz[ii] = zz[ii] + delts*vz[ii];
		}
		sss = sss + pss*delts + (sumv2-ggs*treq)*delts_sqby2/qqs;
		pss = pss + (sumv2-ggs*treq)*deltsby2/qqs;

		// intra forces, short ranged
		rafrc();

		sumv2 = 0.0;
		for (ii=0; ii<natom; ii++)
		{
			// Get specie ID and relative atom ID
			get_specie_and_relative_atom_id(ii, &iSpecie, &iAtom);
			r_mass = 1.0/sample_mole[iSpecie].aw[iAtom];
			vxn[ii] = vx[ii];
			vyn[ii] = vy[ii];
			vzn[ii] = vz[ii];
			sumv2 = sumv2 + (vxn[ii]*vxn[ii]+vyn[ii]*vyn[ii]+vzn[ii]*vzn[ii])/r_mass;
		}
		pssn = pss;
		ready = false;
		iter = 0;
		while (ready==false && iter<100)
		{
			iter++;
			psso = pssn;
			delpss = 0.0;
			for (ii=0; ii<natom; ii++)
			{
				// Get specie ID and relative atom ID
				get_specie_and_relative_atom_id(ii, &iSpecie, &iAtom);
				r_mass = 1.0/sample_mole[iSpecie].aw[iAtom];
				vxo[ii] = vxn[ii];
				vyo[ii] = vyn[ii];
				vzo[ii] = vzn[ii];

				bx[ii] = -deltsby2*(fxs[ii]*r_mass - psso*vxo[ii])
						- (vx[ii]-vxo[ii]);
				ri = vxo[ii]*delts/qqs/r_mass;
				delpss = delpss + ri*bx[ii];

				by[ii] = -deltsby2*(fys[ii]*r_mass - psso*vyo[ii])
						- (vy[ii]-vyo[ii]);
				ri = vyo[ii]*delts/qqs/r_mass;
				delpss = delpss + ri*by[ii];

				bz[ii] = -deltsby2*(fzs[ii]*r_mass - psso*vzo[ii])
						- (vz[ii]-vzo[ii]);
				ri = vzo[ii]*delts/qqs/r_mass;
				delpss = delpss + ri*bz[ii];
			}
			di = -(psso*deltsby2 + 1.0);
			delpss = delpss - di*((-sumv2+ggs*treq)*deltsby2/qqs - (pss
					-psso));
			delpss = delpss/(-delts*deltsby2*sumv2/qqs + di);

			sumv2 = 0.0;
			for (ii=0; ii<natom; ii++)
			{
				// Get specie ID and relative atom ID
				get_specie_and_relative_atom_id(ii, &iSpecie, &iAtom);
				r_mass = 1.0/sample_mole[iSpecie].aw[iAtom];
				vxn[ii] = vxn[ii] + (bx[ii] + deltsby2*vxo[ii]*delpss)/di;
				vyn[ii] = vyn[ii] + (by[ii] + deltsby2*vyo[ii]*delpss)/di;
				vzn[ii] = vzn[ii] + (bz[ii] + deltsby2*vzo[ii]*delpss)/di;
				sumv2 = sumv2 + (vxn[ii]*vxn[ii]+vyn[ii]*vyn[ii]+vzn[ii]*vzn[ii])/r_mass;
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
					{
						ready = false;
					}
					if (fabs((vyn[ii]-vyo[ii])/vyn[ii]) > err)
					{
						ready = false;
					}
					if (fabs((vzn[ii]-vzo[ii])/vzn[ii]) > err)
					{
						ready = false;
					}
				}
				else if (fabs((pssn-psso)/pssn) > err)
				{
					ready = false;
				}
			}
		} // end of while

		for (ii=0; ii<natom; ii++)
		{
			vx[ii] = vxn[ii];
			vy[ii] = vyn[ii];
			vz[ii] = vzn[ii];
		}
		pss = pssn;

		// energy of inner thermostat
		unhtss = (pss*pss*qqs)/2.0 + ggs*treq*sss;
	}
	// ps has unit of  1/(femto second) =  1/fs
	// qq has unit of kg*m^2/mol*1.0*e-30
	// ss has no unit
	ss = ss + ps*delt + (sumv2-gg*treq)*delt_sqby2/qq;
	ps = ps + (sumv2-gg*treq)*deltby2/qq;

	// inter forces, long ranged
	erfrc();

	sumv2 = 0.0;
	for (ii=0; ii<natom; ii++)
	{
		// Get specie ID and relative atom ID
		get_specie_and_relative_atom_id(ii, &iSpecie, &iAtom);
		r_mass = 1.0/sample_mole[iSpecie].aw[iAtom];
		vxn[ii] = vx[ii];
		vyn[ii] = vy[ii];
		vzn[ii] = vz[ii];
		sumv2 = sumv2 + (vxn[ii]*vxn[ii]+vyn[ii]*vyn[ii]+vzn[ii]*vzn[ii])/r_mass;
	}
	psn = ps;
	ready = false;
	iter = 0;
	while (ready==false && iter<100)
	{
		iter++;
		pso = psn;
		delps = 0.0;
		for (ii=0; ii<natom; ii++)
		{
			// Get specie ID and relative atom ID
			get_specie_and_relative_atom_id(ii, &iSpecie, &iAtom);
			r_mass = 1.0/sample_mole[iSpecie].aw[iAtom];

			vxo[ii] = vxn[ii];
			vyo[ii] = vyn[ii];
			vzo[ii] = vzn[ii];

			bx[ii] = -deltby2*(fxl[ii]*r_mass - pso*vxo[ii]) - (vx[ii]-vxo[ii]);
			ri = vxo[ii]*delt/qq/r_mass;
			delps = delps + ri*bx[ii];

			by[ii] = -deltby2*(fyl[ii]*r_mass - pso*vyo[ii]) - (vy[ii]-vyo[ii]);
			ri = vyo[ii]*delt/qq/r_mass;
			delps = delps + ri*by[ii];

			bz[ii] = -deltby2*(fzl[ii]*r_mass - pso*vzo[ii]) - (vz[ii]-vzo[ii]);
			ri = vzo[ii]*delt/qq/r_mass;
			delps = delps + ri*bz[ii];
		}
		di = -(pso*deltby2 + 1.0);
		delps = delps - di*((-sumv2+gg*treq)*deltby2/qq - (ps-pso));
		delps = delps/(-delt*deltby2*sumv2/qq + di);

		sumv2 = 0.0;
		for (ii=0; ii<natom; ii++)
		{
			// Get specie ID and relative atom ID
			get_specie_and_relative_atom_id(ii, &iSpecie, &iAtom);
			r_mass = 1.0/sample_mole[iSpecie].aw[iAtom];

			vxn[ii] = vxn[ii] + (bx[ii] + deltby2*vxo[ii]*delps)/di;
			vyn[ii] = vyn[ii] + (by[ii] + deltby2*vyo[ii]*delps)/di;
			vzn[ii] = vzn[ii] + (bz[ii] + deltby2*vzo[ii]*delps)/di;
			sumv2 = sumv2 + (vxn[ii]*vxn[ii]+vyn[ii]*vyn[ii]+vzn[ii]*vzn[ii])/r_mass;
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
				{
					ready = false;
				}
				if (fabs((vyn[ii]-vyo[ii])/vyn[ii]) > err)
				{
					ready = false;
				}
				if (fabs((vzn[ii]-vzo[ii])/vzn[ii]) > err)
				{
					ready = false;
				}
			}
			else if (fabs((psn-pso)/psn) > err)
			{
				ready = false;
			}
		}
	} // end of first while

	for (ii=0; ii<natom; ii++)
	{
		vx[ii] = vxn[ii];
		vy[ii] = vyn[ii];
		vz[ii] = vzn[ii];
	}
	ps = psn;
	ukin = sumv2/2.0;
	// H = ukin + upot + (ps*ps*qq)/2 + gg*treq*RGAS*ss;

	// energy of thermostat
	unhts = (ps*ps*qq)/2.0 + gg*treq*ss;

	// calculate instant temperature
	tinst = 2.0*ukin/nfree;

	return 0;
}

// Following codes are modified from Dr. Maginn's sample codes.
int nvt_nh_operator()
{
	int ii;
	double mvsq;
	double AA;

	mvsq = 0.0;
	for (ii=0; ii<natom; ii++)
		mvsq += aw[ii]*(vx[ii]*vx[ii]+vy[ii]*vy[ii]+vz[ii]*vz[ii]);

	// compute the driving force for the thermostat
	Gts = (mvsq - nfree*RGAS*treq)/Qts;

	// advance the thermostat velocity 1/4 time step
	vts = vts + dt_outer4*Gts;

	// advance the thermostat position 1/2 timestep
	rts = rts + dt_outer2*vts;

	// advance velocities of atoms @ 1/2 time step
	AA = exp(-dt_outer2*vts);

	for (ii=0; ii<natom; ii++)
	{
		vx[ii] *= AA;
		vy[ii] *= AA;
		vz[ii] *= AA;
	}

	// compute new driving force for the thermostat and advance velocities @ 1/4 time step
	mvsq = mvsq*AA*AA;

	Gts = (mvsq - nfree*RGAS*treq)/Qts;

	vts = vts + dt_outer4*Gts;

	// extra energy for the thermostat for conserving energy
	unhts = 0.5*Gts*vts*vts + treq*RGAS*rts*nfree;

	// calculate kinetic energy and temperature
	ukin = 0.5*mvsq;
	tinst = 2.0*ukin*R_RGAS/nfree;

	return 0;
}

int nvt_respa() // velocity verlet
{
	int ii, ll;

	nvt_nh_operator();

	for (ii=0; ii<natom; ii++)
	{
		// the factor of 1.0e-5 is based on Angstrom (from the force)
		// and femto second (from delt)
		vx[ii] += (deltby2*fxl[ii]*1.0e-5/aw[ii]);
		vy[ii] += (deltby2*fyl[ii]*1.0e-5/aw[ii]);
		vz[ii] += (deltby2*fzl[ii]*1.0e-5/aw[ii]);
	}

	for (ll=0; ll<nstep_inner; ll++)
	{
		for (ii=0; ii<natom; ii++)
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
		for (ii=0; ii<natom; ii++)
		{
			vx[ii] += (deltsby2*fxs[ii]*1.0e-5/aw[ii]);
			vy[ii] += (deltsby2*fys[ii]*1.0e-5/aw[ii]);
			vz[ii] += (deltsby2*fzs[ii]*1.0e-5/aw[ii]);
		}
	}

	// inter forces, long ranged
	erfrc();

	// use the new forces to calculate the new velocities at t+delt
	for (ii=0; ii<natom; ii++)
	{
		vx[ii] += (deltby2*fxl[ii]*1.0e-5/aw[ii]);
		vy[ii] += (deltby2*fyl[ii]*1.0e-5/aw[ii]);
		vz[ii] += (deltby2*fzl[ii]*1.0e-5/aw[ii]);
	}

	nvt_nh_operator();

	return 0;
}
