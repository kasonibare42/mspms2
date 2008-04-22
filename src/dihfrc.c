/**
 * Project: mspms2
 * File: dihfrc.c
 * 
 * Copyright (C) 2008    Yang Wang <ywangd@gmail.com>
 * Created @ Apr 23, 2008
 * 
 * Description:
 * 
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include "mspms2.h"

int dihfrc(int iSpecie, int iDih, int iabs, double *uij, double *virial_ij)
{
	PSAMPLE_MOLECULE pSampleMole;
	int mm1, mm2, mm3, mm4;
	int iia, iib, iic, iid;
	double xdab, ydab, zdab, xdbc, ydbc, zdbc, xdcd, ydcd, zdcd;
	double xab, yab, zab, xbc, ybc, zbc, xcd, ycd, zcd, rrbc;
	double xac, yac, zac, pbx, pby, pbz, pb2, rpb1, rpb2;
	double pcx, pcy, pcz, pc2, rpc1, rpc2, pbpc, cost, sint;
	double theta, rsint, v1, v2, v3, v4, gamma;
	double fax, fay, faz, fcx, fcy, fcz;
	double fb1x, fb1y, fb1z, fd1x, fd1y, fd1z;
	double vopls, dopls;

	int ii1, ii2, ii3, ii4;
	double kphi, nperiod, delta0;
	double rxjk, ryjk, rzjk, rxij, ryij, rzij;
	double rxjl, ryjl, rzjl, rxik, ryik, rzik;
	double rxkl, rykl, rzkl;
	double aa1, bb1, cc1;
	double aa2, bb2, cc2;
	double temp, temp1, temp2, temp3;
	double phi, alpha_temp;
	double fij;
	
	pSampleMole = sample_mole + iSpecie;

	mm1 = pSampleMole->dih_idx[iDih][0]; // relative atom id 1
	mm2 = pSampleMole->dih_idx[iDih][1]; // Relative aotm id 2
	mm3 = pSampleMole->dih_idx[iDih][2]; // Relative atom id 3
	mm4 = pSampleMole->dih_idx[iDih][3]; // relative atom id 4
	// The atom list of dihedrals must be ordered as 0-1-2-3.
	// Otherwise, the results won't be correct.
	iia = iabs + mm1; // absolute atom id 1
	iib = iabs + mm2; // absolute atom id 2
	iic = iabs + mm3; // absolute atom id 3
	iid = iabs + mm4; // absolute atom id 4
	switch (pSampleMole->dih_type[iDih])
	// switch for different dihedral type
	{
	case DIH_NONE:
		break;
	case DIH_OPLS_COSIN:
		xdab = xx[iia] - xx[iib];
		ydab = yy[iia] - yy[iib];
		zdab = zz[iia] - zz[iib];

		xdbc = xx[iib] - xx[iic];
		ydbc = yy[iib] - yy[iic];
		zdbc = zz[iib] - zz[iic];

		xdcd = xx[iic] - xx[iid];
		ydcd = yy[iic] - yy[iid];
		zdcd = zz[iic] - zz[iid];

		xab = xdab;
		yab = ydab;
		zab = zdab;

		xbc = xdbc;
		ybc = ydbc;
		zbc = zdbc;
		rrbc = 1.0/sqrt(xbc*xbc+ybc*ybc+zbc*zbc);

		xcd = xdcd;
		ycd = ydcd;
		zcd = zdcd;

		xac = xab+xbc;
		yac = yab+ybc;
		zac = zab+zbc;

		// construct first dihedral vector
		pbx = yab*zbc - zab*ybc;
		pby = zab*xbc - xab*zbc;
		pbz = xab*ybc - yab*xbc;
		pb2 = pbx*pbx+pby*pby+pbz*pbz;
		rpb1 = 1.0/sqrt(pb2);
		rpb2 = rpb1*rpb1;

		// construct second dihedral vector
		pcx = ybc*zcd - zbc*ycd;
		pcy = zbc*xcd - xbc*zcd;
		pcz = xbc*ycd - ybc*xcd;
		pc2 = pcx*pcx+pcy*pcy+pcz*pcz;
		rpc1 = 1.0/sqrt(pc2);
		rpc2 = rpc1*rpc1;

		// determine diehdral angle
		pbpc = pbx*pcx+pby*pcy+pbz*pcz;
		cost = pbpc*rpb1*rpc1;
		sint = (xbc*(pcy*pbz-pcz*pby)+ybc*(pbx*pcz-pbz*pcx)+zbc*(pcx*pby
				-pcy*pbx)) *(rpb1*rpc1*rrbc);

		theta = atan2(sint, cost);

		// avoid singularity in sint
		sint = copysign(fmax(1.0e-8, fabs(sint)), sint);
		rsint = 1.0/sint;

		// phase angle? theta = theta - prmdih(kk,7) ?
		v1 = pSampleMole->c1[iDih];
		v2 = pSampleMole->c2[iDih];
		v3 = pSampleMole->c3[iDih];
		v4 = pSampleMole->c4[iDih];
		// Calculate the energy
		vopls = v1+0.5*(v2*(1.0+cos(theta))+v3*(1.0-cos(2.0*theta))+v4*(1.0
				+cos(3.0*theta)));
		*uij = vopls;
		// NOTE: Set dihedral virial to zero
		*virial_ij = 0.0;
		// forces
		dopls = -0.5*(v2*sin(theta)-2.0*v3*sin(2.0*theta)+3.0*v4*sin(3.0
				*theta));
		gamma = dopls*rpb1*rpc1*rsint;

		fax = gamma*((-pcy*zbc+pcz*ybc)-pbpc*rpb2*(-pby*zbc+pbz*ybc));
		fay = gamma*((pcx*zbc-pcz*xbc)-pbpc*rpb2*(pbx*zbc-pbz*xbc));
		faz = gamma*((-pcx*ybc+pcy*xbc)-pbpc*rpb2*(-pbx*ybc+pby*xbc));

		fcx = gamma*((-pcy*zab+pcz*yab)-pbpc*rpb2*(-pby*zab+pbz*yab));
		fcy = gamma*((pcx*zab-pcz*xab)-pbpc*rpb2*(pbx*zab-pbz*xab));
		fcz = gamma*((-pcx*yab+pcy*xab)-pbpc*rpb2*(-pbx*yab+pby*xab));

		fb1x= gamma*((-pby*zcd+pbz*ycd)-pbpc*rpc2*(-pcy*zcd+pcz*ycd));
		fb1y= gamma*((pbx*zcd-pbz*xcd)-pbpc*rpc2*(pcx*zcd-pcz*xcd));
		fb1z= gamma*((-pbx*ycd+pby*xcd)-pbpc*rpc2*(-pcx*ycd+pcy*xcd));

		fd1x= gamma*((-pby*zbc+pbz*ybc)-pbpc*rpc2*(-pcy*zbc+pcz*ybc));
		fd1y= gamma*((pbx*zbc-pbz*xbc)-pbpc*rpc2*(pcx*zbc-pcz*xbc));
		fd1z= gamma*((-pbx*ybc+pby*xbc)-pbpc*rpc2*(-pcx*ybc+pcy*xbc));

		fxs[iia] += fax;
		fys[iia] += fay;
		fzs[iia] += faz;

		fxs[iib] += (-fax-fcx+fb1x);
		fys[iib] += (-fay-fcy+fb1y);
		fzs[iib] += (-faz-fcz+fb1z);

		fxs[iic] += fcx-fb1x-fd1x;
		fys[iic] += fcy-fb1y-fd1y;
		fzs[iic] += fcz-fb1z-fd1z;

		fxs[iid] += fd1x;
		fys[iid] += fd1y;
		fzs[iid] += fd1z;
		break;
	case DIH_CHARMM: // future implementation 
		// modified from Shi Wei's code
		// NOT tested.
		// calculate the dihedral angle energy 
		// Note that add count_dih_multiple due to multiple matching parameters of the 
		// dihedral angle part
		// a-b-c-d  is i-j-k-l is 1-2-3-4
		// c1 is kphi, c2 is nperiod, c3 is delta0;
		kphi = pSampleMole->c1[iDih];
		nperiod = pSampleMole->c2[iDih];
		delta0 = pSampleMole->c3[iDih];
		// ** the following calculate the angle i-j-k-l 
		// ** between i-j-k plane and j-k-l plane
		if (fabs(kphi) >= TOLERANCE)
		{
			ii1 = iia;
			ii2 = iib;
			ii3 = iic;
			ii4 = iid;
			rxjk = xx[ii3] - xx[ii2];
			ryjk = yy[ii3] - yy[ii2];
			rzjk = zz[ii3] - zz[ii2];
			rxij = xx[ii1] - xx[ii2];
			ryij = yy[ii1] - yy[ii2];
			rzij = zz[ii1] - zz[ii2];
			// ** calculate a1,b1,c1 for the vector in x, y , and z direction.
			// ** This vector is perpendicular to the plane of i-j-K and obtained from
			// ** the cross product of jk and ij, ie. jk*ji
			// ** rotate from jk vector to ji vector
			aa1 = ryjk*rzij - rzjk*ryij;
			bb1 = rzjk*rxij - rxjk*rzij;
			cc1 = rxjk*ryij - rxij*ryjk;
			rxjl = xx[ii4] - xx[ii2];
			ryjl = yy[ii4] - yy[ii2];
			rzjl = zz[ii4] - zz[ii2];
			// ** calculate a2,b2, and c2 in x, y, z direction
			// ** this vector is the cross product of jl*jk 
			// ** rotate from jl to jk
			aa2 = ryjl*rzjk - rzjl*ryjk;
			bb2 = rzjl*rxjk - rxjl*rzjk;
			cc2 = rxjl*ryjk - rxjk*ryjl;
			temp1 = sqrt(aa1*aa1+bb1*bb1+cc1*cc1);
			temp2 = sqrt(aa2*aa2+bb2*bb2+cc2*cc2);
			temp3 = (aa1*aa2+bb1*bb2+cc1*cc2);
			phi = temp3/(temp1*temp2);
			alpha_temp = phi;
			// ** the angle between the two planes is below
			phi = 180.0 - (acos(phi)*180.0/pi);
			*uij = kphi*(1+cos((nperiod*phi-delta0)/180*pi));
			// NOTE: Set dihedral virial to zero
			*virial_ij = 0.0;
			// ** calculate forces due to dihedral angle part
			// ** the following make sure that it is not divided by 0, if it is, it can
			// ** be reduced from 0/0 by mathematical tricks. Note that the 
			// ** delta0(ii) is eother 0 or 180 from the parameter data file.  
			if (fabs(sin(phi*pi/180)) <= 1e-10)
			{
				temp=-nperiod*nperiod*kphi*cos((nperiod*phi-delta0)*pi/180)
						/alpha_temp;
			}
			else
			{
				temp = nperiod*kphi*sin((nperiod*phi-delta0)*pi/180)
						/sin(phi*pi/180);
			}
			// ** calculate forces on atom 1
			fij = temp/(temp2*temp1*temp1);
			fxs[ii1] = fxs[ii1] + fij*((bb2*rzjk-cc2*ryjk)*temp1-temp3
					/temp1*(bb1*rzjk-cc1*ryjk));
			fys[ii1] = fys[ii1] + fij*((-aa2*rzjk+cc2*rxjk)*temp1-temp3
					/temp1*(-aa1*rzjk+cc1*rxjk));
			fzs[ii1] = fzs[ii1] + fij*((aa2*ryjk-bb2*rxjk)*temp1-temp3
					/temp1*(aa1*ryjk-bb1*rxjk));
			// ** calculate forces on atom 2
			fij = temp/(temp1*temp1*temp2*temp2);
			rxik = xx[ii1] - xx[ii3];
			ryik = yy[ii1] - yy[ii3];
			rzik = zz[ii1] - zz[ii3];
			rxkl = xx[ii4] - xx[ii3];
			rykl = yy[ii4] - yy[ii3];
			rzkl = zz[ii4] - zz[ii3];
			fxs[ii2] = fxs[ii2] + fij
					*((bb2*rzik-bb1*rzkl-cc2*ryik+cc1*rykl)*temp1*temp2
							-temp3*(temp2/temp1*(bb1*rzik-cc1*ryik)+temp1
									/temp2*(-bb2*rzkl+cc2*rykl)));
			fys[ii2] = fys[ii2] + fij*((-aa2*rzik+aa1*rzkl+cc2*rxik-cc1
					*rxkl)*temp1*temp2 -temp3*(temp2/temp1*(-aa1*rzik+cc1
					*rxik)+temp1/temp2*(aa2*rzkl-cc2*rxkl)));
			fzs[ii2] = fzs[ii2] + fij
					*((aa2*ryik-aa1*rykl-bb2*rxik+bb1*rxkl)*temp1*temp2
							-temp3*(temp2/temp1*(aa1*ryik-bb1*rxik)+temp1
									/temp2*(-aa2*rykl+bb2*rxkl)));
			// ** calculate the forces on atom 3
			fxs[ii3] = fxs[ii3] + fij*((-bb2*rzij+bb1*rzjl+cc2*ryij-cc1
					*ryjl)*temp1*temp2 -temp3*(temp2/temp1*(-bb1*rzij+cc1
					*ryij)+temp1/temp2*(bb2*rzjl-cc2*ryjl)));
			fys[ii3] = fys[ii3] + fij
					*((aa2*rzij-aa1*rzjl-cc2*rxij+cc1*rxjl)*temp1*temp2
							-temp3*(temp2/temp1*(aa1*rzij-cc1*rxij)+temp1
									/temp2*(-aa2*rzjl+cc2*rxjl)));
			fzs[ii3]=fzs[ii3] + fij*((-aa2*ryij+aa1*ryjl+bb2*rxij-bb1*rxjl)
					*temp1*temp2 -temp3*(temp2/temp1*(-aa1*ryij+bb1*rxij)
					+temp1/temp2*(aa2*ryjl-bb2*rxjl)));
			// ** calculate the forces on atom 4
			fij = temp/(temp1*temp2*temp2);
			fxs[ii4] = fxs[ii4] + fij*((-bb1*rzjk+cc1*ryjk)*temp2-temp3
					/temp2*(-bb2*rzjk+cc2*ryjk));
			fys[ii4] = fys[ii4] + fij*((aa1*rzjk-cc1*rxjk)*temp2-temp3
					/temp2*(aa2*rzjk-cc2*rxjk));
			fzs[ii4] = fzs[ii4] + fij*((-aa1*ryjk+bb1*rxjk)*temp2-temp3
					/temp2*(-aa2*ryjk+bb2*rxjk));
		}
		break;
	default:
		printf("unknown dihedral type.\n");
		exit(1);
	} // switch for different dihedral type

	return 0;
}

