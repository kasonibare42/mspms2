/* Calculate intra energies and forces
 include bond, angle, dihedral, impropers
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "vars.h"

int bndfrc(int iMole)
{
	int ii;
	int ii1, ii2;
	double rxij, ryij, rzij;
	double rijsq, rij;
	double delta_r;
	double ubond_temp;
	double virij;
	double fij, fxij, fyij, fzij;
	double De;
	double exp_term, one_minus_exp_term;
	double rmax2;
	int iiStart, iiEnd;
	int iPhysMoleID;
	int iSpecieID;

	// Assign the lower and upper limits for ii loop
	if (iMole == -1)
	{
		iiStart = 0;
		iiEnd = nbond_hist_max;
	}
	else if (iMole >= 0)
	{
		iiStart = mole_first_bond_idx[iMole];
		iiEnd = mole_last_bond_idx[iMole];
	}
	else
	{
		fprintf(stderr, "Error: Invalid parameters for bndfrc(%d).\n", iMole);
		fprintf(fpouts, "Error: Invalid parameters for bndfrc(%d).\n", iMole);
		exit(1);
	}

	// bond stretching calculations
	for (ii=iiStart; ii<iiEnd; ii++)
	{
		iPhysMoleID = bnd2mole[ii];
		if (mole_status[iPhysMoleID] == MOLE_STATUS_VACANCY)
		{
			continue;
		}
		
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
			gUbondSession += ubond_temp;
			fij = -Kb[ii]*delta_r/rij;
			fxij = fij*rxij;
			fyij = fij*ryij;
			fzij = fij*rzij;
			// contribution to virial
			gVirialIntraSession = gVirialIntraSession + fxij*rxij + fyij*ryij + fzij*rzij;
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
			gUbondSession += ubond_temp;
			fij = -2.0*De*alpha[ii]*one_minus_exp_term*exp_term/rij;
			fxij = fij*rxij;
			fyij = fij*ryij;
			fzij = fij*rzij;
			// contribution to virial
			gVirialIntraSession = gVirialIntraSession + fxij*rxij + fyij*ryij + fzij*rzij;
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
			gUbondSession += ubond_temp;
			fij = -2.0*Kb[ii]*rmax2/(rmax2-rijsq);
			fxij = fij*rxij;
			fyij = fij*ryij;
			fzij = fij*rzij;
			// contribution to virial
			gVirialIntraSession = gVirialIntraSession + fxij*rxij + fyij*ryij + fzij*rzij;
			// force on atom 1
			fxs[ii1] += fxij;
			fys[ii1] += fyij;
			fzs[ii1] += fzij;
			// force on atom 2
			fxs[ii2] -= fxij;
			fys[ii2] -= fyij;
			fzs[ii2] -= fzij;
			break;
		case BOND_SPRING_PATH_INTEGRAL: //4
			iSpecieID = mole2specie[iPhysMoleID];
            ubond_temp = 0.5*spring[iSpecieID]*rijsq;
            gUbondSession += ubond_temp;
			virij = -spring[iSpecieID]*rijsq;
			gVirialIntraSession += virij;
			// virbond
			fij = -spring[iSpecieID];
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

	return 0;
}

int aglfrc(int iMole)
{
	int ii;
	int iia, iib, iic;
	double xdab, ydab, zdab, xdbc, ydbc, zdbc;
	double rab, rrab, xab, yab, zab;
	double rbc, rrbc, xbc, ybc, zbc;
	double cost, theta, sint;
	double delta_theta;
	double uangle_temp;
	double gamma, fxa, fya, fza, fxc, fyc, fzc;
	double fxb, fyb, fzb;

	double vec_13_x, vec_13_y, vec_13_z;
	double vec_12_x, vec_12_y, vec_12_z;
	double vec_23_x, vec_23_y, vec_23_z;
	double r12, r13, r23;
	double delta_r1, delta_r2, delta_r3;
	double k_theta, k_r_theta, k_r_rprime;
	double u123;
	double frc_term_1;
	double frc_term_2_1, frc_term_2_2;
	double frc_term_3, frc_term_3_1, frc_term_3_2;

	int iiStart, iiEnd;
	int iPhysMoleID;

	// Assign the lower and upper limits for ii loop
	if (iMole == -1)
	{
		iiStart = 0;
		iiEnd = nangle_hist_max;
	}
	else if (iMole >= 0)
	{
		iiStart = mole_first_angle_idx[iMole];
		iiEnd = mole_last_angle_idx[iMole];
	}
	else
	{
		fprintf(stderr, "Error: Invalid parameters for aglfrc(%d).\n", iMole);
		fprintf(fpouts, "Error: Invalid parameters for aglfrc(%d).\n", iMole);
		exit(1);
	}

	// angle bending calculations
	for (ii=iiStart; ii<iiEnd; ii++)
	{
		iPhysMoleID = agl2mole[ii];
		if (mole_status[iPhysMoleID] == MOLE_STATUS_VACANCY)
		{
			continue;
		}
		
		// The atom list of the angle must be list as 0-1-2
		// otherwise, the calculations are not right
		iia = angle_idx[ii][0];
		iib = angle_idx[ii][1];
		iic = angle_idx[ii][2];
		switch (angle_type[ii])
		{
		case angle_none:
			break;
		case angle_harmonic:
			// first bond vector
			xdab = xx[iia] - xx[iib];
			ydab = yy[iia] - yy[iib];
			zdab = zz[iia] - zz[iib];
			// second bond vector
			xdbc = xx[iic] - xx[iib];
			ydbc = yy[iic] - yy[iib];
			zdbc = zz[iic] - zz[iib];
			// define components of first bond vector
			rab = sqrt(xdab*xdab + ydab*ydab + zdab*zdab);
			rrab = 1.0/rab;
			xab = xdab*rrab;
			yab = ydab*rrab;
			zab = zdab*rrab;
			// define components of second bond vector
			rbc = sqrt(xdbc*xdbc + ydbc*ydbc + zdbc*zdbc);
			rrbc = 1.0/rbc;
			xbc = xdbc*rrbc;
			ybc = ydbc*rrbc;
			zbc = zdbc*rrbc;
			// cosin angle
			cost = xab*xbc + yab*ybc+ zab*zbc;
			// angle in radian
			theta = acos(cost);
			sint = sqrt(1.0-cost*cost);
			// must use c99 standard, otherwise fmax wont work
			sint = fmax(1.0e-8, sqrt(1.0-cost*cost)); // for 180 degree
			// angle difference
			delta_theta = theta - Thetaeq[ii];
			uangle_temp = 0.5*Ktheta[ii]*delta_theta*delta_theta;
			gUangleSession += uangle_temp;
			// forces
			gamma = Ktheta[ii]*delta_theta/sint;
			fxa = gamma*(xbc-xab*cost)*rrab; // vector 1
			fya = gamma*(ybc-yab*cost)*rrab;
			fza = gamma*(zbc-zab*cost)*rrab;
			fxc = gamma*(xab-xbc*cost)*rrbc; // vector 2
			fyc = gamma*(yab-ybc*cost)*rrbc;
			fzc = gamma*(zab-zbc*cost)*rrbc;
			// for angle potentials without r as variable
			// the net contribution to virial is zero
			fxs[iia] += fxa;
			fys[iia] += fya;
			fzs[iia] += fza;
			fxs[iib] = fxs[iib] - fxa - fxc;
			fys[iib] = fys[iib] - fya - fyc;
			fzs[iib] = fzs[iib] - fza - fzc;
			fxs[iic] += fxc;
			fys[iic] += fyc;
			fzs[iic] += fzc;
			break;
		case angle_TRwater: // Toukan & Rhaman water potential
			// H-H
			vec_13_x = xx[iia] - xx[iic];
			vec_13_y = yy[iia] - yy[iic];
			vec_13_z = zz[iia] - zz[iic];
			// H-O
			vec_12_x = xx[iia] - xx[iib];
			vec_12_y = yy[iia] - yy[iib];
			vec_12_z = zz[iia] - zz[iib];
			// O-H
			vec_23_x = xx[iib] - xx[iic];
			vec_23_y = yy[iib] - yy[iic];
			vec_23_z = zz[iib] - zz[iic];
			// r_HH
			r13 = sqrt(vec_13_x*vec_13_x + vec_13_y*vec_13_y + vec_13_z
					*vec_13_z);
			// r_HO
			r12 = sqrt(vec_12_x*vec_12_x + vec_12_y*vec_12_y + vec_12_z
					*vec_12_z);
			// r_OH
			r23 = sqrt(vec_23_x*vec_23_x + vec_23_y*vec_23_y + vec_23_z
					*vec_23_z);
			delta_r3 = r13 - agl_para_4[ii]; // H-H stretch
			delta_r1 = r12 - agl_para_5[ii]; // O-H stretch
			delta_r2 = r23 - agl_para_5[ii]; // O-H stretch

			k_theta = Ktheta[ii];
			k_r_theta = Thetaeq[ii];
			k_r_rprime = agl_para_3[ii];
			u123 = 0.5*k_theta*delta_r3*delta_r3 + k_r_theta*delta_r3*(delta_r1
					+delta_r2) + k_r_rprime*delta_r1*delta_r2;
			gUangleSession += u123;
			// forces
			// atom 1 (a)
			frc_term_1 = -k_theta*delta_r3/r13;
			frc_term_2_1 = -k_r_theta*delta_r3/r12;
			frc_term_2_2 = -k_r_theta*(delta_r1+delta_r2)/r13;
			frc_term_3 = -k_r_rprime*delta_r2/r12;
			fxa = (frc_term_1*vec_13_x + frc_term_2_1*vec_12_x + frc_term_2_2
					*vec_13_x + frc_term_3*vec_12_x);
			fya = (frc_term_1*vec_13_y + frc_term_2_1*vec_12_y + frc_term_2_2
					*vec_13_y + frc_term_3*vec_12_y);
			fza = (frc_term_1*vec_13_z + frc_term_2_1*vec_12_z + frc_term_2_2
					*vec_13_z + frc_term_3*vec_12_z);
			fxs[iia] += fxa;
			fys[iia] += fya;
			fzs[iia] += fza;

			// atom 2 (b)
			// frc_term_1 = 0.0;
			frc_term_2_1 = -k_r_theta*delta_r3/r23;
			frc_term_2_2 = k_r_theta*delta_r3/r12;
			frc_term_3_1 = -k_r_rprime*delta_r1/r23;
			frc_term_3_2 = k_r_rprime*delta_r2/r12;
			fxb = (frc_term_2_1*vec_23_x + frc_term_2_2*vec_12_x + frc_term_3_1
					*vec_23_x + frc_term_3_2*vec_12_x);
			fyb = (frc_term_2_1*vec_23_y + frc_term_2_2*vec_12_y + frc_term_3_1
					*vec_23_y + frc_term_3_2*vec_12_y);
			fzb = (frc_term_2_1*vec_23_z + frc_term_2_2*vec_12_z + frc_term_3_1
					*vec_23_z + frc_term_3_2*vec_12_z);
			fxs[iib] += fxb;
			fys[iib] += fyb;
			fzs[iib] += fzb;

			// atom 3 (c)
			frc_term_1 = k_theta*delta_r3/r13;
			frc_term_2_1 = k_r_theta*delta_r3/r23;
			frc_term_2_2 = k_r_theta*(delta_r1+delta_r2)/r13;
			frc_term_3 = k_r_rprime*delta_r1/r23;
			fxc = (frc_term_1*vec_13_x + frc_term_2_1*vec_23_x + frc_term_2_2
					*vec_13_x + frc_term_3*vec_23_x);
			fyc = (frc_term_1*vec_13_y + frc_term_2_1*vec_23_y + frc_term_2_2
					*vec_13_y + frc_term_3*vec_23_y);
			fzc = (frc_term_1*vec_13_z + frc_term_2_1*vec_23_z + frc_term_2_2
					*vec_13_z + frc_term_3*vec_23_z);
			fxs[iic] += fxc;
			fys[iic] += fyc;
			fzs[iic] += fzc;

			// contribution to virial
			// Not sure. need double check
			gVirialIntraSession = gVirialIntraSession + fxa*vec_12_x + fya*vec_12_y + fza
					*vec_12_z;
			gVirialIntraSession = gVirialIntraSession - fxc*vec_23_x - fyc*vec_23_y - fzc
					*vec_23_z;

			break;
		default:
			printf("unknown angle type.\n");
			exit(1);
		} // ending of switch for different angle types
	} // end of angle bending calculations

	return 0;
}

int dihfrc(int iMole)
{
	int ii;
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

	int iiStart, iiEnd;
	int iPhysMoleID;

	// Assign the lower and upper limits for ii loop
	if (iMole == -1)
	{
		iiStart = 0;
		iiEnd = ndih_hist_max;
	}
	else if (iMole >= 0)
	{
		iiStart = mole_first_dih_idx[iMole];
		iiEnd = mole_last_dih_idx[iMole];
	}
	else
	{
		fprintf(stderr, "Error: Invalid parameters for dihfrc(%d).\n", iMole);
		fprintf(fpouts, "Error: Invalid parameters for dihfrc(%d).\n", iMole);
		exit(1);
	}

	// dihedral torsion calculations
	// dihedral has no contribution to virial
	for (ii=iiStart; ii<iiEnd; ii++)
	{
		iPhysMoleID = dih2mole[ii];
		if (mole_status[iPhysMoleID] == MOLE_STATUS_VACANCY)
		{
			continue;
		}
		
		// The atom list of dihedrals must be ordered as 0-1-2-3
		// Otherwise, the results wont be correct
		iia = dih_idx[ii][0];
		iib = dih_idx[ii][1];
		iic = dih_idx[ii][2];
		iid = dih_idx[ii][3];
		switch (dih_type[ii])
		// switch for different dihedral type
		{
		case dih_none:
			break;
		case dih_opls_cosin:
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
			v1 = c1[ii];
			v2 = c2[ii];
			v3 = c3[ii];
			v4 = c4[ii];
			// calculate the energy
			vopls = v1+0.5*(v2*(1.0+cos(theta))+v3*(1.0-cos(2.0*theta))+v4*(1.0
					+cos(3.0*theta)));
			// add into the session energy
			gUdihSession += vopls;
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
		case dih_charmm: // future implementation 
			// modified from Shi Wei's code
			// NOT tested.
			// calculate the dihedral angle energy 
			// Note that add count_dih_multiple due to multiple matching parameters of the 
			// dihedral angle part
			// a-b-c-d  is i-j-k-l is 1-2-3-4
			// c1 is kphi, c2 is nperiod, c3 is delta0;
			kphi = c1[ii];
			nperiod = c2[ii];
			delta0 = c3[ii];
			// ** the following calculate the angle i-j-k-l 
			// ** between i-j-k plane and j-k-l plane
			if (fabs(kphi) >= 1e-5)
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
				gUdihSession += kphi*(1+cos((nperiod*phi-delta0)/180*pi));
				// ** calculate forces due to dihedral angle part
				// ** the following make sure that it is not divided by 0, if it is, it can
				// ** be reduced from 0/0 by mathematical tricks. Note that the 
				// ** delta0(ii) is eother 0 or 180 from the parameter data file.  
				if (fabs(sin(phi*pi/180)) <= 1e-10)
				{
					// write(*,*) ii
					// write(*,*) ii1,ii2,ii3,ii4
					// write(*,*) 'alpha=',acosd(alpha_temp)
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
	} // end of dihedral list loop

	return 0;
}

int impfrc(int iMole) // improper energy/force calculations
{
	int ii;
	int ii1, ii2, ii3, ii4;
	double rxjk, ryjk, rzjk;
	double rxij, ryij, rzij;
	double rxjl, ryjl, rzjl;
	double rxik, ryik, rzik;
	double rxkl, rykl, rzkl;
	double aa1, bb1, cc1;
	double aa2, bb2, cc2;
	double temp, temp1, temp2, temp3;
	double phi, alpha_temp;
	double fij;

	int iiStart, iiEnd;
	int iPhysMoleID;

	// Assign the lower and upper limits for ii loop
	if (iMole == -1)
	{
		iiStart = 0;
		iiEnd = nimp;
	}
	else if (iMole >= 0)
	{
		iiStart = mole_first_imp_idx[iMole];
		iiEnd = mole_last_imp_idx[iMole];
	}
	else
	{
		fprintf(stderr, "Error: Invalid parameters for impfrc(%d).\n", iMole);
		fprintf(fpouts, "Error: Invalid parameters for impfrc(%d).\n", iMole);
		exit(1);
	}

	// improper dihedral makes no contribution
	// to virial according to DL_POLY manual
	for (ii=iiStart; ii<iiEnd; ii++)
	{
		iPhysMoleID = imp2mole[ii];
		if (mole_status[iPhysMoleID] == MOLE_STATUS_VACANCY)
		{
			continue;
		}
		
		switch (imp_type[ii])
		{
		case imp_none:
			break;
		case imp_charmm:
			// following code is modified from shi wei's program
			// NOT tested.
			if (fabs(komega[ii]) >= 1e-5)
			{
				// ** the following calculate the angle i-j-k-l (A-B-C-D) 
				// ** where A is the centeral atom B, C, and D are bound to
				// ** improper angle is defined as the angle between 
				// ** the planes defined by ABC (ijk) and BCD (jkl) 
				// ** between i-j-k plane and j-k-l plane
				ii1 = imp_idx[ii][0];
				ii2 = imp_idx[ii][1];
				ii3 = imp_idx[ii][2];
				ii4 = imp_idx[ii][3];
				rxjk = xx[ii3]-xx[ii2];
				ryjk = yy[ii3]-yy[ii2];
				rzjk = zz[ii3]-zz[ii2];
				rxij = xx[ii1]-xx[ii2];
				ryij = yy[ii1]-yy[ii2];
				rzij = zz[ii1]-zz[ii2];
				// ** calculate a1,b1,c1 for the vector in x, y , and z direction.
				// ** This vector is perpendicular to the plane of i-j-K and obtained from
				// ** the cross product of jk and ij, ie. jk*ji
				// ** rotate from jk vector to ji vector
				aa1 = ryjk*rzij-rzjk*ryij;
				bb1 = rzjk*rxij-rxjk*rzij;
				cc1 = rxjk*ryij-rxij*ryjk;
				rxjl = xx[ii4]-xx[ii2];
				ryjl = yy[ii4]-yy[ii2];
				rzjl = zz[ii4]-zz[ii2];
				// ** calculate a2,b2, and c2 in x, y, z direction
				// ** this vector is the cross product of jl*jk 
				// ** rotate from jl to jk
				aa2 = ryjl*rzjk-rzjl*ryjk;
				bb2 = rzjl*rxjk-rxjl*rzjk;
				cc2 = rxjl*ryjk-rxjk*ryjl;
				temp1 = sqrt(aa1*aa1+bb1*bb1+cc1*cc1);
				temp2 = sqrt(aa2*aa2+bb2*bb2+cc2*cc2);
				temp3 = (aa1*aa2+bb1*bb2+cc1*cc2);
				phi = temp3/(temp1*temp2);
				alpha_temp = phi;
				// ** the angle between the two planes is below
				phi = 180 - (acos(phi)*180/pi);
				gUimpSession += komega[ii]*(phi-omega0[ii])*(phi-omega0[ii]);
				// ** the following is used to avoid dividing by 0
				// ** note that omega0[ii] is either 0 or 180 from the database
				// ** parameter file for improper
				if ((fabs(omega0[ii]) <= 1e-5)
						&& (fabs(acos(alpha_temp)*180/pi) <= 1e-20))
				{
					fprintf(stderr,"omega0 %lf and angle %lf do not match\n",omega0[ii],acos(alpha_temp)*180/pi);
					fprintf(fpouts, "omega0 %lf and angle %lf do not match\n",
							omega0[ii], acos(alpha_temp)*180/pi);
					exit(1);
				}
				if ((fabs(omega0[ii]-180) <= 1e-5) && (fabs(acos(alpha_temp)
						*180/pi-180) <= 1e-20))
				{
					fprintf(stderr,"omega0 %lf and angle %lf do not match\n",omega0[ii],acos(alpha_temp)*180/pi);
					fprintf(fpouts, "omega0 %lf and angle %lf do not match\n",
							omega0[ii], acos(alpha_temp)*180/pi);
					exit(1);
				}

				// ** calculate the forces due to improper 
				// ** the following is extremely important to use the 
				// ** mathematics tricks
				if ((fabs(omega0[ii]) <= 1e-5) && (fabs(acos(alpha_temp)*180/pi
						-180) <= 1e-20))
				{
					fprintf(stderr,"index=%d  %d-%d-%d-%d\n",ii,ii1,ii2,ii3,ii4);
					fprintf(stderr,"alpha=%lf\n",acos(alpha_temp)*180/pi);
					fprintf(fpouts, "index=%d  %d-%d-%d-%d\n", ii, ii1, ii2,
							ii3, ii4);
					fprintf(fpouts, "alpha=%lf\n", acos(alpha_temp)*180/pi);
					temp = 2.0*komega[ii]/alpha_temp;
					exit(1);
				}
				else if ((fabs(omega0[ii]-180) <= 1e-5)
						&& (fabs(acos(alpha_temp)*180/pi) <= 1e-20))
				{
					fprintf(stderr,"index=%d  %d-%d-%d-%d\n",ii,ii1,ii2,ii3,ii4);
					fprintf(stderr,"alpha=%lf\n",acos(alpha_temp)*180/pi);
					fprintf(fpouts, "index=%d  %d-%d-%d-%d\n", ii, ii1, ii2,
							ii3, ii4);
					fprintf(fpouts, "alpha=%lf\n", acos(alpha_temp)*180/pi);
					temp = 2.0*komega[ii]/alpha_temp;
					exit(1);
				}
				else
				{
					temp = -2.0*komega[ii]*(phi-omega0[ii])*pi/180/sin(phi*pi
							/180);
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
				rxik = xx[ii1]-xx[ii3];
				ryik = yy[ii1]-yy[ii3];
				rzik = zz[ii1]-zz[ii3];
				rxkl = xx[ii4]-xx[ii3];
				rykl = yy[ii4]-yy[ii3];
				rzkl = zz[ii4]-zz[ii3];
				fxs[ii2] = fxs[ii2] +fij*((bb2*rzik-bb1*rzkl-cc2*ryik+cc1*rykl)
						*temp1*temp2 -temp3*(temp2/temp1*(bb1*rzik-cc1*ryik)
						+temp1/temp2*(-bb2*rzkl+cc2*rykl)));
				fys[ii2] = fys[ii2] +fij
						*((-aa2*rzik+aa1*rzkl+cc2*rxik-cc1*rxkl)*temp1*temp2
								-temp3*(temp2/temp1*(-aa1*rzik+cc1*rxik)+temp1
										/temp2*(aa2*rzkl-cc2*rxkl)));
				fzs[ii2] = fzs[ii2] +fij*((aa2*ryik-aa1*rykl-bb2*rxik+bb1*rxkl)
						*temp1*temp2 -temp3*(temp2/temp1*(aa1*ryik-bb1*rxik)
						+temp1/temp2*(-aa2*rykl+bb2*rxkl)));
				// ** calculate the forces on atom 3
				fxs[ii3] = fxs[ii3] +fij
						*((-bb2*rzij+bb1*rzjl+cc2*ryij-cc1*ryjl)*temp1*temp2
								-temp3*(temp2/temp1*(-bb1*rzij+cc1*ryij)+temp1
										/temp2*(bb2*rzjl-cc2*ryjl)));
				fys[ii3] = fys[ii3] +fij*((aa2*rzij-aa1*rzjl-cc2*rxij+cc1*rxjl)
						*temp1*temp2 -temp3*(temp2/temp1*(aa1*rzij-cc1*rxij)
						+temp1/temp2*(-aa2*rzjl+cc2*rxjl)));
				fzs[ii3] = fzs[ii3] +fij
						*((-aa2*ryij+aa1*ryjl+bb2*rxij-bb1*rxjl)*temp1*temp2
								-temp3*(temp2/temp1*(-aa1*ryij+bb1*rxij)+temp1
										/temp2*(aa2*ryjl-bb2*rxjl)));
				// ** calculate the forces on atom 4
				fij = temp/(temp1*temp2*temp2);
				fxs[ii4]=fxs[ii4] + fij*((-bb1*rzjk+cc1*ryjk)*temp2-temp3/temp2
						*(-bb2*rzjk+cc2*ryjk));
				fys[ii4]=fys[ii4] + fij*((aa1*rzjk-cc1*rxjk)*temp2-temp3/temp2
						*(aa2*rzjk-cc2*rxjk));
				fzs[ii4]=fzs[ii4] + fij*((-aa1*ryjk+bb1*rxjk)*temp2-temp3/temp2
						*(-aa2*ryjk+bb2*rxjk));
			} // if (fabs(komega[ii]) >= 1e-5)
			break;
		default:
			printf("Error: unknown improper type.\n");
			exit(1);
		} // switch for improper types
	} // loop through all impropers

	if (imp_type[ii]==imp_charmm)
		gUimpSession *= (pi/180)*(pi/180);

	return 0;
}

int fnRafrcSession(int iMole)
{
	int ii;

	// energies
	gUintraSession = 0.0;
	gUbondSession = 0.0;
	gUangleSession = 0.0;
	gUdihSession = 0.0;
	gUimpSession = 0.0;

	// pressure related
	gVirialIntraSession = 0.0;

	for (ii=0; ii<natom_hist_max; ii++)
	{
		fxs[ii] = fys[ii] = fzs[ii] = 0.0;
	}

	if (nbond > 0)
	{
		bndfrc(iMole);
		gUintraSession += gUbondSession;
	}

	if (nangle > 0)
	{
		aglfrc(iMole);
		gUintraSession += gUangleSession;
	}

	if (ndih > 0)
	{
		dihfrc(iMole);
		gUintraSession += gUdihSession;
	}

	if (nimp > 0)
	{
		impfrc(iMole);
		gUintraSession += gUimpSession;
	}

	return 0;
}

int rafrc()
{
	// calculate the session energies
	fnRafrcSession(ENTIRE_SYSTEM);

	// energies
	uintra = gUintraSession;
	ubond = gUbondSession;
	uangle = gUangleSession;
	udih = gUdihSession;
	uimp = gUimpSession;

	// pressure related
	virial_intra = gVirialIntraSession;

	return 0;
}
