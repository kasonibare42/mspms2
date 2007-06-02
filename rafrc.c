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
    float gamma, fxa, fya, fza, fxc, fyc, fzc;

    float  vec_13_x, vec_13_y, vec_13_z;
    float  vec_12_x, vec_12_y, vec_12_z;
    float  vec_23_x, vec_23_y, vec_23_z;
    float  r12, r13, r23;
    float  delta_r1, delta_r2, delta_r3;
    float  k_theta, k_r_theta, k_r_rprime;
    float  u123;
    float  frc_term_1;
    float  frc_term_2_1, frc_term_2_2;
    float  frc_term_3, frc_term_3_1, frc_term_3_2;


    // angle bending calculations
    for (ii=0;ii<nangle;ii++)
    {
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
		if (sint<1.0e-8)
		    sint = 1.0e-8;
		// not sure why fmax wont work here
		// sint = fmaxf(0.1, sqrt(1.0-cost*cost)); // for 180 degree
		// angle difference
		delta_theta = theta - Thetaeq[ii];
		uangle_temp = 0.5*Ktheta[ii]*delta_theta*delta_theta;
		uangle += uangle_temp;
		// forces
		gamma = Ktheta[ii]*delta_theta/sint;
		fxa = gamma*(xbc-xab*cost)*rrab; // vector 1
		fya = gamma*(ybc-yab*cost)*rrab;
		fza = gamma*(zbc-zab*cost)*rrab;
		fxc = gamma*(xab-xbc*cost)*rrbc; // vector 2
		fyc = gamma*(yab-ybc*cost)*rrbc;
		fzc = gamma*(zab-zbc*cost)*rrbc;
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
		r13 = sqrt(vec_13_x*vec_13_x + vec_13_y*vec_13_y + vec_13_z*vec_13_z);
		// r_HO
		r12 = sqrt(vec_12_x*vec_12_x + vec_12_y*vec_12_y + vec_12_z*vec_12_z);
		// r_OH
		r23 = sqrt(vec_23_x*vec_23_x + vec_23_y*vec_23_y + vec_23_z*vec_23_z);
		delta_r3 = r13 - agl_para_4[ii]; // H-H stretch
		delta_r1 = r12 - agl_para_5[ii]; // O-H stretch
		delta_r2 = r23 - agl_para_5[ii]; // O-H stretch

		k_theta = Ktheta[ii];
	       	k_r_theta = Thetaeq[ii];
		k_r_rprime = agl_para_3[ii];
		u123 = 0.5*k_theta*delta_r3*delta_r3 + k_r_theta*delta_r3*(delta_r1+delta_r2) + k_r_rprime*delta_r1*delta_r2;
		uangle = uangle + u123;
		// forces
		// atom 1 (a)
		frc_term_1 = -k_theta*delta_r3/r13;
		frc_term_2_1 = -k_r_theta*delta_r3/r12;
		frc_term_2_2 = -k_r_theta*(delta_r1+delta_r2)/r13;
		frc_term_3 = -k_r_rprime*delta_r2/r12;
		fxs[iia] += (frc_term_1*vec_13_x + frc_term_2_1*vec_12_x + frc_term_2_2*vec_13_x + frc_term_3*vec_12_x);
		fys[iia] += (frc_term_1*vec_13_y + frc_term_2_1*vec_12_y + frc_term_2_2*vec_13_y + frc_term_3*vec_12_y);
		fzs[iia] += (frc_term_1*vec_13_z + frc_term_2_1*vec_12_z + frc_term_2_2*vec_13_z + frc_term_3*vec_12_z);

		// atom 2 (b)
		// frc_term_1 = 0.0;
		frc_term_2_1 = -k_r_theta*delta_r3/r23;
		frc_term_2_2 = k_r_theta*delta_r3/r12;
		frc_term_3_1 = -k_r_rprime*delta_r1/r23;
		frc_term_3_2 = k_r_rprime*delta_r2/r12;
		fxs[iib] += (frc_term_2_1*vec_23_x + frc_term_2_2*vec_12_x + frc_term_3_1*vec_23_x + frc_term_3_2*vec_12_x);
		fys[iib] += (frc_term_2_1*vec_23_y + frc_term_2_2*vec_12_y + frc_term_3_1*vec_23_y + frc_term_3_2*vec_12_y);
		fzs[iib] += (frc_term_2_1*vec_23_z + frc_term_2_2*vec_12_z + frc_term_3_1*vec_23_z + frc_term_3_2*vec_12_z);

		// atom 3 (c)
		frc_term_1 = k_theta*delta_r3/r13;
		frc_term_2_1 = k_r_theta*delta_r3/r23;
		frc_term_2_2 = k_r_theta*(delta_r1+delta_r2)/r13;
		frc_term_3 = k_r_rprime*delta_r1/r23;
		fxs[iic] += (frc_term_1*vec_13_x + frc_term_2_1*vec_23_x + frc_term_2_2*vec_13_x + frc_term_3*vec_23_x);
		fys[iic] += (frc_term_1*vec_13_y + frc_term_2_1*vec_23_y + frc_term_2_2*vec_13_y + frc_term_3*vec_23_y);
		fzs[iic] += (frc_term_1*vec_13_z + frc_term_2_1*vec_23_z + frc_term_2_2*vec_13_z + frc_term_3*vec_23_z);
		break;
	    default:
		printf("unknown angle type.\n");
		exit(1);
	} // ending of switch for different angle types
    } // end of angle bending calculations

}

int dihfrc()
{
    int ii;
    int iia, iib, iic, iid;
    float xdab, ydab, zdab, xdbc, ydbc, zdbc, xdcd, ydcd, zdcd;
    float xab, yab, zab, xbc, ybc, zbc, xcd, ycd, zcd, rrbc;
    float xac, yac, zac, pbx, pby, pbz, pb2, rpb1, rpb2;
    float pcx, pcy, pcz, pc2, rpc1, rpc2, pbpc, cost, sint;
    float theta, rsint, v1, v2, v3, v4, gamma;
    float fax, fay, faz, fcx, fcy, fcz;
    float fb1x, fb1y, fb1z, fd1x, fd1y, fd1z;
    float vopls, dopls;

    // dihedral torsion calculations
    for (ii=0;ii<ndih;ii++)
    {
	// The atom list of dihedrals must be ordered as 0-1-2-3
	// Otherwise, the results wont be correct
	iia = dih_idx[ii][0];
	iib = dih_idx[ii][1];
	iic = dih_idx[ii][2];
	iid = dih_idx[ii][3];
	switch (dih_type[ii]) // switch for different dihedral type
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
		sint = (xbc*(pcy*pbz-pcz*pby)+ybc*(pbx*pcz-pbz*pcx)+zbc*(pcx*pby-pcy*pbx))
		    *(rpb1*rpc1*rrbc);

		theta = atan2(sint, cost);

		// avoid singularity in sint
		sint = copysign(fmax(1.0e-8,fabs(sint)), sint);
		rsint = 1.0/sint;

		// phase angle? theta = theta - prmdih(kk,7) ?
		v1 = c1[ii];
		v2 = c2[ii];
		v3 = c3[ii];
		v4 = c4[ii];
		// calculate the energy
		vopls = v1+0.5*(v2*(1.0+cos(theta))+v3*(1.0-cos(2.0*theta))+v4*(1.0+cos(3.0*theta)));
		udih = udih + vopls;
		// forces
		dopls = -0.5*(v2*sin(theta)-2.0*v3*sin(2.0*theta)+3.0*v4*sin(3.0*theta));
		gamma = dopls*rpb1*rpc1*rsint;

		fax = gamma*((-pcy*zbc+pcz*ybc)-pbpc*rpb2*(-pby*zbc+pbz*ybc));
		fay = gamma*(( pcx*zbc-pcz*xbc)-pbpc*rpb2*( pbx*zbc-pbz*xbc));
		faz = gamma*((-pcx*ybc+pcy*xbc)-pbpc*rpb2*(-pbx*ybc+pby*xbc));

		fcx = gamma*((-pcy*zab+pcz*yab)-pbpc*rpb2*(-pby*zab+pbz*yab));
		fcy = gamma*(( pcx*zab-pcz*xab)-pbpc*rpb2*( pbx*zab-pbz*xab));
		fcz = gamma*((-pcx*yab+pcy*xab)-pbpc*rpb2*(-pbx*yab+pby*xab));

		fb1x= gamma*((-pby*zcd+pbz*ycd)-pbpc*rpc2*(-pcy*zcd+pcz*ycd));
		fb1y= gamma*(( pbx*zcd-pbz*xcd)-pbpc*rpc2*( pcx*zcd-pcz*xcd));
		fb1z= gamma*((-pbx*ycd+pby*xcd)-pbpc*rpc2*(-pcx*ycd+pcy*xcd));

		fd1x= gamma*((-pby*zbc+pbz*ybc)-pbpc*rpc2*(-pcy*zbc+pcz*ybc));
		fd1y= gamma*(( pbx*zbc-pbz*xbc)-pbpc*rpc2*( pcx*zbc-pcz*xbc));
		fd1z= gamma*((-pbx*ybc+pby*xbc)-pbpc*rpc2*(-pcx*ybc+pcy*xbc));

		fxs[iia] += fax;
		fys[iia] += fay;
		fzs[iia] += faz;

		fxs[iib] += (-fax-fcx+fb1x);
		fys[iib] += (-fay-fcy+fb1y);
		fzs[iib] += (-faz-fcz+fb1z);

		fxs[iic] +=  fcx-fb1x-fd1x;
		fys[iic] +=  fcy-fb1y-fd1y;
		fzs[iic] +=  fcz-fb1z-fd1z;

		fxs[iid] +=  fd1x;
		fys[iid] +=  fd1y;
		fzs[iid] +=  fd1z;
		break;
	    case dih_charmm: // future implementation
		break;
	    default:
		printf("unknown dihedral type.\n");
		exit(1);
	} // switch for different dihedral type
    } // end of dihedral list loop
}

int impfrc() // improper energy/force calculations
{}

int rafrc()
{
    int ii;

    uintra = 0.0;
    ubond = 0.0;
    uangle = 0.0;
    udih = 0.0;
    uimp = 0.0;

    for (ii=0;ii<natom;ii++)
	fxs[ii] = fys[ii] = fzs[ii] = 0.0;

    // if (nbond > 0)
// 	bndfrc();

    if (nangle > 0)
	aglfrc();

  //   if (ndih > 0)
//	dihfrc();

    if (nimp > 0)
	impfrc();

    // printf("%f %f %f %f\n", ubond, uangle, udih, uimp);

    // total intra energy
    uintra = ubond + uangle + udih + uimp;
}

