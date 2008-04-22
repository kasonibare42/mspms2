/**
 * Project: mspms2
 * File: aglfrc.c
 * 
 * Copyright (C) 2008    Yang Wang <ywangd@gmail.com>
 * Created @ Apr 22, 2008
 * 
 * Description:
 * 
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include "mspms2.h"

int aglfrc(int iSpecie, int iAngle, int iabs, double *uij, double *virial_ij)
{
	int mm1, mm2, mm3, iia, iib, iic;
	PSAMPLE_MOLECULE pSampleMole;
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
	
	pSampleMole = sample_mole + iSpecie;
	
	mm1 = pSampleMole->agl_idx[iAngle][0]; // relative atom id 1
	mm2 = pSampleMole->agl_idx[iAngle][1]; // Relative atom id 2
	mm3 = pSampleMole->agl_idx[iAngle][2]; // relative atom id 3
	// Calculate absolute atom id
	// The atom list of the angle must be list as 0-1-2.
	// Otherwise, the calculations are not right.
	iia = iabs + mm1; // absolute atom id 1
	iib = iabs + mm2; // absolute atom id 2
	iic = iabs + mm3; // absolute atom id 3
	switch (pSampleMole->agl_type[iAngle])
	{
	case ANGLE_NONE:
		break;
	case ANGLE_HARMONIC:
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
		delta_theta = theta - pSampleMole->Thetaeq[iAngle];
		*uij = 0.5*pSampleMole->Ktheta[iAngle]*delta_theta*delta_theta;
		// forces
		gamma = pSampleMole->Ktheta[iAngle]*delta_theta/sint;
		fxa = gamma*(xbc-xab*cost)*rrab; // vector 1
		fya = gamma*(ybc-yab*cost)*rrab;
		fza = gamma*(zbc-zab*cost)*rrab;
		fxc = gamma*(xab-xbc*cost)*rrbc; // vector 2
		fyc = gamma*(yab-ybc*cost)*rrbc;
		fzc = gamma*(zab-zbc*cost)*rrbc;
		// For angle potentials without r as variable, the net contribution to virial is zero.
		*virial_ij = 0.0; 
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
	case ANGLE_TR_WATER: // Toukan & Rhaman water potential
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
		delta_r3 = r13 - pSampleMole->agl_para_4[iAngle]; // H-H stretch
		delta_r1 = r12 - pSampleMole->agl_para_5[iAngle]; // O-H stretch
		delta_r2 = r23 - pSampleMole->agl_para_5[iAngle]; // O-H stretch

		k_theta = pSampleMole->Ktheta[iAngle];
		k_r_theta = pSampleMole->Thetaeq[iAngle];
		k_r_rprime = pSampleMole->agl_para_3[iAngle];
		*uij = 0.5*k_theta*delta_r3*delta_r3 + k_r_theta*delta_r3*(delta_r1
				+delta_r2) + k_r_rprime*delta_r1*delta_r2;
		// Forces
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
		// Contribution to virial. Not sure. Need double check
		*virial_ij = fxa*vec_12_x + fya*vec_12_y + fza*vec_12_z;
		*virial_ij = *virial_ij - fxc*vec_23_x - fyc*vec_23_y - fzc*vec_23_z;
		break;
	default:
		printf("unknown angle type.\n");
		exit(1);
	} // ending of switch for different angle types
	
	return 0;
}

