/**
 * Project: mspms2
 * File: impfrc.c
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

int impfrc(int iSpecie, int iImp, int iabs, double *uij, double *virial_ij)
{
	PSAMPLE_MOLECULE pSampleMole;
	int mm1, mm2, mm3, mm4;
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

	pSampleMole = sample_mole + iSpecie;

	mm1 = pSampleMole->imp_idx[iImp][0]; // relative atom id 1
	mm2 = pSampleMole->imp_idx[iImp][1]; // Relative aotm id 2
	mm3 = pSampleMole->imp_idx[iImp][2]; // Relative atom id 3
	mm4 = pSampleMole->imp_idx[iImp][3]; // relative atom id 4
	// The atom list of dihedrals must be ordered as 0-1-2-3.
	// Otherwise, the results won't be correct.
	ii1 = iabs + mm1; // absolute atom id 1
	ii2 = iabs + mm2; // absolute atom id 2
	ii3 = iabs + mm3; // absolute atom id 3
	ii4 = iabs + mm4; // absolute atom id 4
	switch (pSampleMole->imp_type[iImp])
	{
	case IMP_NONE:
		break;
	case IMP_CHARMM:
		// following code is modified from shi wei's program
		// NOT tested.
		if (fabs(pSampleMole->komega[iImp]) >= 1e-5)
		{
			// ** the following calculate the angle i-j-k-l (A-B-C-D) 
			// ** where A is the centeral atom B, C, and D are bound to
			// ** improper angle is defined as the angle between 
			// ** the planes defined by ABC (ijk) and BCD (jkl) 
			// ** between i-j-k plane and j-k-l plane
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
			// Calculate the energy
			*uij = pSampleMole->komega[iImp]*(phi-pSampleMole->omega0[iImp])*(phi-pSampleMole->omega0[iImp]);
			*uij *= (pi/180)*(pi/180);
			// Set virial to be zero
			*virial_ij = 0.0;
			// ** the following is used to avoid dividing by 0
			// ** note that omega0[ii] is either 0 or 180 from the database
			// ** parameter file for improper
			if ((fabs(pSampleMole->omega0[iImp])<=1e-5) && (fabs(acos(alpha_temp)*180/pi)<=1e-20))
			{
				fprintf(stderr,"omega0 %lf and angle %lf do not match\n",pSampleMole->omega0[iImp],
						acos(alpha_temp)*180/pi);
				fprintf(fpouts, "omega0 %lf and angle %lf do not match\n", pSampleMole->omega0[iImp], 
						acos(alpha_temp)*180/pi);
				exit(1);
			}
			if ((fabs(pSampleMole->omega0[iImp]-180)<=1e-5) && (fabs(acos(alpha_temp)*180/pi-180)<=1e-20))
			{
				fprintf(stderr,"omega0 %lf and angle %lf do not match\n",pSampleMole->omega0[iImp],
						acos(alpha_temp)*180/pi);
				fprintf(fpouts, "omega0 %lf and angle %lf do not match\n", 
						pSampleMole->omega0[iImp], acos(alpha_temp)*180/pi);
				exit(1);
			}

			// ** calculate the forces due to improper 
			// ** the following is extremely important to use the 
			// ** mathematics tricks
			if ((fabs(pSampleMole->omega0[iImp])<=1e-5) && (fabs(acos(alpha_temp)*180/pi-180)<=1e-20))
			{
				fprintf(stderr,"index=%d  %d-%d-%d-%d\n",iImp,ii1,ii2,ii3,ii4);
				fprintf(stderr,"alpha=%lf\n",acos(alpha_temp)*180/pi);
				fprintf(fpouts, "index=%d  %d-%d-%d-%d\n",iImp,ii1, ii2, ii3,ii4);
				fprintf(fpouts, "alpha=%lf\n", acos(alpha_temp)*180/pi);
				temp = 2.0*pSampleMole->komega[iImp]/alpha_temp;
				exit(1);
			}
			else if ((fabs(pSampleMole->omega0[iImp]-180) <= 1e-5) && (fabs(acos(alpha_temp)*180/pi) <= 1e-20))
			{
				fprintf(stderr,"index=%d  %d-%d-%d-%d\n",iImp,ii1,ii2,ii3,ii4);
				fprintf(stderr,"alpha=%lf\n",acos(alpha_temp)*180/pi);
				fprintf(fpouts, "index=%d  %d-%d-%d-%d\n", iImp, ii1, ii2, ii3,
						ii4);
				fprintf(fpouts, "alpha=%lf\n", acos(alpha_temp)*180/pi);
				temp = 2.0*pSampleMole->komega[iImp]/alpha_temp;
				exit(1);
			}
			else
			{
				temp = -2.0*pSampleMole->komega[iImp]*(phi-pSampleMole->omega0[iImp])*pi/180/sin(phi*pi/180);
			}
			// ** calculate forces on atom 1
			fij = temp/(temp2*temp1*temp1);
			fxs[ii1] = fxs[ii1] + fij*((bb2*rzjk-cc2*ryjk)*temp1-temp3 /temp1
					*(bb1*rzjk-cc1*ryjk));
			fys[ii1] = fys[ii1] + fij*((-aa2*rzjk+cc2*rxjk)*temp1-temp3 /temp1
					*(-aa1*rzjk+cc1*rxjk));
			fzs[ii1] = fzs[ii1] + fij*((aa2*ryjk-bb2*rxjk)*temp1-temp3 /temp1
					*(aa1*ryjk-bb1*rxjk));
			// ** calculate forces on atom 2
			fij = temp/(temp1*temp1*temp2*temp2);
			rxik = xx[ii1]-xx[ii3];
			ryik = yy[ii1]-yy[ii3];
			rzik = zz[ii1]-zz[ii3];
			rxkl = xx[ii4]-xx[ii3];
			rykl = yy[ii4]-yy[ii3];
			rzkl = zz[ii4]-zz[ii3];
			fxs[ii2] = fxs[ii2] +fij*((bb2*rzik-bb1*rzkl-cc2*ryik+cc1*rykl)
					*temp1*temp2 -temp3*(temp2/temp1*(bb1*rzik-cc1*ryik) +temp1
					/temp2*(-bb2*rzkl+cc2*rykl)));
			fys[ii2] = fys[ii2] +fij *((-aa2*rzik+aa1*rzkl+cc2*rxik-cc1*rxkl)
					*temp1*temp2 -temp3*(temp2/temp1*(-aa1*rzik+cc1*rxik)+temp1
					/temp2*(aa2*rzkl-cc2*rxkl)));
			fzs[ii2] = fzs[ii2] +fij*((aa2*ryik-aa1*rykl-bb2*rxik+bb1*rxkl)
					*temp1*temp2 -temp3*(temp2/temp1*(aa1*ryik-bb1*rxik) +temp1
					/temp2*(-aa2*rykl+bb2*rxkl)));
			// ** calculate the forces on atom 3
			fxs[ii3] = fxs[ii3] +fij *((-bb2*rzij+bb1*rzjl+cc2*ryij-cc1*ryjl)
					*temp1*temp2 -temp3*(temp2/temp1*(-bb1*rzij+cc1*ryij)+temp1
					/temp2*(bb2*rzjl-cc2*ryjl)));
			fys[ii3] = fys[ii3] +fij*((aa2*rzij-aa1*rzjl-cc2*rxij+cc1*rxjl)
					*temp1*temp2 -temp3*(temp2/temp1*(aa1*rzij-cc1*rxij) +temp1
					/temp2*(-aa2*rzjl+cc2*rxjl)));
			fzs[ii3] = fzs[ii3] +fij *((-aa2*ryij+aa1*ryjl+bb2*rxij-bb1*rxjl)
					*temp1*temp2 -temp3*(temp2/temp1*(-aa1*ryij+bb1*rxij)+temp1
					/temp2*(aa2*ryjl-bb2*rxjl)));
			// ** calculate the forces on atom 4
			fij = temp/(temp1*temp2*temp2);
			fxs[ii4]=fxs[ii4] + fij*((-bb1*rzjk+cc1*ryjk)*temp2-temp3/temp2
					*(-bb2*rzjk+cc2*ryjk));
			fys[ii4]=fys[ii4] + fij*((aa1*rzjk-cc1*rxjk)*temp2-temp3/temp2
					*(aa2*rzjk-cc2*rxjk));
			fzs[ii4]=fzs[ii4] + fij*((-aa1*ryjk+bb1*rxjk)*temp2-temp3/temp2
					*(-aa2*ryjk+bb2*rxjk));
		} // if (fabs(pSampleMole->komega[iImp]) >= 1e-5)
		break;
	default:
		printf("Error: unknown improper type.\n");
		exit(1);
	} // switch for improper types

	return 0;
}

