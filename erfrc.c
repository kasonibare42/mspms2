/* Calculate inter atom and inter molecule energies and forces.
 These include different atoms on the same molecule (exclude 1-2,1-3,1-4)
 1-4 will be still be calculated after
 and atoms on different molecules.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_deriv.h>
#include "vars.h"

extern int sffrc();

double deriv_inc_gamma(double x, void *params)
{
	return gsl_sf_gamma_inc(0.0, x);
}

// calculate the center of mass for all the molecules
// and the inner coordinates
int cal_com_and_inner_coords()
{
	int ii, jj;

	for (ii=0; ii<nmole; ii++)
	{
		mole_xx[ii] = 0.0;
		mole_yy[ii] = 0.0;
		mole_zz[ii] = 0.0;
		// calculate the center of mass
		for (jj=mole_first_atom_idx[ii]; jj<mole_first_atom_idx[ii+1]; jj++)
		{
			mole_xx[ii] += xx[jj]*aw[jj];
			mole_yy[ii] += yy[jj]*aw[jj];
			mole_zz[ii] += zz[jj]*aw[jj];
		}
		mole_xx[ii] /= mw[ii];
		mole_yy[ii] /= mw[ii];
		mole_zz[ii] /= mw[ii];
		// calculate the relative position of the atoms to the center of mass
		for (jj=mole_first_atom_idx[ii]; jj<mole_first_atom_idx[ii+1]; jj++)
		{
			ex[jj] = xx[jj] - mole_xx[ii];
			fy[jj] = yy[jj] - mole_yy[ii];
			gz[jj] = zz[jj] - mole_zz[ii];
		}
	}
	return 0;
}

// calculate the atom position according to the PBC'd center of mass
// the coordinates are saved to ex, fy, gz
int reconstruct_coords_from_com()
{
	int ii, jj;

	for (ii=0; ii<nmole; ii++)
	{
		for (jj=mole_first_atom_idx[ii]; jj<mole_first_atom_idx[ii+1]; jj++)
		{
			ex[jj] = ex[jj] + mole_xx[ii];
			fy[jj] = fy[jj] + mole_yy[ii];
			gz[jj] = gz[jj] + mole_zz[ii];
		}
	}
	return 0;
}

// calculate the PBC'd atom coordinates
// the new coordinates are saved in ex,fy,gz
int cal_coords_PBC()
{
	int ii;
	for (ii=0; ii<natom; ii++)
		;
	{
		ex[ii] = xx[ii] - boxlx*rint(xx[ii]/boxlx);
		fy[ii] = yy[ii] - boxlx*rint(yy[ii]/boxlx);
		gz[ii] = zz[ii] - boxlx*rint(zz[ii]/boxlx);
	}
	return 0;
}

// calcualte interaction between atom ii and jj, exclude those from excluding list
int loop_ij()
{
	int ii, jj, kk;
	double xxi, yyi, zzi;
	double fxi, fyi, fzi;
	double rxij, ryij, rzij;
	double rijsq, rij, r_rijsq;
	double r_r6, r_r12, r_r12_minus_r_r6;
	double epsiloni, sigmai, chargei;
	double epsilonij, sigmaij;
	int isNotexcl[natom_max];
	double uij_vdw, uij_vdw_temp, uij_real, uij_real_temp;
	double uij_wolf_temp;
	double fij, fxij, fyij, fzij;
	double temp1, temp2;
	double rhoklsq, kapparhokl_sq;
	double uij_Gz0; // 1D ewald
	gsl_function FF;

	FF.function = &deriv_inc_gamma;
	FF.params = 0;

	uij_vdw = 0.0;
	uij_real = 0.0;
	uij_Gz0 = 0.0;

	// atom ii
	for (ii=0; ii<natom-1; ii++)
	{
		if (isghost[ii] == all_ghost) // move to ii+1 if atom ii is full ghost atom
			continue;
		xxi = xx[ii];
		yyi = yy[ii];
		zzi = zz[ii];
		fxi = fxl[ii];
		fyi = fyl[ii];
		fzi = fzl[ii];
		epsiloni = epsilon[ii];
		sigmai = sigma[ii];
		chargei = charge[ii];
		// set the exclude list for atom ii
		for (kk=0; kk<natom; kk++)
			isNotexcl[kk] = true;
		for (kk=pointexcl[ii]; kk<pointexcl[ii+1]; kk++)
			isNotexcl[excllist[kk]] = false;
		// atom jj
		for (jj=ii+1; jj<natom; jj++)
		{
			if (isghost[jj] == all_ghost) // move to jj+1 if atom jj is full ghost atom
				continue;
			if (isNotexcl[jj]) // calculate interactions if it is not exclusion pair
			{
				rxij = xxi - xx[jj];
				ryij = yyi - yy[jj];
				rzij = zzi - zz[jj];

				// calculate Gz=0 term for 1D ewald summation
				// this has to be done before the minimum image convention
				// not sure if cutoff check should be applied
				// If the cutoff is larger than the size in x,y dimension,
				// this should not be a problem
				if (isEwaldOn && fEwald_Dim==ewald_1D)
				{
					rhoklsq = rxij*rxij + ryij*ryij;
					if (rhoklsq>parallel_to_z_err)
					{
						kapparhokl_sq = kappasq*rhoklsq;
						temp1 = Euler_const + gsl_sf_gamma_inc(0.0,kapparhokl_sq) + log(kapparhokl_sq);
						uij_Gz0 = uij_Gz0 -chargei*charge[jj]*temp1;
						// forces
						// derivative of incomplete gamma function
						// use numerical differential
						// This may not be a good way to do it
						// temp1 is the value, temp2 is the absolute std. err.
						gsl_deriv_central(&FF, kapparhokl_sq, 1.0e-8, &temp1,
								&temp2);
						temp2 = kappasq*temp1 + 1.0/rhoklsq;
						fij = chargei*charge[jj]*temp2*const_columb/boxlz;
						fxij = fij*rxij;
						fyij = fij*ryij;
						// force on atom ii
						fxi += fxij;
						fyi += fyij;
						// force on atom jj
						fxl[jj] -= fxij;
						fyl[jj] -= fyij;
					} // if rhoklsq != 0.0 (>1.0e-10)
				} // if 1D ewald is used

				// minimum image convention
				rxij = rxij - boxlx*rint(rxij/boxlx);
				ryij = ryij - boxly*rint(ryij/boxly);
				rzij = rzij - boxlz*rint(rzij/boxlz);
				rijsq = rxij*rxij + ryij*ryij + rzij*rzij;
				rij = sqrt(rijsq);

				// LJ part, also does ghost atom check
				if (isghost[ii]!=lj_ghost && isghost[jj]!=lj_ghost && rijsq
						<rcutoffsq)
				{
					// if switch potential for LJ is on, then calculate the switch
					if (isLJswitchOn)
					{
						if (rijsq<rcutonsq)
							LJswitch = 1.0;
						else
							LJswitch = (rcutoffsq-rijsq)*(rcutoffsq-rijsq)
									*(rcutoffsq+2.0*rijsq-3.0*rcutonsq)
									/roff2_minus_ron2_cube;
					}
					sigmaij = 0.5*(sigmai+sigma[jj]);
					epsilonij = sqrt(epsiloni*epsilon[jj]);
					r_rijsq = sigmaij*sigmaij/rijsq;
					r_r6 = r_rijsq*r_rijsq*r_rijsq;
					r_r12 = r_r6*r_r6;
					r_r12_minus_r_r6 = r_r12 - r_r6;
					uij_vdw_temp = epsilonij*r_r12_minus_r_r6; // still need *4.0
					if (isLJswitchOn) // if switch is used
						uij_vdw += uij_vdw_temp*LJswitch; // still need 4.0
					else
						uij_vdw += uij_vdw_temp; // still need *4.0
					// force calculations
					if (isLJswitchOn)
					{
						if (rijsq<rcutonsq)
							fij = 24.0*epsilonij*(r_r12_minus_r_r6+r_r12)/rijsq;
						else
							fij = 24.0*epsilonij*(r_r12_minus_r_r6+r_r12)/rijsq
									*LJswitch // 4.0 for the real energy
									-4.0*uij_vdw_temp*12.0*(rcutoffsq-rijsq)
											*(rcutonsq-rijsq)
											/roff2_minus_ron2_cube;
					}
					else
						fij = 24.0*epsilonij*(r_r12_minus_r_r6+r_r12)/rijsq;
					fxij = fij*rxij;
					fyij = fij*ryij;
					fzij = fij*rzij;
					// contribution to the virial
					virial_inter = virial_inter + fxij*rxij + fyij*ryij + fzij
							*rzij;
					// force on atom ii
					fxi += fxij;
					fyi += fyij;
					fzi += fzij;
					// force on atom jj
					fxl[jj] -= fxij;
					fyl[jj] -= fyij;
					fzl[jj] -= fzij;
				}

				// electrostatic part
				if (isEwaldOn && rijsq<rcutoffelecsq)
				{
					uij_real_temp = chargei*charge[jj]/rij;
					uij_real += uij_real_temp*erfc(kappa*rij); // real part ewald energy, still need 1/4*pi*epsilon0
					temp1 = kappa*rij;
					temp2 = uij_real_temp*erfc(temp1);
					// real part force calculation
					fij = (temp2 + uij_real_temp*2.0*temp1*exp(-temp1*temp1)
							/sqrt(pi))*const_columb/rijsq;
					fxij = fij*rxij;
					fyij = fij*ryij;
					fzij = fij*rzij;
					// force on atom ii
					fxi += fxij;
					fyi += fyij;
					fzi += fzij;
					// force on atom jj
					fxl[jj] -= fxij;
					fyl[jj] -= fyij;
					fzl[jj] -= fzij;
				} // is Charge On check and rcutoffelecsq check

				// wolf method for electrostatic
				if (isWolfOn && rijsq<rcutoffelecsq)
				{
					uwolf_real = uwolf_real + chargei*charge[jj]
							*(erfc(kappa*rij)/rij + wolfvcon1 + wolfvcon2*(rij
									-rcutoffelec));
					uij_wolf_temp = chargei*charge[jj]/rij;
					fij = const_columb*uij_wolf_temp
					*(erfc(kappa*rij)/rijsq + wolffcon1*exp(-(kappa*rij)*(kappa*rij))/rij + wolffcon2);
					fxij = fij*rxij;
					fyij = fij*ryij;
					fzij = fij*rzij;
					// contribution to the virial
					virial_inter = virial_inter + fxij*rxij + fyij*ryij + fzij
							*rzij;
					// force on atom ii
					fxi += fxij;
					fyi += fyij;
					fzi += fzij;
					// force on atom jj
					fxl[jj] -= fxij;
					fyl[jj] -= fyij;
					fzl[jj] -= fzij;
				}
			} // end of excluding list
		} // end of atom jj
		fxl[ii] = fxi;
		fyl[ii] = fyi;
		fzl[ii] = fzi;
	} // end of atom ii

	uvdw += uij_vdw; // still need 4.0
	ureal += uij_real; // still neeed constant
	uGz0 += uij_Gz0/(2.0*boxlz); // still need constant
	return 0;
}

// calculate the interactions between dihedral ending atom 1 and 4
int loop_14()
{
	int ii;
	int ii1, ii2;
	double rxij, ryij, rzij;
	double rijsq, rij, r_rijsq;
	double r_r6, r_r12, r_r12_minus_r_r6;
	double sigmaij, epsilonij;
	double uij_vdw14, uij_vdw14_temp, uij_real14, uij_real14_temp;
	double uij_wolf14_temp;
	double rijsq_old;
	int fCalculate14_LJ, fCalculate14_elec;
	double fij, fxij, fyij, fzij;
	double temp1, temp2;
	double rhoklsq, kapparhokl_sq;
	double uij_Gz0; // 1D ewald
	double uij_coulomb14_temp;
	gsl_function FF;

	FF.function = &deriv_inc_gamma;
	FF.params = 0;

	uij_vdw14 = 0.0;
	uij_real14 = 0.0;
	uij_Gz0 = 0.0;

	for (ii=0; ii<ndih; ii++)
	{
		if (isDih_unique[ii]) // only calculate interactions for unique 1,4 pairs
		{
			ii1 = dih_idx[ii][0];
			ii2 = dih_idx[ii][3];
			// no ghost check needed here since 1,4 are in the same molecule
			// and ghost check is only for inter molecules or solid-fluid
			rxij = xx[ii1] - xx[ii2];
			ryij = yy[ii1] - yy[ii2];
			rzij = zz[ii1] - zz[ii2];

			// Gz=0 term for 1D ewald
			if (isEwaldOn && fEwald_Dim==ewald_1D)
			{
				rhoklsq = rxij*rxij + ryij*ryij;
				if (rhoklsq>parallel_to_z_err)
				{
					kapparhokl_sq = kappasq*rhoklsq;
					temp1 = Euler_const + gsl_sf_gamma_inc(0.0,kapparhokl_sq) + log(kapparhokl_sq);
					uij_Gz0 = uij_Gz0 -charge[ii1]*charge[ii2]*temp1;
					// forces
					// derivative of incomplete gamma function
					// use numerical differential
					// This may not be a good way to do it
					// temp1 is the value, temp2 is the absolute std. err.
					gsl_deriv_central(&FF, kapparhokl_sq, 1.0e-8, &temp1,
							&temp2);
					temp2 = kappasq*temp1 + 1.0/rhoklsq;
					fij = charge[ii1]*charge[ii2]*temp2*const_columb/boxlz;
					fxij = fij*rxij;
					fyij = fij*ryij;
					// force on atom ii1
					fxl[ii1] += fxij;
					fyl[ii1] += fyij;
					// force on atom ii2
					fxl[ii2] -= fxij;
					fyl[ii2] -= fyij;
				} // if rhoklsq != 0.0 
			} // if 1D ewald is used

			// check for necessary 1,4 parameter modifiers
			// and calculate the corresponding sigmaij and epsilonij
			// OPLS modifier
			sigmaij = 0.5*(sigma[ii1]+sigma[ii2]);
			epsilonij = f0*sqrt(epsilon[ii1]*epsilon[ii2]);
			// can add CHARMM modifier here, different 1,4 parameters
			// must be done before minimum image convention

			// save old rxij, ryij and rzij to check if 1,4 are in same molecule
			// or 1,4' image interaction
			rijsq_old = rxij*rxij + ryij*ryij + rzij*rzij;
			// minimum image convention
			rxij = rxij - boxlx*rint(rxij/boxlx);
			ryij = ryij - boxly*rint(ryij/boxly);
			rzij = rzij - boxlz*rint(rzij/boxlz);
			rijsq = rxij*rxij + ryij*ryij + rzij*rzij;
			rij = sqrt(rijsq);
			// if 1,4 distance is changed after minimum image convention
			// which means 1,4' are used now and they are in two different molecules
			// now ghost check should be applied
			if (fabs(rijsq_old-rijsq)>tolerant_err) // 1,4' found, need check ghost before calculate LJ for them
			{
				fprintf(stderr,"Warning: long 1,4 non-bonded pair %d-%d found...\n",ii1,ii2);
				fprintf(fpouts,
						"Warning: long 1,4 non-bonded pair %d-%d found...\n",
						ii1, ii2);
				// ghost check for 1,4'
				if (isghost[ii1]==all_ghost || isghost[ii2]==all_ghost)
				{
					// full ghost
					fCalculate14_LJ = false;
					fCalculate14_elec = false;
				}
				else if (isghost[ii1]==lj_ghost || isghost[ii2]==lj_ghost)
				{
					// LJ ghost
					fCalculate14_LJ = false;
					fCalculate14_elec = true;
				}
				else
				{
					// not ghost
					fCalculate14_LJ = true;
					fCalculate14_elec = true;
				}
			}
			else // 1,4 found, always calculate interactions for them, in same molecule
			{
				fCalculate14_LJ = true;
				fCalculate14_elec = true;
			}

			// LJ part, 
			if (fCalculate14_LJ==true && rijsq<rcutoffsq)
			{
				if (isLJswitchOn) // if use switch for LJ
				{
					if (rijsq<rcutonsq)
						LJswitch = 1.0;
					else
						LJswitch = (rcutoffsq-rijsq)*(rcutoffsq-rijsq)
								*(rcutoffsq+2.0*rijsq-3.0*rcutonsq)
								/roff2_minus_ron2_cube;
				}
				r_rijsq = sigmaij*sigmaij/rijsq;
				r_r6 = r_rijsq*r_rijsq*r_rijsq;
				r_r12 = r_r6*r_r6;
				r_r12_minus_r_r6 = r_r12 - r_r6;
				uij_vdw14_temp = epsilonij*r_r12_minus_r_r6; // still need *4.0
				if (isLJswitchOn)
					uij_vdw14 += uij_vdw14_temp*LJswitch; // still need 4.0
				else
					uij_vdw14 += uij_vdw14_temp; // still need *4.0
				// force calculations
				if (isLJswitchOn)
				{
					if (rijsq<rcutonsq)
						fij = 24.0*epsilonij*(r_r12_minus_r_r6+r_r12)/rijsq;
					else
						fij = 24.0*epsilonij*(r_r12_minus_r_r6+r_r12)/rijsq
								*LJswitch // 4.0 for the real energy
								-4.0*uij_vdw14_temp*12.0*(rcutoffsq-rijsq)
										*(rcutonsq-rijsq)/roff2_minus_ron2_cube;
				}
				else
					fij = 24.0*epsilonij*(r_r12_minus_r_r6+r_r12)/rijsq;
				fxij = fij*rxij;
				fyij = fij*ryij;
				fzij = fij*rzij;
				// contribution to the virial
				virial_inter = virial_inter + fxij*rxij + fyij*ryij + fzij*rzij;
				// force on atom ii1
				fxl[ii1] += fxij;
				fyl[ii1] += fyij;
				fzl[ii1] += fzij;
				// force on atom ii2
				fxl[ii2] -= fxij;
				fyl[ii2] -= fyij;
				fzl[ii2] -= fzij;
			} // rcutoffsq

			// electrostatic part
			if (fCalculate14_elec && isEwaldOn && rijsq<rcutoffelecsq)
			{
				uij_real14_temp = charge[ii1]*charge[ii2]/rij;
				uij_real14 += uij_real14_temp*erfc(kappa*rij); // real part ewald energy, still need 1/4*pi*epsilon0
				temp1 = kappa*rij;
				temp2 = uij_real14_temp*erfc(temp1);
				// real part force calculation
				fij = (temp2 + uij_real14_temp*2.0*temp1*exp(-temp1*temp1)
						/sqrt(pi))*const_columb/rijsq;
				fxij = fij*rxij;
				fyij = fij*ryij;
				fzij = fij*rzij;
				// force on atom ii
				fxl[ii1] += fxij;
				fyl[ii1] += fyij;
				fzl[ii1] += fzij;
				// force on atom jj
				fxl[ii2] -= fxij;
				fyl[ii2] -= fyij;
				fzl[ii2] -= fzij;
			} // is charge on check and rcutoffelecsq

			// wolf method for electrostatic
			if (fCalculate14_elec && isWolfOn && rijsq<rcutoffelecsq)
			{
				uwolf_real = uwolf_real + charge[ii1]*charge[ii2]*(erfc(kappa
						*rij)/rij + wolfvcon1 + wolfvcon2*(rij-rcutoffelec));
				// need + some scale factor here
				uij_wolf14_temp = charge[ii1]*charge[ii2]/rij;
				fij = const_columb*uij_wolf14_temp
				*(erfc(kappa*rij)/rijsq + wolffcon1*exp(-(kappa*rij)*(kappa*rij))/rij + wolffcon2);
				// need + some scale factor here 
				fxij = fij*rxij;
				fyij = fij*ryij;
				fzij = fij*rzij;
				// contribution to the virial
				virial_inter = virial_inter + fxij*rxij + fyij*ryij + fzij*rzij;
				// force on atom ii1
				fxl[ii1] += fxij;
				fyl[ii1] += fyij;
				fzl[ii1] += fzij;
				// force on atom ii2
				fxl[ii2] -= fxij;
				fyl[ii2] -= fyij;
				fzl[ii2] -= fzij;
			}

			// simple coulomb 
			if (fCalculate14_elec && isSimpleCoulomb && rijsq<rcutoffelecsq)
			{
				uij_coulomb14_temp = charge[ii1]*charge[ii2]/rij; // need constant
				ucoulomb += uij_coulomb14_temp; // need constant
				fij = const_columb*uij_coulomb14_temp/rijsq;
				fxij = fij*rxij;
				fyij = fij*ryij;
				fzij = fij*rzij;
				// contribution to the virial
				virial_inter = virial_inter + fxij*rxij + fyij*ryij + fzij*rzij;
				// force on atom ii1
				fxl[ii1] += fxij;
				fyl[ii1] += fyij;
				fzl[ii1] += fzij;
				// force on atom ii2
				fxl[ii2] -= fxij;
				fyl[ii2] -= fyij;
				fzl[ii2] -= fzij;
			}

		} // unique 1,4 pair check
	} // loop through dihedrals

	// add into total energy
	uvdw += uij_vdw14; // still need 4.0
	ureal += uij_real14; // still need constant
	uGz0 += uij_Gz0/(2.0*boxlz); // still need constant

	return 0;
}

int loop_13()
{
	int ii;
	int ii1, ii2;
	double rxij, ryij, rzij;
	double rijsq, rij, r_rijsq;
	double r_r6, r_r12, r_r12_minus_r_r6;
	double sigmaij, epsilonij;
	double fij, fxij, fyij, fzij;
	double uij_vdw13img, uij_vdw13img_temp;
	double uij_real13;
	double uij_excl_13;
	double rxij_old, ryij_old, rzij_old;
	double temp1, temp2, temp3;
	double rhoklsq, kapparhokl_sq;
	double uij_Gz0; // 1D ewald
	gsl_function FF;

	FF.function = &deriv_inc_gamma;
	FF.params = 0;

	uij_excl_13 = 0.0;
	uij_vdw13img = 0.0;
	uij_real13 = 0.0;
	uij_Gz0 = 0.0;

	for (ii=0; ii<nangle; ii++)
	{
		if (isAngle_unique) // only calculate interactions for unique 1,3 pairs
		{
			ii1 = angle_idx[ii][0];
			ii2 = angle_idx[ii][2];
			// do not need check ghost atom since 1,3 is in the same molecule
			rxij = xx[ii1] - xx[ii2];
			ryij = yy[ii1] - yy[ii2];
			rzij = zz[ii1] - zz[ii2];

			// Gz=0 term for 1D ewald
			if (isEwaldOn && fEwald_Dim==ewald_1D)
			{
				rhoklsq = rxij*rxij + ryij*ryij;
				if (rhoklsq>parallel_to_z_err)
				{
					kapparhokl_sq = kappasq*rhoklsq;
					temp1 = Euler_const + gsl_sf_gamma_inc(0.0,kapparhokl_sq) + log(kapparhokl_sq);
					uij_Gz0 = uij_Gz0 -charge[ii1]*charge[ii2]*temp1;
					// forces
					// derivative of incomplete gamma function
					// use numerical differential
					// This may not be a good way to do it
					// temp1 is the value, temp2 is the absolute std. err.
					gsl_deriv_central(&FF, kapparhokl_sq, 1.0e-8, &temp1,
							&temp2);
					temp2 = kappasq*temp1 + 1.0/rhoklsq;
					fij = charge[ii1]*charge[ii2]*temp2*const_columb/boxlz;
					fxij = fij*rxij;
					fyij = fij*ryij;
					// force on atom ii1
					fxl[ii1] += fxij;
					fyl[ii1] += fyij;
					// force on atom ii2
					fxl[ii2] -= fxij;
					fyl[ii2] -= fyij;
				} // if rhoklsq != 0.0 
			} // if 1D ewald is used


			// save the old positions
			rxij_old = rxij;
			ryij_old = ryij;
			rzij_old = rzij;

			if (isEwaldOn || isWolfOn) // if we need charge interactions
			{
				// calcualte the exculde ewald part, which use the raw distance between ii,jj
				// so no minimum image convention here
				rijsq = rxij*rxij + ryij*ryij + rzij*rzij;
				rij = sqrt(rijsq);
				temp2 = charge[ii1]*charge[ii2]/rij; // still need constant
				uij_excl_13 += temp2; // still need constant
				// forces
				fij = -const_columb*temp2/rijsq; // NOTE: the negative sign here
				// The above negative sign is because these forces are meant to be
				// substracted from the total forces. Not additive.
				fxij = fij*rxij;
				fyij = fij*ryij;
				fzij = fij*rzij;
				if (isWolfOn) // negative contribution to the virial
					virial_inter = virial_inter + fxij*rxij + fyij*ryij + fzij
							*rzij;
				// forces on ii1
				fxl[ii1] += fxij;
				fyl[ii1] += fyij;
				fzl[ii1] += fzij;
				// forces on ii2
				fxl[ii2] -= fxij;
				fyl[ii2] -= fyij;
				fzl[ii2] -= fzij;

				// calculate the real part of ewald summation
				// minimum image convention now
				rxij = rxij - boxlx*rint(rxij/boxlx);
				ryij = ryij - boxly*rint(ryij/boxly);
				rzij = rzij - boxlz*rint(rzij/boxlz);
				rijsq = rxij*rxij + ryij*ryij + rzij*rzij;
				if (rijsq<rcutoffelecsq) // elec cutoff
				{
					rij = sqrt(rijsq);
					temp1 = kappa*rij;
					temp2 = charge[ii1]*charge[ii2]/rij;
					temp3 = temp2*erfc(temp1);
					if (isEwaldOn)
					{
						uij_real13 += temp3; // still need constant
						// forces
						fij
								= (temp3+temp2*2.0*temp1*exp(-temp1*temp1)
										/sqrt(pi))*const_columb/rijsq;
					}
					else if (isWolfOn)
					{
						uwolf_real = uwolf_real + charge[ii1]*charge[ii2]
								*(erfc(temp1)/rij + wolfvcon1 + wolfvcon2*(rij
										-rcutoffelec));
						fij = const_columb*temp2*(erfc(temp1)/rijsq + wolffcon1*exp(-temp1*temp1)/rij + wolffcon2);

					}
					fxij = fij*rxij;
					fyij = fij*ryij;
					fzij = fij*rzij;

					if (isWolfOn) // contribution to the virial
						virial_inter = virial_inter + fxij*rxij + fyij*ryij
								+ fzij*rzij;

					// force on atom ii1
					fxl[ii1] += fxij;
					fyl[ii1] += fyij;
					fzl[ii1] += fzij;
					// force on atom ii2
					fxl[ii2] -= fxij;
					fyl[ii2] -= fyij;
					fzl[ii2] -= fzij;
				} // rcutoffelecsq 
			} // is Ewald or Wolf Charge On


			// check if 1,3 distance < 1,3' distance
			// if it is true, no LJ between 1,3 needed
			// otherwise, calculate LJ interaction between 1 and 3'
			// use the old rxij, ryij and rzij for check
			if (fabs(rxij_old)>boxlx/2.0 || fabs(ryij_old)>boxly/2.0
					|| fabs(rzij_old)>boxlz/2.0)
			{
				fprintf(stderr,"Warning: long 1,3 angle ending pair %d-%d found...\n",ii1,ii2);
				fprintf(fpouts,
						"Warning: long 1,3 angle ending pair %d-%d found...\n",
						ii1, ii2);
				// minimum image convention already applied in ewald part if charge interaction is ON
				// however, charge interaction may not be on
				// So its better to calculate them again
				rxij = rxij - boxlx*rint(rxij/boxlx);
				ryij = ryij - boxly*rint(ryij/boxly);
				rzij = rzij - boxlz*rint(rzij/boxlz);
				rijsq = rxij*rxij + ryij*ryij + rzij*rzij;
				// LJ part, do not need check ghost atom since 1,3 are in the same molecule
				// here we check ghost atom because in this case, 1,3 are actually 1 and 3'
				// which are in two different molecule
				if (isghost[ii1]!=lj_ghost && isghost[ii2]!=lj_ghost && rijsq
						<rcutoffsq) // LJ cutoff
				{
					rij = sqrt(rijsq);
					if (isLJswitchOn) // check if switch for LJ is used
					{
						if (rijsq<rcutonsq)
							LJswitch = 1.0;
						else
							LJswitch = (rcutoffsq-rijsq)*(rcutoffsq-rijsq)
									*(rcutoffsq+2.0*rijsq-3.0*rcutonsq)
									/roff2_minus_ron2_cube;
					}
					sigmaij = 0.5*(sigma[ii1]+sigma[ii2]);
					epsilonij = sqrt(epsilon[ii1]*epsilon[ii2]);
					r_rijsq = sigmaij*sigmaij/rijsq;
					r_r6 = r_rijsq*r_rijsq*r_rijsq;
					r_r12 = r_r6*r_r6;
					r_r12_minus_r_r6 = r_r12 - r_r6;
					uij_vdw13img_temp = epsilonij*r_r12_minus_r_r6; // still need *4.0
					if (isLJswitchOn) // if switch is used
						uij_vdw13img += uij_vdw13img_temp*LJswitch; // still need 4.0
					else
						uij_vdw13img += uij_vdw13img_temp; // still need *4.0
					// calculate LJ forces
					if (isLJswitchOn)
					{
						if (rijsq<rcutonsq)
							fij = 24.0*epsilonij*(r_r12_minus_r_r6+r_r12)/rijsq;
						else
							fij = 24.0*epsilonij*(r_r12_minus_r_r6+r_r12)/rijsq
									*LJswitch // 4.0 for the real energy
									-4.0*uij_vdw13img_temp*12.0*(rcutoffsq
											-rijsq)*(rcutonsq-rijsq)
											/roff2_minus_ron2_cube;
					}
					else
						fij = 24.0*epsilonij*(r_r12_minus_r_r6+r_r12)/rijsq;
					fxij = fij*rxij;
					fyij = fij*ryij;
					fzij = fij*rzij;
					// contribution to the virial
					virial_inter = virial_inter + fxij*rxij + fyij*ryij + fzij
							*rzij;
					// force on atom ii1
					fxl[ii1] += fxij;
					fyl[ii1] += fyij;
					fzl[ii1] += fzij;
					// force on atom ii2
					fxl[ii2] -= fxij;
					fyl[ii2] -= fyij;
					fzl[ii2] -= fzij;
				} // rcutoffsq
			} // if 13 > 13'
		} // unique 1,3 pair check
	} // nangle loop

	// add into total energy
	uexcl += uij_excl_13; // still need constant
	uvdw += uij_vdw13img; // still need 4.0
	ureal += uij_real13; // still need constant
	uGz0 += uij_Gz0/(2.0*boxlz);
	return 0;
}

int loop_12()
{
	int ii;
	int ii1, ii2;
	double rxij, ryij, rzij;
	double rijsq, rij, r_rijsq;
	double r_r6, r_r12, r_r12_minus_r_r6;
	double sigmaij, epsilonij;
	double fij, fxij, fyij, fzij;
	double uij_vdw12img, uij_vdw12img_temp;
	double uij_real12;
	double uij_excl_12;
	double rxij_old, ryij_old, rzij_old;
	double temp1, temp2, temp3;
	double rhoklsq, kapparhokl_sq;
	double uij_Gz0; // 1D ewald
	gsl_function FF;

	FF.function = &deriv_inc_gamma;
	FF.params = 0;

	uij_excl_12 = 0.0;
	uij_vdw12img = 0.0;
	uij_real12 = 0.0;
	uij_Gz0 = 0.0;

	for (ii=0; ii<nbond; ii++)
	{
		ii1 = bond_idx[ii][0];
		ii2 = bond_idx[ii][1];
		rxij = xx[ii1] - xx[ii2];
		ryij = yy[ii1] - yy[ii2];
		rzij = zz[ii1] - zz[ii2];

		// Gz=0 term for 1D ewald
		if (isEwaldOn && fEwald_Dim==ewald_1D)
		{
			rhoklsq = rxij*rxij + ryij*ryij;
			if (rhoklsq>parallel_to_z_err)
			{
				kapparhokl_sq = kappasq*rhoklsq;
				temp1 = Euler_const + gsl_sf_gamma_inc(0.0,kapparhokl_sq) + log(kapparhokl_sq);
				uij_Gz0 = uij_Gz0 -charge[ii1]*charge[ii2]*temp1;
				// forces
				// derivative of incomplete gamma function
				// use numerical differential
				// This may not be a good way to do it
				// temp1 is the value, temp2 is the absolute std. err.
				gsl_deriv_central(&FF, kapparhokl_sq, 1.0e-8, &temp1, &temp2);
				temp2 = kappasq*temp1 + 1.0/rhoklsq;
				fij = charge[ii1]*charge[ii2]*temp2*const_columb/boxlz;
				fxij = fij*rxij;
				fyij = fij*ryij;
				// force on atom ii1
				fxl[ii1] += fxij;
				fyl[ii1] += fyij;
				// force on atom ii2
				fxl[ii2] -= fxij;
				fyl[ii2] -= fyij;
			} // if rhoklsq != 0.0 
		} // if 1D ewald is used

		// save the old positions
		rxij_old = rxij;
		ryij_old = ryij;
		rzij_old = rzij;

		if (isEwaldOn || isWolfOn) // if ewald is needed
		{
			// calculate the exclude ewald part, which uses the raw distance between atom
			// ii and jj, no minimum image convention
			rijsq = rxij*rxij + ryij*ryij + rzij*rzij;
			rij = sqrt(rijsq);
			temp2 = charge[ii1]*charge[ii2]/rij; // still need constant
			uij_excl_12 += temp2; // still need constant
			// forces
			fij = -const_columb*temp2/rijsq; // NOTE: the negative sign here, see comments for 1,3
			fxij = fij*rxij;
			fyij = fij*ryij;
			fzij = fij*rzij;
			if (isWolfOn) // negative contribution to the virial
				virial_inter = virial_inter + fxij*rxij + fyij*ryij + fzij*rzij;
			// forces on ii1
			fxl[ii1] += fxij;
			fyl[ii1] += fyij;
			fzl[ii1] += fzij;
			// forces on ii2
			fxl[ii2] -= fxij;
			fyl[ii2] -= fyij;
			fzl[ii2] -= fzij;

			// calculate the real part of ewald summation
			// minimum image convetion now
			rxij = rxij - boxlx*rint(rxij/boxlx);
			ryij = ryij - boxly*rint(ryij/boxly);
			rzij = rzij - boxlz*rint(rzij/boxlz);
			rijsq = rxij*rxij + ryij*ryij + rzij*rzij;
			if (rijsq<rcutoffelecsq) // elec cutoff
			{
				rij = sqrt(rijsq);
				temp1 = kappa*rij;
				temp2 = charge[ii1]*charge[ii2]/rij;
				temp3 = temp2*erfc(temp1);
				if (isEwaldOn)
				{
					uij_real12 += temp3; // still need constant
					// forces
					fij = (temp3+temp2*2.0*temp1*exp(-temp1*temp1)/sqrt(pi))
							*const_columb/rijsq;
				}
				else if (isWolfOn)
				{
					uwolf_real = uwolf_real + charge[ii1]*charge[ii2]
							*(erfc(temp1)/rij + wolfvcon1 + wolfvcon2*(rij
									-rcutoffelec));
					fij = const_columb*temp2*(erfc(temp1)/rijsq + wolffcon1*exp(-temp1*temp1)/rij + wolffcon2);
				}
				fxij = fij*rxij;
				fyij = fij*ryij;
				fzij = fij*rzij;

				if (isWolfOn) // contribution to the virial
					virial_inter = virial_inter + fxij*rxij + fyij*ryij + fzij
							*rzij;

				// force on atom ii1
				fxl[ii1] += fxij;
				fyl[ii1] += fyij;
				fzl[ii1] += fzij;
				// force on atom ii2
				fxl[ii2] -= fxij;
				fyl[ii2] -= fyij;
				fzl[ii2] -= fzij;
			} // rcutoffelecsq
		} // if ewald or wolf is needed


		// check if 1,2 distance < 1,2' distance
		// if it is true, no LJ needed
		// otherwise, calculate the LJ between 1 and 2'
		// use old rxij_old to make the check since the rxij could be modified in ewald
		// part if ewald is on
		if (fabs(rxij_old)>boxlx/2.0 || fabs(ryij_old)>boxly/2.0
				|| fabs(rzij_old)>boxlz/2.0)
		{
			fprintf(stderr,"Warning: long 1,2 bond ending pair %d-%d found...\n",ii1,ii2);
			fprintf(fpouts,
					"Warning: long 1,2 bond ending pair %d-%d found...\n", ii1,
					ii2);
			// minimum image convention, see more comments in 1,3 calculations
			rxij = rxij - boxlx*rint(rxij/boxlx);
			ryij = ryij - boxly*rint(ryij/boxly);
			rzij = rzij - boxlz*rint(rzij/boxlz);
			rijsq = rxij*rxij + ryij*ryij + rzij*rzij;
			// LJ part, also check ghost atoms
			// here we check ghost atom because in this case, 12 are actually 1 and 2'
			// which are in two different molecule
			if (isghost[ii1]!=lj_ghost && isghost[ii2]!=lj_ghost && rijsq
					<rcutoffsq) // LJ cutoff
			{
				rij = sqrt(rijsq);
				// if switch potential for LJ is on, then calculate the switch
				if (isLJswitchOn)
				{
					if (rijsq<rcutonsq)
						LJswitch = 1.0;
					else
						LJswitch = (rcutoffsq-rijsq)*(rcutoffsq-rijsq)
								*(rcutoffsq+2.0*rijsq-3.0*rcutonsq)
								/roff2_minus_ron2_cube;
				}
				sigmaij = 0.5*(sigma[ii1]+sigma[ii2]);
				epsilonij = sqrt(epsilon[ii1]*epsilon[ii2]);
				r_rijsq = sigmaij*sigmaij/rijsq;
				r_r6 = r_rijsq*r_rijsq*r_rijsq;
				r_r12 = r_r6*r_r6;
				r_r12_minus_r_r6 = r_r12 - r_r6;
				uij_vdw12img_temp = epsilonij*r_r12_minus_r_r6; // still need *4.0
				if (isLJswitchOn) // if switch is used
					uij_vdw12img += uij_vdw12img_temp*LJswitch; // still need 4.0
				else
					uij_vdw12img += uij_vdw12img_temp; // still need *4.0
				// calculate LJ forces
				if (isLJswitchOn)
				{
					if (rijsq<rcutonsq)
						fij = 24.0*epsilonij*(r_r12_minus_r_r6+r_r12)/rijsq;
					else
						fij = 24.0*epsilonij*(r_r12_minus_r_r6+r_r12)/rijsq
								*LJswitch // 4.0 for the real energy
								-4.0*uij_vdw12img_temp*12.0*(rcutoffsq-rijsq)
										*(rcutonsq-rijsq)/roff2_minus_ron2_cube;
				}
				else
					fij = 24.0*epsilonij*(r_r12_minus_r_r6+r_r12)/rijsq;
				fxij = fij*rxij;
				fyij = fij*ryij;
				fzij = fij*rzij;
				// contribution to the virial
				virial_inter = virial_inter + fxij*rxij + fyij*ryij + fzij*rzij;
				// force on atom ii1
				fxl[ii1] += fxij;
				fyl[ii1] += fyij;
				fzl[ii1] += fzij;
				// force on atom ii2
				fxl[ii2] -= fxij;
				fyl[ii2] -= fyij;
				fzl[ii2] -= fzij;
			} // rcutoffsq
		} // if 12 > 12'
	} // nbond loop

	// add into total energy
	uexcl += uij_excl_12; // still need constant
	uvdw += uij_vdw12img; // still need 4.0
	ureal += uij_real12; // still need constant
	uGz0 += uij_Gz0/(2.0*boxlz);
	return 0;
}

// calcualte interaction between nonbonded pairs
int loop_nbp()
{
	int ii;
	int ii1, ii2;
	double rxij, ryij, rzij;
	double rijsq, rij, r_rijsq;
	double r_r6, r_r12, r_r12_minus_r_r6;
	double epsilonij, sigmaij;
	double uij_vdwnbp, uij_vdwnbp_temp, uij_realnbp, uij_realnbp_temp;
	double uij_wolfnbp_temp;
	double fij, fxij, fyij, fzij;
	double temp1, temp2;
	double rhoklsq, kapparhokl_sq;
	double uij_Gz0nbp; // 1D ewald
	double rxij_old, ryij_old, rzij_old;
	double rijsq_old;
	int isSameMole;
	int fCalculate_nbpLJ, fCalculate_nbpelec;
	double uij_coulombnbp_temp;
	gsl_function FF;

	FF.function = &deriv_inc_gamma;
	FF.params = 0;

	uij_vdwnbp = 0.0;
	uij_realnbp = 0.0;
	uij_Gz0nbp = 0.0;

	// loop through all nonbonded pairs
	for (ii=0; ii<nnbp; ii++)
	{
		ii1 = nbp_idx[ii][0];
		ii2 = nbp_idx[ii][1];
		// save old rxij, ryij for possible 1D ewald
		rxij_old = rxij = xx[ii1] - xx[ii2];
		ryij_old = ryij = yy[ii1] - yy[ii2];
		rzij_old = rzij = zz[ii1] - zz[ii2];
		rijsq_old = rxij*rxij + ryij*ryij + rzij*rzij;
		// minimum image convention
		// if distance between ii1 and ii2 is smaller than
		// half box length, which means no change after
		// minimum image convention applied. It means that
		// the two atoms are indeed in the same molecule.
		// In this case, the interaction between them are 
		// always needed to be calculated even they are marked
		// as ghost, since ghost only means inter-molecular
		// interaction.
		// When the distance is larger than half box length,
		// the interaction is actually between ii1 and ii2' (image)
		// from two different molecules. Ghost check will be applied
		// in this case.
		rxij = rxij - boxlx*rint(rxij/boxlx);
		ryij = ryij - boxly*rint(ryij/boxly);
		rzij = rzij - boxlz*rint(rzij/boxlz);
		rijsq = rxij*rxij + ryij*ryij + rzij*rzij;
		// check if the position changed after the minimum image convention
		if (fabs(rijsq-rijsq_old)>tolerant_err)
		{
			fprintf(stderr,"Warning: long nonbonded pair %d-%d found...\n",ii1,ii2);
			fprintf(fpouts, "Warning: long nonbonded pair %d-%d found...\n",
					ii1, ii2);
			// image found, check ghost
			// only calculate if they are not ghost atoms
			isSameMole = false;
			if (isghost[ii1]==all_ghost || isghost[ii2]==all_ghost)
			{
				// if they are full ghost, dont calculate them 
				fCalculate_nbpLJ = false;
				fCalculate_nbpelec = false;
			}
			else if (isghost[ii1]==lj_ghost || isghost[ii2]==lj_ghost)
			{
				// if they are lj ghost, dont calculate LJ, but calculate electric
				fCalculate_nbpLJ = false;
				fCalculate_nbpelec = true;
			}
			else
			{
				// they are not ghost, calculate LJ and elec
				fCalculate_nbpLJ = true;
				fCalculate_nbpelec = true;
			}
		}
		else // they are in the same molecule, no ghost check need, always calculate everything
		{
			isSameMole = true;
			fCalculate_nbpLJ = true;
			fCalculate_nbpelec = true;

		} // end of same molecule check

		// calculate rij
		rij = sqrt(rijsq);

		// if LJ interactions are needed
		if (fCalculate_nbpLJ==true)
		{
			// cut off still needed to be checked
			if (rijsq<rcutoffsq)
			{
				// if switch potential for LJ is on, then calculate the switch
				if (isLJswitchOn)
				{
					if (rijsq<rcutonsq)
						LJswitch = 1.0;
					else
						LJswitch = (rcutoffsq-rijsq)*(rcutoffsq-rijsq)
								*(rcutoffsq+2.0*rijsq-3.0*rcutonsq)
								/roff2_minus_ron2_cube;
				}
				sigmaij = 0.5*(sigma[ii1]+sigma[ii2]);
				epsilonij = sqrt(epsilon[ii1]*epsilon[ii2]);
				r_rijsq = sigmaij*sigmaij/rijsq;
				r_r6 = r_rijsq*r_rijsq*r_rijsq;
				r_r12 = r_r6*r_r6;
				r_r12_minus_r_r6 = r_r12 - r_r6;
				uij_vdwnbp_temp = epsilonij*r_r12_minus_r_r6; // still need *4.0
				if (isLJswitchOn) // if switch is used
					uij_vdwnbp += uij_vdwnbp_temp*LJswitch; // still need 4.0 else 
				uij_vdwnbp += uij_vdwnbp_temp; // still need *4.0
				// force calculations
				if (isLJswitchOn)
				{
					if (rijsq<rcutonsq)
						fij = 24.0*epsilonij*(r_r12_minus_r_r6+r_r12)/rijsq;
					else
						fij = 24.0*epsilonij*(r_r12_minus_r_r6+r_r12)/rijsq
								*LJswitch // 4.0 for the real energy
								-4.0*uij_vdwnbp_temp*12.0*(rcutoffsq-rijsq)
										*(rcutonsq-rijsq)/roff2_minus_ron2_cube;
				}
				else
					fij = 24.0*epsilonij*(r_r12_minus_r_r6+r_r12)/rijsq;
				fxij = fij*rxij;
				fyij = fij*ryij;
				fzij = fij*rzij;
				// contribution to the virial
				virial_inter = virial_inter + fxij*rxij + fyij*ryij + fzij*rzij;
				// force on atom ii1
				fxl[ii1] += fxij;
				fyl[ii1] += fyij;
				fzl[ii1] += fzij;
				// force on atom ii2
				fxl[ii2] -= fxij;
				fyl[ii2] -= fyij;
				fzl[ii2] -= fzij;
			} // cut off check
		}// if LJ is needed

		// if electrostatic interactions are needed
		if (fCalculate_nbpelec == true)
		{
			// calculate Gz=0 term for 1D ewald summation
			// this has to be done before the minimum image convention
			// so rxij_old, ryij_old are used
			// not sure if cutoff check should be applied
			// If the cutoff is larger than the size in x,y dimension,
			// this should not be a problem
			if (isEwaldOn && fEwald_Dim==ewald_1D)
			{
				rhoklsq = rxij_old*rxij_old + ryij_old*ryij_old;
				if (rhoklsq>parallel_to_z_err)
				{
					kapparhokl_sq = kappasq*rhoklsq;
					temp1 = Euler_const + gsl_sf_gamma_inc(0.0,kapparhokl_sq) + log(kapparhokl_sq);
					uij_Gz0nbp = uij_Gz0nbp -charge[ii1]*charge[ii2]*temp1;
					// forces
					// derivative of incomplete gamma function
					// use numerical differential
					// This may not be a good way to do it
					// temp1 is the value, temp2 is the absolute std. err.
					gsl_deriv_central(&FF, kapparhokl_sq, 1.0e-8, &temp1,
							&temp2);
					temp2 = kappasq*temp1 + 1.0/rhoklsq;
					fij = charge[ii1]*charge[ii2]*temp2*const_columb/boxlz;
					fxij = fij*rxij_old;
					fyij = fij*ryij_old;
					// force on atom ii1
					fxl[ii1] += fxij;
					fyl[ii1] += fyij;
					// force on atom ii2
					fxl[ii2] -= fxij;
					fyl[ii2] -= fyij;
				} // if rhoklsq != 0.0 (>1.0e-10)
			} // if 1D ewald is used

			// electrostatic part
			if (isEwaldOn && rijsq<rcutoffelecsq)
			{
				uij_realnbp_temp = charge[ii1]*charge[ii2]/rij;
				uij_realnbp += uij_realnbp_temp*erfc(kappa*rij); // real part ewald energy, still need 1/4*pi*epsilon0
				temp1 = kappa*rij;
				temp2 = uij_realnbp_temp*erfc(temp1);
				// real part force calculation
				fij = (temp2 + uij_realnbp_temp*2.0*temp1*exp(-temp1*temp1)
						/sqrt(pi))*const_columb/rijsq;
				fxij = fij*rxij;
				fyij = fij*ryij;
				fzij = fij*rzij;
				// force on atom ii1
				fxl[ii1] += fxij;
				fyl[ii1] += fyij;
				fzl[ii1] += fzij;
				// force on atom ii2
				fxl[ii2] -= fxij;
				fyl[ii2] -= fyij;
				fzl[ii2] -= fzij;
			} // is Charge On check and rcutoffelecsq check

			// wolf method for electrostatic
			if (isWolfOn && rijsq<rcutoffelecsq)
			{
				uwolf_real = uwolf_real + charge[ii1]*charge[ii2]*(erfc(kappa
						*rij)/rij + wolfvcon1 + wolfvcon2*(rij-rcutoffelec));
				uij_wolfnbp_temp = charge[ii1]*charge[ii2]/rij;
				fij = const_columb*uij_wolfnbp_temp
				*(erfc(kappa*rij)/rijsq + wolffcon1*exp(-(kappa*rij)*(kappa*rij))/rij + wolffcon2);
				fxij = fij*rxij;
				fyij = fij*ryij;
				fzij = fij*rzij;
				// contribution to the virial
				virial_inter = virial_inter + fxij*rxij + fyij*ryij + fzij*rzij;
				// force on atom ii1
				fxl[ii1] += fxij;
				fyl[ii1] += fyij;
				fzl[ii1] += fzij;
				// force on atom ii2
				fxl[ii2] -= fxij;
				fyl[ii2] -= fyij;
				fzl[ii2] -= fzij;
			} // if ewald is needed

			// simple coulomb method
			if (isSimpleCoulomb && rijsq<rcutoffelecsq)
			{
				uij_coulombnbp_temp = charge[ii1]*charge[ii2]/rij; // need constant
				ucoulomb += uij_coulombnbp_temp; // need constant
				fij = const_columb*uij_coulombnbp_temp/rijsq;
				fxij = fij*rxij;
				fyij = fij*ryij;
				fzij = fij*rzij;
				// contribution to the virial
				virial_inter = virial_inter + fxij*rxij + fyij*ryij + fzij*rzij;
				// force on atom ii1
				fxl[ii1] += fxij;
				fyl[ii1] += fyij;
				fzl[ii1] += fzij;
				// force on atom ii2
				fxl[ii2] -= fxij;
				fyl[ii2] -= fyij;
				fzl[ii2] -= fzij;
			}

		} // if electrostatic is needed
	} // nonbonded pair loop

	// add into total energy
	uvdw += uij_vdwnbp; // still need 4.0
	unbp_vdw += uij_vdwnbp; // still need 4.0
	ureal += uij_realnbp; // still neeed constant
	uGz0 += uij_Gz0nbp/(2.0*boxlz); // still need constant
	return 0;
}

// fourier space sum and self interaction correction
int ewald_fourier_and_self()
{
	int ii;
	int kx, ky, kz, ksq;
	double rkx, rky, rkz, rksq;
	double kvec;
	double sr, si;
	double t;
	double fij;

	for (kx=-KMAXX; kx<=KMAXX; kx++) // NOTE: <=
	{
		for (ky=-KMAXY; ky<=KMAXY; ky++) // <= 
		{
			for (kz=-KMAXZ; kz<=KMAXZ; kz++) // <=
			{
				// check for 1D ewald summation
				// skip all Gz=0 terms
				if (fEwald_Dim==ewald_1D && kz==0)
					continue;
				ksq = kx*kx + ky*ky + kz*kz;
				if (ksq<=KSQMAX && ksq!=0)
				{
					rkx = TWOPI_LX*kx;
					rky = TWOPI_LY*ky;
					rkz = TWOPI_LZ*kz;
					rksq = rkx*rkx + rky*rky + rkz*rkz;
					kvec = exp(-Bfactor_ewald*rksq)/rksq;
					// calculate |rho(k)|^2 = sr*sr + si*si
					sr = 0.0;
					si = 0.0;
					// energy
					for (ii=0; ii<natom; ii++)
					{
						// the x,y,z coordinates dont have to be PBC'd
						// before doing the calculation since the existence of the TWOPI_L
						// factor will make the cos, sin functions have the same results
						// with periodical system
						t = rkx*xx[ii] + rky*yy[ii] + rkz*zz[ii];
						sr = sr + charge[ii]*cos(t);
						si = si + charge[ii]*sin(t);
					}
					ufourier += kvec*(sr*sr+si*si);
					// forces
					for (ii=0; ii<natom; ii++)
					{
						t = rkx*xx[ii] + rky*yy[ii] + rkz*zz[ii];
						fij = 2.0*charge[ii]*(sr*sin(t)-si*cos(t))
								*Vfactor_ewald*kvec*const_columb;
						fxl[ii] += fij*rkx;
						fyl[ii] += fij*rky;
						fzl[ii] += fij*rkz;
					}
				} // if (ksq<=KSQMAX && ksq!=0)
			} // KMAXZ
		} // KMAXY
	} // KMAXX


	// total fourier energy part of ewald
	ufourier = ufourier*Vfactor_ewald*const_columb;

	// self interaction corrections, constant, so no forces
	for (ii=0; ii<natom; ii++)
		uself += charge[ii]*charge[ii];
	uself = uself*const_columb*sqrt(kappa*kappa/pi);

	return 0;
}

int ewald_vacuum()
{
	int ii;
	double xxpri, yypri, zzpri; // PBC coords into primal box
	double qrx, qry, qrz;
	double chargei;
	double fij;

	// the atoms have to be grouped adjacent to their respective molecular
	// centers of mass before performing this sum
	// Because any molecules straddling a periodic boundary will cause
	// the total dipole moment of the cell to be greatly exaggerated.
	// -------------------------
	// calculate the center of mass and relative positions
	cal_com_and_inner_coords();
	// calcuate the new molecular center of mass using PBC
	for (ii=0; ii<nmole; ii++)
	{
		mole_xx[ii] = mole_xx[ii] - boxlx*rint(mole_xx[ii]/boxlx);
		mole_yy[ii] = mole_yy[ii] - boxly*rint(mole_yy[ii]/boxly);
		mole_zz[ii] = mole_zz[ii] - boxlz*rint(mole_zz[ii]/boxlz);
	}
	// calculate the new positions to the PBC'd center of mass
	reconstruct_coords_from_com();

	qrx = qry = qrz = 0.0;
	for (ii=0; ii<natom; ii++)
	{
		// use the reconstructed coordinates for this calculations
		// the reconstructed coordinates make sure that the
		// molecular center of mass are in the primal simulation
		// box and the atoms in one molecule are grouped together
		xxpri = ex[ii];
		yypri = fy[ii];
		zzpri = gz[ii];

		chargei = charge[ii];
		qrx = qrx + chargei*xxpri;
		qry = qry + chargei*yypri;
		qrz = qrz + chargei*zzpri;
	}
	// energy
	uvacuum = qrx*qrx + qry*qry + qrz*qrz; // still need constant
	uvacuum = uvacuum*twopi_over_3v*const_columb;

	// forces
	for (ii=0; ii<natom; ii++)
	{
		fij = -2.0*twopi_over_3v*const_columb*charge[ii];
		fxl[ii] += fij*qrx;
		fyl[ii] += fij*qry;
		fzl[ii] += fij*qrz;
	}

	// Does this term have additional contribution to the pressure?
	// Or is it already included in the ewald??

	return 0;
}

// calculate the wolf self interaction part
int wolf_con()
{
	int ii;
	for (ii=0; ii<natom; ii++)
	{
		uwolf_con = uwolf_con + charge[ii]*charge[ii];
	}
	uwolf_con = uwolf_con*const_columb*(kappa/sqrt(pi)+erfc(kappa*rcutoffelec)
			/(2.0*rcutoffelec));
	uwolf_real = uwolf_real*const_columb;
	uexcl *= const_columb;
	// total wolf
	uwolf = uwolf_real - uwolf_con - uexcl;

	return 0;
}

// calculate simple direct coulomb interactions between different molecules
// the cutoff is based on molecular center (center of mass)
int simple_coulomb_inter_mole()
{
	int ii, jj;
	int mm, nn;
	int first_atom_of_mole_mm, first_atom_of_mole_mm_plus_one;
	int first_atom_of_mole_nn, first_atom_of_mole_nn_plus_one;
	double uij_coulomb_temp;
	double rxmn, rymn, rzmn, rmnsq, rmn; // distance between molecular center of mass
	double rxij, ryij, rzij, rijsq, rij; // distance between atoms
	double fij, fxij, fyij, fzij;

	// first calculate all molecules center of mass
	cal_com_and_inner_coords();

	// calculate coulomb energy
	// for cross molecules
	for (mm=0; mm<nmole-1; mm++)
	{
		for (nn=mm+1; nn<nmole; nn++)
		{
			// calculate the distance between molecular center of mass
			rxmn = mole_xx[mm] - mole_xx[nn];
			rymn = mole_yy[mm] - mole_yy[nn];
			rzmn = mole_zz[mm] - mole_zz[nn];
			// minimum image convention
			rxmn = rxmn - boxlx*rint(rxmn/boxlx);
			rymn = rymn - boxly*rint(rymn/boxly);
			rzmn = rzmn - boxlz*rint(rzmn/boxlz);
			rmnsq = rxmn*rxmn + rymn*rymn + rzmn*rzmn;
			// the cut off check is based on the molecular center of mass
			// there could be different types of the check
			// e.g. based on distances of any sites between 2 molecules
			//      based on distances of any groups between 2 molecules
			// cutoff check
			if (rmnsq < rcutoffelecsq)
			{
				// get index for the atoms from each molecules
				first_atom_of_mole_mm = mole_first_atom_idx[mm];
				first_atom_of_mole_mm_plus_one = mole_first_atom_idx[mm+1];
				first_atom_of_mole_nn = mole_first_atom_idx[nn];
				first_atom_of_mole_nn_plus_one = mole_first_atom_idx[nn+1];
				// distance between two center of mass
				rmn = sqrt(rmnsq);
				// first atom
				for (ii=first_atom_of_mole_mm; ii
						<first_atom_of_mole_mm_plus_one; ii++)
				{
					// ghost check
					if (isghost[ii]==all_ghost)
						continue;
					// second atom
					for (jj=first_atom_of_mole_nn; jj
							<first_atom_of_mole_nn_plus_one; jj++)
					{
						// ghost check
						if (isghost[jj]==all_ghost)
							continue;
						// i         j
						//  \       /
						//   m-----n
						//  
						//  (i<-j) = (i<-m) + (m<-n) + (n<-j)
						//
						//  where x<-y = x - y, x and y are vectors
						// 
						// calculate the distance between atoms
						// based on the minimum image convention'd centers of mass
						rxij = ex[ii] + rxmn - ex[jj];
						ryij = fy[ii] + rymn - fy[jj];
						rzij = gz[ii] + rzmn - gz[jj];
						rijsq = rxij*rxij + ryij*ryij + rzij*rzij;
						rij = sqrt(rijsq);
						uij_coulomb_temp = charge[ii]*charge[jj]/rij; // need constant
						ucoulomb += uij_coulomb_temp; // need constant
						fij = const_columb*uij_coulomb_temp/rijsq;
						fxij = fij*rxij;
						fyij = fij*ryij;
						fzij = fij*rzij;
						// contribution to the virial
						virial_inter = virial_inter + fxij*rxij + fyij*ryij
								+ fzij*rzij;
						// force on atom ii
						fxl[ii] += fxij;
						fyl[ii] += fyij;
						fzl[ii] += fzij;
						// force on atom jj
						fxl[jj] -= fxij;
						fyl[jj] -= fyij;
						fzl[jj] -= fzij;
					} // second atom jj
				} // first atom ii
			} // cut off check
		} // inner layer molecules
	} // outer layer molecules

	// coulomb interactions from 1,4 and nonbonded atom pairs
	// are calculated via loop_14 and loop_nbp

	// constant for coulomb energy
	ucoulomb *= const_columb;

	return 0;
}

int erfrc()
{
	int ii;

	// zero energies
	uvdw = 0.0;
	unbp_vdw = 0.0;
	ureal = 0.0;
	uexcl = 0.0;
	ufourier = 0.0;
	uself = 0.0;
	uwolf = 0.0;
	uwolf_real = 0.0;
	uwolf_con = 0.0;
	usflj = 0.0;
	uvacuum = 0.0;
	uGz0 = 0.0;
	ucoulomb = 0.0;

	// pressure related
	virial_inter = 0.0;

	// zero forces
	for (ii=0; ii<natom; ii++)
		fxl[ii] = fyl[ii] = fzl[ii] = 0.0;

	// i-j loop, nonboned loop, 1-4 loop, 1-3 loop and 1-2 loop
	loop_ij();
	loop_nbp();
	loop_14();
	loop_13();
	loop_12();

	// 4.0 factor for LJ energies
	uvdw *= 4.0;
	unbp_vdw *= 4.0;

	if (isEwaldOn) // if ewald is on, calculate the fourier and self correction parts
	{
		// calculate fourier space sum for ewald and self interaction corrections
		ewald_fourier_and_self();
		// constant for ewald energies
		ureal *= const_columb;
		uexcl *= const_columb;
		// total 3D ewald energy with tinfoil boundary condition
		uewald = ureal + ufourier - uself - uexcl;

		// calculate vacuum boundary condition
		// energy/force term is needed
		if (fEwald_BC == ewald_bc_vacuum)
		{
			ewald_vacuum();
			uewald += uvacuum; // add into total ewald energy
		}
		// if 1D ewald is used
		if (fEwald_Dim==ewald_1D)
		{
			uGz0 *= const_columb;
			uewald += uGz0; // add Gz=0 term into total ewald energy for 1D ewald
		}
		// add the virial contribution from ewald
		virial_inter = virial_inter + uewald;
	}
	else if (isWolfOn) // if wolf is on
		wolf_con(); // calculate the self interaction term for wolf
	else if (isSimpleCoulomb)
		simple_coulomb_inter_mole(); // if simple coulomb is used, calculate inter-mole coulomb

	// calculate solid fluid energy if necessary
	if (isSFon)
		sffrc();

	// total inter molecule energy add up everything
	// they should be zero if they are not used
	uinter = uvdw + uewald + uwolf + ucoulomb + usflj;

	return 0;
}

