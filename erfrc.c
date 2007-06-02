/* Calculate inter atom and inter molecule energies and forces.
   These include different atoms on the same molecule (exclude 1-2,1-3,1-4)
   1-4 will be still be calculated after
   and atoms on different molecules.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "vars.h"

// calcualte interaction between atom ii and jj, exclude those from excluding list
int loop_ij()
{
    int ii, jj, kk;
    float xxi, yyi, zzi;
    float fxi, fyi, fzi;
    float rxij, ryij, rzij;
    float rijsq, rij, r_rijsq;
    float r_r6, r_r12, r_r12_minus_r_r6;
    float epsiloni, sigmai, chargei;
    float epsilonij, sigmaij;
    int isNotexcl[natom_max];
    float uij_vdw, uij_vdw_temp, uij_real, uij_real_temp;
    float fij, fxij, fyij, fzij;
    float temp1, temp2;

    uij_vdw = 0.0;
    uij_real = 0.0;

    // atom ii
    for (ii=0;ii<natom-1;ii++)
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
	for (kk=0;kk<natom;kk++)
	    isNotexcl[kk] = true;
	for (kk=pointexcl[ii];kk<pointexcl[ii+1];kk++)
	    isNotexcl[excllist[kk]] = false;
	// atom jj
	for (jj=ii+1;jj<natom;jj++)
	{
	    if (isghost[jj] == all_ghost) // move to jj+1 if atom jj is full ghost atom
		continue;
	    if (isNotexcl[jj]) // calculate interactions if it is not exclusion pair
	    {
		rxij = xxi - xx[jj];
		ryij = yyi - yy[jj];
		rzij = zzi - zz[jj];
		// minimum image convention
		rxij = rxij - boxlx*rint(rxij/boxlx);
		ryij = ryij - boxly*rint(ryij/boxly);
		rzij = rzij - boxlz*rint(rzij/boxlz);
		rijsq = rxij*rxij + ryij*ryij + rzij*rzij;
		rij = sqrt(rijsq);

		// LJ part, still need ghost atom check
		if (rijsq<rcutoffsq)
		{
		    sigmaij = 0.5*(sigmai+sigma[jj]);
		    epsilonij = sqrt(epsiloni*epsilon[jj]);
		    r_rijsq = sigmaij*sigmaij/rijsq;
		    r_r6 = r_rijsq*r_rijsq*r_rijsq;
		    r_r12 = r_r6*r_r6;
		    r_r12_minus_r_r6 = r_r12 - r_r6;
		    uij_vdw_temp = epsilonij*r_r12_minus_r_r6; // still need *4.0
		    uij_vdw += uij_vdw_temp; // still need *4.0
		    // force calculations
		    fij = epsilonij*(r_r12_minus_r_r6+r_r12)/rijsq;
		    fxij = 24.0*fij*rxij;
		    fyij = 24.0*fij*ryij;
		    fzij = 24.0*fij*rzij;
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
		    fij = (temp2 + uij_real_temp*2.0*temp1*exp(-temp1*temp1)/sqrt(pi))*const_columb/rijsq;
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


	    } // end of excluding list
	} // end of atom jj
	fxl[ii] = fxi;
	fyl[ii] = fyi;
	fzl[ii] = fzi;
    } // end of atom ii

    uvdw += uij_vdw; // still need 4.0
    ureal += uij_real; // still neeed constant
}

// calculate the interactions between dihedral ending atom 1 and 4
int loop_14()
{
    int ii;
    int ii1, ii2;
    float rxij, ryij, rzij;
    float rijsq, rij, r_rijsq;
    float r_r6, r_r12, r_r12_minus_r_r6;
    float sigmaij, epsilonij;
    float uij_vdw14, uij_vdw14_temp, uij_real14, uij_real14_temp;
    float fij, fxij, fyij, fzij;
    float temp1, temp2;

    uij_vdw14 = 0.0;
    uij_real14 = 0.0;

    for (ii=0;ii<ndih;ii++)
    {
	if (isDih_unique[ii]) // only calculate interactions for unique 1,4 pairs
	{
	    ii1 = dih_idx[ii][0];
	    ii2 = dih_idx[ii][3];
	    rxij = xx[ii1] - xx[ii2];
	    ryij = yy[ii1] - yy[ii2];
	    rzij = zz[ii1] - zz[ii2];
	    // check for necessary 1,4 parameter modifiers
	    // and calculate the corresponding sigmaij and epsilonij;
	    // OPLS modifier
	    sigmaij = 0.5*(sigma[ii1]+sigma[ii2]);
	    epsilonij = f0*sqrt(epsilon[ii1]*epsilon[ii2]);
	    // can add CHARMM modifier here, different 1,4 parameters
	    // must be done before minimum image convention

	    // minimum image convention
	    rxij = rxij - boxlx*rint(rxij/boxlx);
	    ryij = ryij - boxly*rint(ryij/boxly);
	    rzij = rzij - boxlz*rint(rzij/boxlz);
	    rijsq = rxij*rxij + ryij*ryij + rzij*rzij;
	    rij = sqrt(rijsq);

	    // LJ part
	    if (rijsq<rcutoffsq)
	    {
		r_rijsq = sigmaij*sigmaij/rijsq;
		r_r6 = r_rijsq*r_rijsq*r_rijsq;
		r_r12 = r_r6*r_r6;
		r_r12_minus_r_r6 = r_r12 - r_r6;
		uij_vdw14_temp = epsilonij*r_r12_minus_r_r6; // still need *4.0
		uij_vdw14 += uij_vdw14_temp; // still need *4.0
		// force calculations
		fij = epsilonij*(r_r12_minus_r_r6+r_r12)/rijsq;
		fxij = 24.0*fij*rxij;
		fyij = 24.0*fij*ryij;
		fzij = 24.0*fij*rzij;
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
	    if (isEwaldOn && rijsq<rcutoffelecsq)
	    {
		uij_real14_temp = charge[ii1]*charge[ii2]/rij;
		uij_real14 += uij_real14_temp*erfc(kappa*rij); // real part ewald energy, still need 1/4*pi*epsilon0
		temp1 = kappa*rij;
		temp2 = uij_real14_temp*erfc(temp1);
		// real part force calculation
		fij = (temp2 + uij_real14_temp*2.0*temp1*exp(-temp1*temp1)/sqrt(pi))*const_columb/rijsq;
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
	} // unique 1,4 pair check
    } // loop through dihedrals

    // add into total energy
    uvdw += uij_vdw14; // still need 4.0
    ureal += uij_real14; // still need constant
}

int loop_13()
{
    int ii;
    int ii1, ii2;
    float rxij, ryij, rzij;
    float rijsq, rij, r_rijsq;
    float r_r6, r_r12, r_r12_minus_r_r6;
    float sigmaij, epsilonij;
    float fij, fxij, fyij, fzij;
    float uij_vdw13img, uij_vdw13img_temp;
    float uij_real13, uij_real13_temp;
    float uij_excl_13;
    float rxij_old, ryij_old, rzij_old;
    float temp1, temp2, temp3;

    uij_excl_13 = 0.0;
    uij_vdw13img = 0.0;
    uij_real13 = 0.0;

    for (ii=0;ii<nangle;ii++)
    {
	if (isAngle_unique) // only calculate interactions for unique 1,3 pairs
	{
	    ii1 = angle_idx[ii][0];
	    ii2 = angle_idx[ii][2];
	    rxij = xx[ii1] - xx[ii2];
	    ryij = yy[ii1] - yy[ii2];
	    rzij = zz[ii1] - zz[ii2];
	    // save the old positions
	    rxij_old = rxij;
	    ryij_old = ryij;
	    rzij_old = rzij;

	    if (isEwaldOn) // if we need charge interactions
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
		fxl[ii1] += fxij;
		fyl[ii1] += fyij;
		fzl[ii1] += fzij;
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
		    uij_real13 += temp3; // still need constant
		    // forces
		    fij = (temp3+temp2*2.0*temp1*exp(-temp1*temp1)/sqrt(pi))*const_columb/rijsq;
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
		} // rcutoffelecsq 
	    } // is Charge On

	    // check if 1,3 distance < 1,3' distance
	    // if it is true, no LJ between 1,3 needed
	    // otherwise, calculate LJ interaction between 1 and 3'
	    // use the old rxij, ryij and rzij for check
	    if (fabs(rxij_old)>boxlx/2.0 || fabs(ryij_old)>boxly/2.0 || fabs(rzij_old)>boxlz/2.0)
	    {
		// minimum image convention already applied in ewald part if charge interaction is ON
		// however, charge interaction may not be on
		// So its better to calculate them again
		rxij = rxij - boxlx*rint(rxij/boxlx);
		ryij = ryij - boxly*rint(ryij/boxly);
		rzij = rzij - boxlz*rint(rzij/boxlz);
		rijsq = rxij*rxij + ryij*ryij + rzij*rzij;
		if (rijsq<rcutoffsq) // LJ cutoff
		{
		    rij = sqrt(rijsq);
		    sigmaij = 0.5*(sigma[ii1]+sigma[ii2]);
		    epsilonij = sqrt(epsilon[ii1]*epsilon[ii2]);
		    r_rijsq = sigmaij*sigmaij/rijsq;
		    r_r6 = r_rijsq*r_rijsq*r_rijsq;
		    r_r12 = r_r6*r_r6;
		    r_r12_minus_r_r6 = r_r12 - r_r6;
		    uij_vdw13img_temp = epsilonij*r_r12_minus_r_r6; // still need *4.0
		    uij_vdw13img += uij_vdw13img_temp; // still need *4.0
		    // calculate LJ forces
		    fij = epsilonij*(r_r12_minus_r_r6+r_r12)/rijsq;
		    fxij = 24.0*fij*rxij;
		    fyij = 24.0*fij*ryij;
		    fzij = 24.0*fij*rzij;
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
}

int loop_12()
{
    int ii;
    int ii1, ii2;
    float rxij, ryij, rzij;
    float rijsq, rij, r_rijsq;
    float r_r6, r_r12, r_r12_minus_r_r6;
    float sigmaij, epsilonij;
    float fij, fxij, fyij, fzij;
    float uij_vdw12img, uij_vdw12img_temp;
    float uij_real12, uij_real12_temp;
    float uij_excl_12;
    float rxij_old, ryij_old, rzij_old;
    float temp1, temp2, temp3;

    uij_excl_12 = 0.0;
    uij_vdw12img = 0.0;
    uij_real12 = 0.0;

    for (ii=0;ii<nbond;ii++)
    {
	ii1 = bond_idx[ii][0];
	ii2 = bond_idx[ii][1];
	rxij = xx[ii1] - xx[ii2];
	ryij = yy[ii1] - yy[ii2];
	rzij = zz[ii1] - zz[ii2];
	// save the old positions
	rxij_old = rxij;
	ryij_old = ryij;
	rzij_old = rzij;

	if (isEwaldOn) // if ewald is needed
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
	    fxl[ii1] += fxij;
	    fyl[ii1] += fyij;
	    fzl[ii1] += fzij;
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
		uij_real12 += temp3; // still need constant
		// forces
		fij = (temp3+temp2*2.0*temp1*exp(-temp1*temp1)/sqrt(pi))*const_columb/rijsq;
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
	    } // rcutoffelecsq
	} // if ewald is needed

	// check if 1,2 distance < 1,2' distance
	// if it is true, no LJ needed
	// otherwise, calculate the LJ between 1 and 2'
	// use old rxij_old to make the check since the rxij could be modified in ewald
	// part if ewald is on
	if (fabs(rxij_old)>boxlx/2.0 || fabs(ryij_old)>boxly/2.0 || fabs(rzij_old)>boxlz/2.0)
	{
	    // minimum image convention, see more comments in 1,3 calculations
	    rxij = rxij - boxlx*rint(rxij/boxlx);
	    ryij = ryij - boxly*rint(ryij/boxly);
	    rzij = rzij - boxlz*rint(rzij/boxlz);
	    rijsq = rxij*rxij + ryij*ryij + rzij*rzij;
	    if (rijsq<rcutoffsq) // LJ cutoff
	    {
		rij = sqrt(rijsq);
		sigmaij = 0.5*(sigma[ii1]+sigma[ii2]);
		epsilonij = sqrt(epsilon[ii1]*epsilon[ii2]);
		r_rijsq = sigmaij*sigmaij/rijsq;
		r_r6 = r_rijsq*r_rijsq*r_rijsq;
		r_r12 = r_r6*r_r6;
		r_r12_minus_r_r6 = r_r12 - r_r6;
		uij_vdw12img_temp = epsilonij*r_r12_minus_r_r6; // still need *4.0
		uij_vdw12img += uij_vdw12img_temp; // still need *4.0
		// calculate LJ forces
		fij = epsilonij*(r_r12_minus_r_r6+r_r12)/rijsq;
		fxij = 24.0*fij*rxij;
		fyij = 24.0*fij*ryij;
		fzij = 24.0*fij*rzij;
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
    } // nbond loop

    // add into total energy
    uexcl += uij_excl_12; // still need constant
    uvdw += uij_vdw12img; // still need 4.0
    ureal += uij_real12; // still need constant
}

int ewald_fourier_and_self()
{
    int ii;
    int kx, ky, kz, ksq;
    float rkx, rky, rkz, rksq;
    float kvec;
    float sr, si;
    float t;
    float fij;

    for (kx=-KMAXX;kx<=KMAXX;kx++) // NOTE: <=
    {
	for (ky=-KMAXY;ky<=KMAXY;ky++) // <= 
	{
	    for (kz=-KMAXZ;kz<=KMAXZ;kz++) // <=
	    {
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
		    for (ii=0;ii<natom;ii++)
		    {
		       	t = rkx*xx[ii] + rky*yy[ii] + rkz*zz[ii];   
			sr = sr + charge[ii]*cos(t);
		       	si = si + charge[ii]*sin(t);
		    }
		    ufourier += kvec*(sr*sr+si*si);
		    // forces
		    for (ii=0;ii<natom;ii++)
		    {
			t = rkx*xx[ii] + rky*yy[ii] + rkz*zz[ii];
			fij = 2.0*charge[ii]*(sr*sin(t)-si*cos(t))*Vfactor_ewald*kvec*const_columb;
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
    for (ii=0;ii<natom;ii++)
	uself += charge[ii]*charge[ii];
    uself = uself*const_columb*sqrt(kappa*kappa/pi);
}

int erfrc()
{
    int ii;

    // zero energies
    uvdw = 0.0;
    ureal = 0.0;
    uexcl = 0.0;
    ufourier = 0.0;
    uself = 0.0;

    // zero forces
    for (ii=0;ii<natom;ii++)
	fxl[ii] = fyl[ii] = fzl[ii] = 0.0;

    // i-j loop, 1-4 loop, 1-3 loop and 1-2 loop
    loop_ij();
    loop_14();
    loop_13();
    loop_12();
    // add 4.0 and constant to LJ and real part ewald energy
    uvdw *= 4.0;
    ureal *= const_columb;
    uexcl *= const_columb;

    if (isEwaldOn) // if ewald is on, calculate the fourier and self correction parts
	ewald_fourier_and_self();

    // total ewald energy
    uewald = ureal + ufourier - uself - uexcl;

    // total inter molecule energy
    uinter = uvdw + uewald;
}


