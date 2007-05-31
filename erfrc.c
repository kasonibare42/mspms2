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
    float uij_vdw, uij_vdw_temp, uij_elec_temp;
    float fij, fxij, fyij, fzij;
    float temp1, temp2;

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
	    if (isNotexcl) // calculate interactions if it is not exclusion pair
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

		// LJ part
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
		if (rijsq<rcutoffelecsq)
		{
		    uij_elec_temp = chargei*charge[jj]/rij;
		    ureal = ureal + uij_elec_temp*erfc(kappa*rij); // real part energy, still need 1/4*pi*epsilon0
		    temp1 = kappa*rij;
		    temp2 = uij_elec_temp*erfc(temp1);
		    // real part force calculation
		    fij = (temp2 + uij_elec_temp*2.0*temp1*exp(-temp1*temp1)/sqrt(pi))*const_columb/rijsq;
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
		}


	    } // end of excluding list
	} // end of atom jj
	fxl[ii] = fxi;
	fyl[ii] = fyi;
	fzl[ii] = fzi;
    } // end of atom ii

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
    float uij_vdw14, uij_vdw14_temp, uij_elec14_temp;
    float fij, fxij, fyij, fzij;

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
	    }






	}
    }
}

int erfrc()
{
    int ii;

    // zero energies
    uvdw = 0.0;
    ureal = 0.0;

    // zero forces
    for (ii=0;ii<natom;ii++)
	fxl[ii] = fyl[ii] = fzl[ii] = 0.0;

    // make different functions for ii,jj loop
    // 1-4 loop, 1-3 loop and 1-2 loop

    loop_ij();
    loop_14();
}
