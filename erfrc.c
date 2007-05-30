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

int erfrc()
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
    float uij_vdw_temp, uij_elec_temp;

    // zero forces
    for (ii=0;ii<natom;ii++)
	fxl[ii] = fyl[ii] = fzl[ii] = 0.0;

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

		if (rijsq<rcutoffsq)
		{
		    rij = sqrt(rijsq);
		}
		sigmaij = 0.5*(sigmai+sigma[jj]);
		epsilonij = sqrt(epsiloni*epsilon[jj]);
		r_rijsq = sigmaij*sigmaij/rijsq;
		r_r6 = r_rijsq*r_rijsq*r_rijsq;
		r_r12 = r_r6*r_r6;
		r_r12_minus_r_r6 = r_r12 - r_r6;
		// uij_vdw_temp = epsilonij*(r_r12_minus_r_r6-r_r6);

		uij_elec_temp = chargei*charge[jj]/rij;

		
	    }
	}
    }


}
