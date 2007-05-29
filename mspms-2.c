/*
 * Maintainable Simplex Purpose Molecular Simulator 2
 * Rewritten with standard C language
 * The goal is simpler, quicker, easier, better.
 * No class and other complex data structures.
 * Use common names for input files. So the every job must
 * have its own working directory.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "vars.h"
#include "random.h"


/* Initiate variables */
int init_vars()
{
    int ii;
    /* change atom weight unit from g/mol to kg/mol for future calculations */
    for (ii=0;ii<natom;ii++)
	aw[ii] *= 0.001;

    /* long range corrections */

    // initialize random number generator
    rmarin(ij, jk);

    // calculate the degree of freedom
    nfree = 3*natom - nconstraint;
}

/* Read in input and config files */
int readins()
{
    int ii;
    const int datalen = 200;
    char buffer[200];

    /* read input file */
    fpins = fopen(INPUT,"r");
    
    sscanf(fgets(buffer,datalen,fpins), "%d %d", &ij, &jk);
    sscanf(fgets(buffer,datalen,fpins), "%f", &treq);
    sscanf(fgets(buffer,datalen,fpins), "%d", &nconstraint);
    

    fclose(fpins);

    /* read config file */
    fpcfg = fopen(CONFIG,"r");
    // read the first line of comments
    // fscanf(fpcfg, "%[^\n]", buffer); /* read in everything to the buffer except the return */
    fgets(sysname,datalen,fpcfg);

    // read in atom list
    fscanf(fpcfg, "%d atoms\n", &natom);
    assert(natom<=natom_max);
    for (ii=0;ii<natom;ii++)
	fscanf(fpcfg,"%f %f %f %f %d %d\n",&aw[ii],&epsilon[ii],&sigma[ii],
		&charge[ii],&isghost[ii],&tasostype[ii]);

    // read in bond list
    fscanf(fpcfg, "%d bonds\n", &nbond);
    assert(nbond<=nbond_max);
    for (ii=0;ii<nbond;ii++)
	fscanf(fpcfg,"%d %d %d %f %f\n",&bond_idx[ii][0],&bond_idx[ii][1],
		&bond_type[ii],&Kb[ii],&Req[ii]);

    // read in angle list
    fscanf(fpcfg, "%d angles\n", &nangle);
    assert(nangle<=nangle_max);
    for (ii=0;ii<nangle;ii++)
    {
	fscanf(fpcfg,"%d %d %d %d %f %f\n",&angle_idx[ii][0],&angle_idx[ii][1],
		&angle_idx[ii][2],&angle_type[ii],&Ktheta[ii],&Thetaeq[ii]);
    }

    // read in dihedral list
    fscanf(fpcfg, "%d dihedrals\n", &ndih);
    assert(ndih<=ndih_max);
    for (ii=0;ii<ndih;ii++)
    {
	fscanf(fpcfg,"%d %d %d %d %d %f %f %f %f\n",&dih_idx[ii][0],&dih_idx[ii][1],
		&dih_idx[ii][2],&dih_idx[ii][3],&dih_type[ii],&c1[ii],&c2[ii],&c3[ii],&c4[ii]);
    }

    // read in improper list
    fscanf(fpcfg, "%d impropers\n", &nimp);
    assert(nimp<=nimp_max);

    // read in nonbonded pair list
    fscanf(fpcfg, "%d nonbonded\n", &nnbp);
    assert(nnbp<=nnbp_max);


    fclose(fpcfg);
}

int echo()
{
    fprintf(stderr,"natom=%d\n",natom);
    fprintf(stderr,"nconstraint=%d\n",nconstraint);

    fprintf(stderr,"%s",sysname);
    fprintf(stderr,"%d atoms.\n",natom);
    fprintf(stderr,"%d bonds.\n",nbond);
    fprintf(stderr,"%d angles.\n",nangle);
    fprintf(stderr,"%d dihedrals.\n",ndih);
    fprintf(stderr,"%d impropers.\n",nimp);
    fprintf(stderr,"%d nonbonded pairs.\n",nnbp);
}

int make_exclude_list()
{}

int velinit()
{
    int ii, jj, kk;
    float px, py, pz;
    float stdvtmp, stdv;
    float totalmass;
    float scaling;

    totalmass = 0.0;
    px = py = pz = 0.0;
    // mv^2=kT, when m is kg/mol, the equation becomes to mv^2=RT
    // p = sqrt(RTm), v = sqrt(RT/m)
    stdvtmp = sqrt(Rgas*treq);

    // Temperature is K. R is J/K/mol. m is kg/mol. v is m/s
    for (ii=0;ii<natom;ii++)
    {
	stdv = stdvtmp/sqrt(aw[ii]);
	vx[ii] = stdv*gaussran();
	vy[ii] = stdv*gaussran();
	vz[ii] = stdv*gaussran();
	px += vx[ii]*aw[ii];
	py += vy[ii]*aw[ii];
	pz += vz[ii]*aw[ii];
	totalmass += aw[ii];
    }
    // zero the momentum
    px /= totalmass;
    py /= totalmass;
    pz /= totalmass;
    for (ii=0;ii<natom;ii++)
    {
	vx[ii] -= px;
	vy[ii] -= py;
	vz[ii] -= pz;
    }
    // rescale velocity for required temperature
    ukin = 0.0;
    for (ii=0;ii<natom;ii++)
	ukin += aw[ii]*(vx[ii]*vx[ii]+vy[ii]*vy[ii]+vz[ii]*vz[ii]);
    ukin = 0.5*ukin;
    tinst = 2.0*ukin/(Rgas*nfree);
    scaling = sqrt(treq/tinst);
    for (ii=0;ii<natom;ii++)
    {
	vx[ii] *= scaling;
	vy[ii] *= scaling;
	vz[ii] *= scaling;
    }
    // recalculate the kinetic energy and instantaneous temperature
    // should be exactly the set tempature
    ukin = 0.0;
    for (ii=0;ii<natom;ii++)
	ukin += aw[ii]*(vx[ii]*vx[ii]+vy[ii]*vy[ii]+vz[ii]*vz[ii]);
    ukin = 0.5*ukin;
    tinst = 2.0*ukin/(Rgas*nfree);
}

int md()
{
    /*
       velinit();
     */
}

int main (int argc, char *argv[])
{
    readins();
    init_vars();
    echo();

    velinit();

}

