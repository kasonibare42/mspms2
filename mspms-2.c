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

extern int erfrc();
extern int rafrc();
extern int vver_nh();
extern int vver_nh_1();
extern int vver_nh_2();
extern int init_tasos_grid();

/* Initiate variables */
int init_vars()
{
    int ii, jj;
    FILE *fpcoords;
    char buffer[200];

    // initiate files
    fplog = fopen(LOGFILE,"w");
    fptrj = fopen(MOVIE,"wb");
    fpouts = fopen(OUTPUT,"w");

    nframe = 0; // number of frames in trajectory file

    for (ii=0;ii<num_counter_max;ii++)
    {
	icounter[ii] = 0;
	for (jj=0;jj<5;jj++)
	    accumulator[ii][jj] = 0.0;
    }

    /* change atom weight unit from g/mol to kg/mol for future calculations */
    for (ii=0;ii<natom;ii++)
	aw[ii] *= 0.001;

    /* long range corrections */

    // initialize random number generator
    rmarin(ij, jk);

    // calculate the degree of freedom
    nfree = 3*natom - nconstraint;

    // cutoff related
    rcutoffsq = rcutoff*rcutoff;
    rcutoffelecsq = rcutoffelec*rcutoffelec;
    rcutonsq = rcuton*rcuton;
    roff2_minus_ron2_cube = (rcutoffsq-rcutonsq)*(rcutoffsq-rcutonsq)*(rcutoffsq-rcutonsq);

    // delt
    deltby2 = delt/2.0;
    delts = delt/nstep_inner;
    deltsby2 = delts/2.0;

    // nose hoover
    // following for Dr. Maginn's nose hoover
    // they are not used anymore
    // dt_outer2 = deltby2;
    // dt_outer4 = dt_outer2/2.0;
    // ukin_nhts = 0.0;
    // upot_nhts = 0.0;
    // NRT = Rgas*treq*nfree;
    // Gts = 0.0;
    // vts = 0.0;
    // rts = 0.0;
    /*
     * Qts = Rgas*treq*nfree/Omega
     * where Omega is a parameter related to the mass of the thermostat
     * for this program, we read in the Qts directly.
     */
    delt_sqby2 = delt*delt/2.0;
    delts_sqby2 = delts*delts/2.0;
    unhts = 0.0;
    gg = nfree; // need double check
    ss = 0.0;
    ps = 0.0;
    ggs = nfree;
    sss = 0.0;
    pss = 0.0;
    unhtss = 0.0;


    // sf energy, tasos initiate part
    if (isSFon && sf_type==nanotube_tasos)
	init_tasos_grid();

    // check unique for dihedrals
    // this is for possible ring structures where 1,4 atoms can form multiple dihedrals
    // e.g. 1-2-3-4
    //       \5-6/
    // 1234 and 1564
    // The 14 pair should only be calculated once for energy/force
    // Thats the unique check for
    for (ii=0;ii<ndih;ii++)
	isDih_unique[ii] = true;
    for (ii=0;ii<ndih-1;ii++)
    {
	if (isDih_unique[ii])
	{
	    for (jj=ii;jj<ndih;jj++)
	    {
		if (isDih_unique[jj])
		{
		    if (dih_idx[ii][0]==dih_idx[jj][0] && dih_idx[ii][3]==dih_idx[jj][3])
			isDih_unique[jj] = false;
		    else if (dih_idx[ii][0]==dih_idx[jj][3] && dih_idx[ii][3]==dih_idx[jj][0])
			isDih_unique[jj] = false;
		}
	    }
	}
    }

    // following codes make sure 14 and 13 do not share same ending pairs
    // this is also for ring kind structures
    for (ii=0;ii<ndih;ii++)
    {
	if (isDih_unique[ii])
	{
	    for (jj=0;jj<nangle;jj++)
	    {
		if (dih_idx[ii][0]==angle_idx[jj][0] && dih_idx[ii][3]==angle_idx[jj][2])
		    isDih_unique[jj] = false;
		else if (dih_idx[ii][0]==angle_idx[jj][2] && dih_idx[ii][3]==angle_idx[jj][0])
		    isDih_unique[jj] = false;
	    }
	}
    }

    // following codes make sure 14 and 12 do not share the same ending pairs
    for (ii=0;ii<ndih;ii++)
    {
	if (isDih_unique[ii])
	{
	    for (jj=0;jj<nbond;jj++)
	    {
		if (dih_idx[ii][0]==bond_idx[jj][0] && dih_idx[ii][3]==bond_idx[jj][1])
		    isDih_unique[jj] = false;
		else if (dih_idx[ii][0]==bond_idx[jj][1] && dih_idx[ii][3]==bond_idx[jj][0])
		    isDih_unique[jj] = false;
	    }
	}
    }

    // check unique for angles
    // see above comments for dihedrals
    for (ii=0;ii<nangle;ii++)
	isAngle_unique[ii] = true;
    for (ii=0;ii<nangle-1;ii++)
    {
	if (isAngle_unique[ii])
	{
	    for (jj=ii;jj<nangle;jj++)
	    {
		if (isAngle_unique[jj])
		{
		    if (angle_idx[ii][0]==angle_idx[jj][0] && angle_idx[ii][2]==angle_idx[jj][2])
			isAngle_unique[jj] = false;
		    else if (angle_idx[ii][0]==angle_idx[jj][2] && angle_idx[ii][2]==angle_idx[jj][0])
			isAngle_unique[jj] = false;
		}
	    }
	}
    } // end of checking unique angles

    // following codes make sure 13 and 12 do not share the same ending pairs
    for (ii=0;ii<nangle;ii++)
    {
	if (isAngle_unique[ii])
	{
	    for (jj=0;jj<nbond;jj++)
	    {
		if (angle_idx[ii][0]==bond_idx[jj][0] && angle_idx[ii][2]==bond_idx[jj][1])
		    isAngle_unique[jj] = false;
		else if (angle_idx[ii][0]==bond_idx[jj][1] && angle_idx[ii][2]==bond_idx[jj][0])
		    isAngle_unique[jj] = false;
	    }
	}
    }

    // set ewald parameters
    Bfactor_ewald = 1.0/(4.0*kappa*kappa);
    Vfactor_ewald = 2.0*pi/(boxlx*boxly*boxlz);
    TWOPI_LX = 2.0*pi/boxlx;
    TWOPI_LY = 2.0*pi/boxly;
    TWOPI_LZ = 2.0*pi/boxlz;

    // read in coordinates
    fpcoords = fopen(coords_file,"r");
    fscanf(fpcoords, "%[^\n]", buffer);
    fscanf(fpcoords, "%[^\n]", buffer);
    for (ii=0;ii<natom;ii++)
    {
	fscanf(fpcoords,"%s %lf %lf %lf\n",buffer,&xx[ii],&yy[ii],&zz[ii]);
    }
    fclose(fpcoords);
}

int ending()
{
    fprintf(fpouts,"=========================================================\n");
    fprintf(fpouts,"Total energy                %15.6le\n",accumulator[0][5]);
    fprintf(fpouts,"   std. dev.                %15.4lf\n",accumulator[0][6]);
    fprintf(fpouts,"   fluctuation              %15.4lf\n",accumulator[0][7]);

    fprintf(fpouts,"Potentail energy            %15.6le\n",accumulator[1][5]);
    fprintf(fpouts,"   std. dev.                %15.4lf\n",accumulator[1][6]);
    fprintf(fpouts,"   fluctuation              %15.4lf\n",accumulator[1][7]);

    fprintf(fpouts,"Kinetic energy              %15.6le\n",accumulator[2][5]);
    fprintf(fpouts,"   std. dev.                %15.4lf\n",accumulator[2][6]);
    fprintf(fpouts,"   fluctuation              %15.4lf\n",accumulator[2][7]);

    fprintf(fpouts,"Inter potential energy      %15.6le\n",accumulator[3][5]);
    fprintf(fpouts,"   std. dev.                %15.4lf\n",accumulator[3][6]);
    fprintf(fpouts,"   fluctuation              %15.4lf\n",accumulator[3][7]);

    fprintf(fpouts,"Intra potential energy      %15.6le\n",accumulator[4][5]);
    fprintf(fpouts,"   std. dev.                %15.4lf\n",accumulator[4][6]);
    fprintf(fpouts,"   fluctuation              %15.4lf\n",accumulator[4][7]);

    fprintf(fpouts,"LJ energy                   %15.6le\n",accumulator[5][5]);
    fprintf(fpouts,"   std. dev.                %15.4lf\n",accumulator[5][6]);
    fprintf(fpouts,"   fluctuation              %15.4lf\n",accumulator[5][7]);

    fprintf(fpouts,"Bond energy                 %15.6le\n",accumulator[6][5]);
    fprintf(fpouts,"   std. dev.                %15.4lf\n",accumulator[6][6]);
    fprintf(fpouts,"   fluctuation              %15.4lf\n",accumulator[6][7]);

    fprintf(fpouts,"Angle energy                %15.6le\n",accumulator[7][5]);
    fprintf(fpouts,"   std. dev.                %15.4lf\n",accumulator[7][6]);
    fprintf(fpouts,"   fluctuation              %15.4lf\n",accumulator[7][7]);

    fprintf(fpouts,"Dihedral energy             %15.6le\n",accumulator[8][5]);
    fprintf(fpouts,"   std. dev.                %15.4lf\n",accumulator[8][6]);
    fprintf(fpouts,"   fluctuation              %15.4lf\n",accumulator[8][7]);

    fprintf(fpouts,"Improper energy             %15.6le\n",accumulator[9][5]);
    fprintf(fpouts,"   std. dev.                %15.4lf\n",accumulator[9][6]);
    fprintf(fpouts,"   fluctuation              %15.4lf\n",accumulator[9][7]);

    fprintf(fpouts,"Ewald energy                %15.6le\n",accumulator[10][5]);
    fprintf(fpouts,"   std. dev.                %15.4lf\n",accumulator[10][6]);
    fprintf(fpouts,"   fluctuation              %15.4lf\n",accumulator[10][7]);

    fprintf(fpouts,"Real part energy            %15.6le\n",accumulator[11][5]);
    fprintf(fpouts,"   std. dev.                %15.4lf\n",accumulator[11][6]);
    fprintf(fpouts,"   fluctuation              %15.4lf\n",accumulator[11][7]);

    fprintf(fpouts,"Fourier part energy         %15.6le\n",accumulator[12][5]);
    fprintf(fpouts,"   std. dev.                %15.4lf\n",accumulator[12][6]);
    fprintf(fpouts,"   fluctuation              %15.4lf\n",accumulator[12][7]);

    fprintf(fpouts,"Self part energy            %15.6le\n",accumulator[13][5]);
    fprintf(fpouts,"   std. dev.                %15.4lf\n",accumulator[13][6]);
    fprintf(fpouts,"   fluctuation              %15.4lf\n",accumulator[13][7]);

    fprintf(fpouts,"Solid fluid energy          %15.6le\n",accumulator[14][5]);
    fprintf(fpouts,"   std. dev.                %15.4lf\n",accumulator[14][6]);
    fprintf(fpouts,"   fluctuation              %15.4lf\n",accumulator[14][7]);

    fprintf(fpouts,"Temperature                 %15.6le\n",accumulator[15][5]);
    fprintf(fpouts,"   std. dev.                %15.4lf\n",accumulator[15][6]);
    fprintf(fpouts,"   fluctuation              %15.4lf\n",accumulator[15][7]);

    // close files
    fclose(fplog);
    fclose(fptrj);
    fclose(fpouts);
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
    sscanf(fgets(buffer,datalen,fpins), "%lf", &treq);
    sscanf(fgets(buffer,datalen,fpins), "%lf %lf %lf", &boxlx, &boxly, &boxlz);
    sscanf(fgets(buffer,datalen,fpins), "%lf %lf %lf", &rcuton, &rcutoff, &rcutoffelec);
    sscanf(fgets(buffer,datalen,fpins), "%s", coords_file);
    sscanf(fgets(buffer,datalen,fpins), "%d %d", &nstep, &nstep_start);
    sscanf(fgets(buffer,datalen,fpins), "%d %d %d %d %d", 
	    &nstep_ave, &nstep_print, &nstep_save, &nstep_ss, &nstep_trj);
    sscanf(fgets(buffer,datalen,fpins), "%lf %d", &delt, &nstep_inner);
    sscanf(fgets(buffer,datalen,fpins), "%lf", &f0);
    sscanf(fgets(buffer,datalen,fpins), "%d", &isNVTnh);
    sscanf(fgets(buffer,datalen,fpins), "%lf %lf", &qq, &qqs);
    sscanf(fgets(buffer,datalen,fpins), "%d", &isLJswitchOn);
    sscanf(fgets(buffer,datalen,fpins), "%d", &isEwaldOn);
    sscanf(fgets(buffer,datalen,fpins), "%d", &isWolfOn);
    sscanf(fgets(buffer,datalen,fpins), "%d %d %d", &KMAXX, &KMAXY, &KMAXZ);
    sscanf(fgets(buffer,datalen,fpins), "%d", &KSQMAX);
    sscanf(fgets(buffer,datalen,fpins), "%lf", &kappa);

    sscanf(fgets(buffer,datalen,fpins), "%d", &nconstraint);

    sscanf(fgets(buffer,datalen,fpins), "%d", &isSFon);
    sscanf(fgets(buffer,datalen,fpins), "%d", &sf_type);



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
	fscanf(fpcfg,"%lf %lf %lf %lf %d %d\n",&aw[ii],&epsilon[ii],&sigma[ii],
		&charge[ii],&isghost[ii],&tasostype[ii]);

    // read in bond list
    fscanf(fpcfg, "%d bonds\n", &nbond);
    assert(nbond<=nbond_max);
    for (ii=0;ii<nbond;ii++)
	fscanf(fpcfg,"%d %d %d %lf %lf %lf\n",&bond_idx[ii][0],&bond_idx[ii][1],
		&bond_type[ii],&Kb[ii],&Req[ii],&alpha[ii]);

    // read in angle list
    fscanf(fpcfg, "%d angles\n", &nangle);
    assert(nangle<=nangle_max);
    for (ii=0;ii<nangle;ii++)
    {
	fscanf(fpcfg,"%d %d %d %d %lf %lf %lf %lf %lf\n",&angle_idx[ii][0],&angle_idx[ii][1],
		&angle_idx[ii][2],&angle_type[ii],&Ktheta[ii],&Thetaeq[ii],
		&agl_para_3[ii], &agl_para_4[ii], &agl_para_5[ii]); // these 3 parameters only for TRwater
    }

    // read in dihedral list
    fscanf(fpcfg, "%d dihedrals\n", &ndih);
    assert(ndih<=ndih_max);
    for (ii=0;ii<ndih;ii++)
    {
	fscanf(fpcfg,"%d %d %d %d %d %lf %lf %lf %lf\n",&dih_idx[ii][0],&dih_idx[ii][1],
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
    fprintf(fpouts,"=========================================================\n");
    fprintf(fpouts,"Units:\n");
    fprintf(fpouts,"  Energy(J/mol) Temperature(K) Velocity(m/s) Distance(Angstrom)\n");
    fprintf(fpouts,"  Force(J/Angstrom/mol)\n");
    fprintf(fpouts,"=========================================================\n");
    fprintf(fpouts,"natom=%d\n",natom);
    fprintf(fpouts,"nconstraint=%d\n",nconstraint);

    fprintf(fpouts,"%s",sysname);
    fprintf(fpouts,"%d atoms.\n",natom);
    fprintf(fpouts,"%d bonds.\n",nbond);
    fprintf(fpouts,"%d angles.\n",nangle);
    fprintf(fpouts,"%d dihedrals.\n",ndih);
    fprintf(fpouts,"%d impropers.\n",nimp);
    fprintf(fpouts,"%d nonbonded pairs.\n",nnbp);
    fprintf(fpouts,"%d %d %d %d KMAX etc.\n",KMAXX,KMAXY,KMAXZ,KSQMAX);
    fprintf(fpouts,"qq = %le\n",qq);
}

int make_exclude_list()
{
    int ii, jj, nexcllist;
    nexcllist = 0;

    for (ii=0;ii<natom;ii++)
    {
	pointexcl[ii] = nexcllist;
	// exclude bonded atoms
	for (jj=0;jj<nbond;jj++)
	{
	    if (bond_idx[jj][0] == ii)
	    {
		excllist[nexcllist] = bond_idx[jj][1];
		nexcllist++;
	    }
	    else if (bond_idx[jj][1] == ii)
	    {
		excllist[nexcllist] = bond_idx[jj][0];
		nexcllist++;
	    }
	    assert(nexcllist<exclude_max);
	}
	// exclude angled atoms
	for (jj=0;jj<nangle;jj++)
	{
	    if (angle_idx[jj][0] == ii)
	    {
		excllist[nexcllist] = angle_idx[jj][2];
		nexcllist++;
	    }
	    else if (angle_idx[jj][2] == ii)
	    {
		excllist[nexcllist] = angle_idx[jj][0];
		nexcllist++;
	    }
	    assert(nexcllist<exclude_max);
	}
	// exclude dihedraled atoms
	for (jj=0;jj<ndih;jj++)
	{
	    if (dih_idx[jj][0] == ii)
	    {
		excllist[nexcllist] = dih_idx[jj][3];
		nexcllist++;
	    }
	    else if (dih_idx[jj][3] == ii)
	    {
		excllist[nexcllist] = dih_idx[jj][0];
		nexcllist++;
	    }
	    assert(nexcllist<exclude_max);
	}
    }
    pointexcl[natom] = nexcllist;

    // check the exclude list
    /*
       printf("n = %d\n",nexcllist);
       int tmp = 0;
       for (ii=0;ii<natom;ii++)
       {
       printf("%d ",ii);
       for (jj=pointexcl[ii];jj<pointexcl[ii+1];jj++)
       {
       printf("%d ",excllist[tmp]);
       tmp++;
       }
       printf("\n");
       }
     */
}

int velinit()
{
    int ii, jj, kk;
    double px, py, pz;
    double stdvtmp, stdv;
    double totalmass;
    double scaling;

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

int printit()
{
    upot = uinter + uintra;
    // add solid-fluid energy to the potential energy, will be zero if no usf
    upot += usflj;
    utot = upot + ukin;
    // add energy of thermostat, if nose hoover is not used, they will just be zero
    utot = utot + upot_nhts + unhts + unhtss;
    fprintf(stderr,"%10d %10.4le %10.4le %10.4le\n",istep,utot,upot,ukin);
    fprintf(fplog,"%10d %10.4le %10.4le %10.4le %10.4le %10.4le %10.4le %10.4le %10.4le %10.4le %10.4le %10.4le %10.4le %10.4le %10.4le\n",
	    istep,utot,upot,ukin,tinst,uinter,uintra,uvdw,ubond,uangle,udih,uimp,uewald,unhts,unhtss);
}

int vver() // velocity verlet
{
    int ii, ll;

    for (ii=0;ii<natom;ii++)
    {
	// the factor of 1.0e-5 is based on Angstrom (from the force)
	// and femto second (from delt)
	vx[ii] += (deltby2*fxl[ii]*1.0e-5/aw[ii]);
	vy[ii] += (deltby2*fyl[ii]*1.0e-5/aw[ii]);
	vz[ii] += (deltby2*fzl[ii]*1.0e-5/aw[ii]);
    }

    for (ll=0;ll<nstep_inner;ll++)
    {
	for (ii=0;ii<natom;ii++)
	{
	    vx[ii] += (deltsby2*fxs[ii]*1.0e-5/aw[ii]);
	    vy[ii] += (deltsby2*fys[ii]*1.0e-5/aw[ii]);
	    vz[ii] += (deltsby2*fzs[ii]*1.0e-5/aw[ii]);
	    xx[ii] = xx[ii] + delts*vx[ii]*1.0e-5;
	    yy[ii] = yy[ii] + delts*vy[ii]*1.0e-5;
	    zz[ii] = zz[ii] + delts*vz[ii]*1.0e-5;
	}

	// intra forces, short ranged
	rafrc();

	// compute the pseudo velocity at delts
	for (ii=0;ii<natom;ii++)
	{
	    vx[ii] += (deltsby2*fxs[ii]*1.0e-5/aw[ii]);
	    vy[ii] += (deltsby2*fys[ii]*1.0e-5/aw[ii]);
	    vz[ii] += (deltsby2*fzs[ii]*1.0e-5/aw[ii]);
	}
    }

    // inter forces, long ranged
    erfrc();

    // use the new forces to calculate the new velocities at t+delt
    ukin = 0.0;
    for (ii=0;ii<natom;ii++)
    {
	vx[ii] += (deltby2*fxl[ii]*1.0e-5/aw[ii]);
	vy[ii] += (deltby2*fyl[ii]*1.0e-5/aw[ii]);
	vz[ii] += (deltby2*fzl[ii]*1.0e-5/aw[ii]);
	ukin += aw[ii]*(vx[ii]*vx[ii]+vy[ii]*vy[ii]+vz[ii]*vz[ii]);
    }
    ukin = ukin/2.0;

    // calculate instant temperature
    tinst = 2.0*ukin*rRgas/nfree;
}

int snapshot()
{
    int ii;

    fpss = fopen(SNAPSHOT,"w");
    fprintf(fpss,"%d\n",natom);
    fprintf(fpss,"%d %lf %lf %lf\n",istep,boxlx,boxly,boxlz);
    for (ii=0;ii<natom;ii++)
	fprintf(fpss,"C  %lf  %lf  %lf\n",xx[ii],yy[ii],zz[ii]);
    fclose(fpss);
}

int trajectory()
{
    nframe++;
    fwrite(xx,sizeof(double),natom,fptrj);
    fwrite(yy,sizeof(double),natom,fptrj);
    fwrite(zz,sizeof(double),natom,fptrj);
}

int saveit()
{
    int ii;

    fpsave = fopen(SAVEFILE,"wb");

    fwrite(xx,sizeof(double),natom,fpsave);
    fwrite(yy,sizeof(double),natom,fpsave);
    fwrite(zz,sizeof(double),natom,fpsave);
    fwrite(vx,sizeof(double),natom,fpsave);
    fwrite(vy,sizeof(double),natom,fpsave);
    fwrite(vz,sizeof(double),natom,fpsave);

    fwrite(&istep,sizeof(int),1,fpsave);
    fwrite(icounter,sizeof(int),num_counter_max,fpsave);
    fwrite(accumulator,sizeof(double),num_counter_max,fpsave);

    fclose(fpsave);
}

int averages()
{
    int ii;
    double temp1;

    for (ii=0;ii<num_counter_max;ii++)
    {
	temp1 = accumulator[ii][0]/nstep_ave;
	accumulator[ii][2] += temp1;
	accumulator[ii][3] += accumulator[ii][1]/nstep_ave;
	accumulator[ii][4] += temp1*temp1;
	// rezero
	accumulator[ii][0] = 0.0;
	accumulator[ii][1] = 0.0;
    }
    icounter[10]++; // number of average cycles
}

int calres()
{
    int ii;
    double ave_of_square, ave_of_ave_square;
    double ave, err, fluc;
    for (ii=0;ii<num_counter_max;ii++)
    {
	ave = accumulator[ii][5] = accumulator[ii][2]/icounter[10]; // ave
	ave_of_square = accumulator[ii][3]/icounter[10];
	ave_of_ave_square = accumulator[ii][4]/icounter[10];
	err = accumulator[ii][6] = sqrt(fabs(ave_of_ave_square-ave*ave)); // err
	fluc = accumulator[ii][7] = sqrt(fabs(ave_of_square-ave*ave)); // fluc
	// rezero?
    }
}

int main (int argc, char *argv[])
{
    readins();
    init_vars();
    echo();
    make_exclude_list();

    // start of the MD simulation
    // initialize velocities
    velinit(); 
    
    // loadit();

    // calculate total energies
    erfrc(); 
    rafrc();
    // print initial properties
    printit();
    // make snapshots & movies
    trajectory();

    for (istep=nstep_start;istep<=nstep;istep++) // NOTE: start from 1 and <=
    {
	if (isNVTnh)
	    vver_nh(); // vv with nose hoover
	else 
	    vver(); // velocity verlet

	// accumulators
	accumulator[0][0] += utot;
	accumulator[0][1] += utot*utot;
	accumulator[1][0] += upot;
	accumulator[1][1] += upot*upot;
	accumulator[2][0] += ukin;
	accumulator[2][1] += ukin*upot;
	accumulator[3][0] += uinter;
	accumulator[3][1] += uinter*uinter;
	accumulator[4][0] += uintra;
	accumulator[4][1] += uintra*uintra;
	accumulator[5][0] += uvdw;
	accumulator[5][1] += uvdw*uvdw;
	accumulator[6][0] += ubond;
	accumulator[6][1] += ubond*ubond;
	accumulator[7][0] += uangle;
	accumulator[7][1] += uangle*uangle;
	accumulator[8][0] += udih;
	accumulator[8][1] += udih*udih;
	accumulator[9][0] += uimp;
	accumulator[9][1] += uimp*uimp;
	accumulator[10][0] += uewald;
	accumulator[10][1] += uewald*uewald;
	accumulator[11][0] += ureal;
	accumulator[11][1] += ureal*ureal;
	accumulator[12][0] += ufourier;
	accumulator[12][1] += ufourier*ufourier;
	accumulator[13][0] += uself;
	accumulator[13][1] += uself*uself;
	accumulator[14][0] += usflj;
	accumulator[14][1] += usflj*usflj;
	accumulator[15][0] += tinst;
	accumulator[15][1] += tinst*tinst;

	if (istep%nstep_print == 0) printit();

	if (nstep_ss && istep%nstep_ss == 0) snapshot();

	if (nstep_trj && istep%nstep_trj==0) trajectory();

	if (istep%nstep_ave==0) averages();

	if (istep%nstep_save==0) saveit();

    }

    snapshot();

    calres();

    fprintf(stderr,"%d frames in the trajectory file.\n",nframe);

    // write out results and clean up
    ending();

}

