/*
 * Maintainable Simplex Purpose Molecular Simulator 2
 * Rewritten with standard C language
 * The goal is simpler, quicker, easier, better.
 * No class and other complex data structures.
 * Use common names for input files. So the every job must
 * have its own working directory.
 * 
 * 
 * Written by Yang Wang 2007
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
extern int vver_nh_1();
extern int vver_nh_2();
extern int vver_nh_3();
extern int init_sf_hypergeo();
extern int init_sf_atom_explicit();
extern int init_tasos_grid();
extern int init_my_interp();


/* Initiate variables */
int init_vars()
{
    int ii, jj;
    FILE *fpcoords;
    char buffer[200];
    const int datalen = 200;

    // initiate files
    // output file is initialized already at the very beginning of the run
    fplog = fopen(LOGFILE,"w");
    fptrj = fopen(MOVIE,"wb");

    fprintf(stderr,"Initializing variables...\n");
    fprintf(fpouts,"Initializing variables...\n");
    for (ii=0;ii<nmole;ii++)
    {
	mw[ii] = 0.0;
	for (jj=mole_first_atom_idx[ii];jj<mole_first_atom_idx[ii+1];jj++)
	    mw[ii] += aw[jj];
	// fprintf(fpouts,"molecule %d: %d-%d :  mw=%lf : 1st bond id=%d : 1st angle id=%d : 1st dih id=%d\n",
	// ii,mole_first_atom_idx[ii],mole_first_atom_idx[ii+1]-1,mw[ii],
	// mole_first_bond_idx[ii],mole_first_angle_idx[ii],mole_first_dih_idx[ii]);
    }

    nframe = 0; // number of frames in trajectory file

    // set the starting step to 1, will be changed by load it if it is continue run
    nstep_start = 1;

    for (ii=0;ii<num_counter_max;ii++)
    {
	icounter[ii] = 0;
	for (jj=0;jj<5;jj++)
	    accumulator[ii][jj] = 0.0;
    }

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
    if (isSFon)
    {
	if (sf_type==nanotube_hypergeo)
	    init_sf_hypergeo();
	else if (sf_type==nanotube_atom_explicit)
	    init_sf_atom_explicit();
	else if (sf_type==nanotube_tasos)
	    init_tasos_grid();
	else if (sf_type==nanotube_my_interp)
	    init_my_interp();
    }

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
	    for (jj=ii+1;jj<ndih;jj++)
	    {
		if (isDih_unique[jj])
		{
		    if (dih_idx[ii][0]==dih_idx[jj][0] && dih_idx[ii][3]==dih_idx[jj][3])
		    {
			isDih_unique[jj] = false;
			fprintf(fpouts,"Dihedral %d and dihedral %d have the same ending pairs.\n",ii,jj);
		    }
		    else if (dih_idx[ii][0]==dih_idx[jj][3] && dih_idx[ii][3]==dih_idx[jj][0])
		    {
			isDih_unique[jj] = false;
			fprintf(fpouts,"Dihedral %d and dihedral %d have the same ending pairs.\n",ii,jj);
		    }
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
		{
		    isDih_unique[jj] = false; 
		    fprintf(fpouts,"Dihedral %d and angle %d have the same ending pairs.\n",ii,jj);
		}
		else if (dih_idx[ii][0]==angle_idx[jj][2] && dih_idx[ii][3]==angle_idx[jj][0])
		{
		    isDih_unique[jj] = false;
		    fprintf(fpouts,"Dihedral %d and angle %d have the same ending pairs.\n",ii,jj);
		}
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
		{
		    isDih_unique[jj] = false;
		    fprintf(fpouts,"Dihedral %d and bond %d have the same ending pairs.\n",ii,jj);
		}
		else if (dih_idx[ii][0]==bond_idx[jj][1] && dih_idx[ii][3]==bond_idx[jj][0])
		{
		    isDih_unique[jj] = false;
		    fprintf(fpouts,"Dihedral %d and bond %d have the same ending pairs.\n",ii,jj);
		}
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
	    for (jj=ii+1;jj<nangle;jj++)
	    {
		if (isAngle_unique[jj])
		{
		    if (angle_idx[ii][0]==angle_idx[jj][0] && angle_idx[ii][2]==angle_idx[jj][2])
		    {
			isAngle_unique[jj] = false;
			fprintf(fpouts,"Angle %d and angle %d have the same ending pairs.\n",ii,jj);
		    }
		    else if (angle_idx[ii][0]==angle_idx[jj][2] && angle_idx[ii][2]==angle_idx[jj][0])
		    {
			isAngle_unique[jj] = false;
			fprintf(fpouts,"Angle %d and angle %d have the same ending pairs.\n",ii,jj);
		    }
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
		{
		    isAngle_unique[jj] = false;
		    fprintf(fpouts,"Angle %d and bond %d have the same ending pairs.\n",ii,jj);
		}
		else if (angle_idx[ii][0]==bond_idx[jj][1] && angle_idx[ii][2]==bond_idx[jj][0])
		{
		    isAngle_unique[jj] = false;
		    fprintf(fpouts,"Angle %d and bond %d have the same ending pairs.\n",ii,jj);
		}
	    }
	}
    }

    // set ewald parameters
    kappasq = kappa*kappa;
    Bfactor_ewald = 1.0/(4.0*kappa*kappa);
    Vfactor_ewald = 2.0*pi/(boxlx*boxly*boxlz);
    TWOPI_LX = 2.0*pi/boxlx;
    TWOPI_LY = 2.0*pi/boxly;
    TWOPI_LZ = 2.0*pi/boxlz;
    // 1D ewald constant
    twopi_over_3v = 2.0*pi/3.0/boxlx/boxly/boxlz;

    // set wolf parameters
    wolfvcon1 = -erfc(kappa*rcutoffelec)/rcutoffelec;
    wolfvcon2 = erfc(kappa*rcutoffelec)/rcutoffelecsq 
	+ 2.0*kappa*exp(-(kappa*rcutoffelec)*(kappa*rcutoffelec))/(sqrt(pi)*rcutoffelec);
    wolffcon1 = 2.0*kappa/sqrt(pi);
    wolffcon2 = -wolfvcon2;


    // read in coordinates
    fpcoords = fopen(COORDSIN,"r");
    // fscanf(fpcoords, "%[^\n]", buffer);
    fgets(buffer,datalen,fpcoords);
    fgets(buffer,datalen,fpcoords);
    for (ii=0;ii<natom;ii++)
    {
	fscanf(fpcoords,"%s %lf %lf %lf\n",buffer,&xx[ii],&yy[ii],&zz[ii]);
    }
    fclose(fpcoords);
}

int ending()
{
    int ii;

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

    fprintf(fpouts,"Vaccum energy               %15.6le\n",accumulator[16][5]);
    fprintf(fpouts,"   std. dev.                %15.4lf\n",accumulator[16][6]);
    fprintf(fpouts,"   fluctuation              %15.4lf\n",accumulator[16][7]);

    fprintf(fpouts,"Solid fluid energy          %15.6le\n",accumulator[14][5]);
    fprintf(fpouts,"   std. dev.                %15.4lf\n",accumulator[14][6]);
    fprintf(fpouts,"   fluctuation              %15.4lf\n",accumulator[14][7]);

    fprintf(fpouts,"Temperature                 %15.6le\n",accumulator[15][5]);
    fprintf(fpouts,"   std. dev.                %15.4lf\n",accumulator[15][6]);
    fprintf(fpouts,"   fluctuation              %15.4lf\n",accumulator[15][7]);

    fprintf(fpouts,"Wolf energy                 %15.6le\n",accumulator[17][5]);
    fprintf(fpouts,"   std. dev.                %15.4lf\n",accumulator[17][6]);
    fprintf(fpouts,"   fluctuation              %15.4lf\n",accumulator[17][7]);

    fprintf(fpouts,"=========================================================\n");

    if (isSFon)
    {
	if (sf_type==nanotube_hypergeo)
	{
	    free(hgntc_xx);
	    free(hgntc_yy);
	    free(hgnt_radius);
	}
	else if (sf_type==nanotube_atom_explicit)
	{
	    free(solid_sigma);
	    free(solid_epsilon);
	    free(solid_charge);
	    free(solid_xx);
	    free(solid_yy);
	    free(solid_zz);
	}
	else if (sf_type==nanotube_my_interp)
	{
	    // free memories
	    free(interp_vector);
	    for (ii=0;ii<nunique_atom_max;ii++)
	    {
		free(ene0[ii]);
		free(ene1[ii]);
		free(ene2[ii]);
		free(ene3[ii]);
		free(ene4[ii]);
		free(ene5[ii]);
		free(ene6[ii]);
		free(ene7[ii]);
		free(ene8[ii]);
		free(ene9[ii]);
		free(ene10[ii]);
		free(ene11[ii]);
		free(ene12[ii]);
		free(ene13[ii]);
		free(ene14[ii]);
		free(ene15[ii]);
		free(ene16[ii]);
		free(ene17[ii]);
		free(ene18[ii]);
		free(ene19[ii]);
		free(ene20[ii]);
		free(ene21[ii]);
		free(ene22[ii]);
		free(ene23[ii]);
		free(ene24[ii]);
		free(ene25[ii]);
		free(ene26[ii]);
		free(ene27[ii]);
		free(ene28[ii]);
		free(ene29[ii]);
		free(ene30[ii]);
		free(ene31[ii]);

		free(fxa0[ii]);
		free(fxa1[ii]);
		free(fxa2[ii]);
		free(fxa3[ii]);
		free(fxa4[ii]);
		free(fxa5[ii]);
		free(fxa6[ii]);
		free(fxa7[ii]);
		free(fxa8[ii]);
		free(fxa9[ii]);
		free(fxa10[ii]);
		free(fxa11[ii]);
		free(fxa12[ii]);
		free(fxa13[ii]);
		free(fxa14[ii]);
		free(fxa15[ii]);
		free(fxa16[ii]);
		free(fxa17[ii]);
		free(fxa18[ii]);
		free(fxa19[ii]);
		free(fxa20[ii]);
		free(fxa21[ii]);
		free(fxa22[ii]);
		free(fxa23[ii]);
		free(fxa24[ii]);
		free(fxa25[ii]);
		free(fxa26[ii]);
		free(fxa27[ii]);
		free(fxa28[ii]);
		free(fxa29[ii]);
		free(fxa30[ii]);
		free(fxa31[ii]);

		free(fya0[ii]);
		free(fya1[ii]);
		free(fya2[ii]);
		free(fya3[ii]);
		free(fya4[ii]);
		free(fya5[ii]);
		free(fya6[ii]);
		free(fya7[ii]);
		free(fya8[ii]);
		free(fya9[ii]);
		free(fya10[ii]);
		free(fya11[ii]);
		free(fya12[ii]);
		free(fya13[ii]);
		free(fya14[ii]);
		free(fya15[ii]);
		free(fya16[ii]);
		free(fya17[ii]);
		free(fya18[ii]);
		free(fya19[ii]);
		free(fya20[ii]);
		free(fya21[ii]);
		free(fya22[ii]);
		free(fya23[ii]);
		free(fya24[ii]);
		free(fya25[ii]);
		free(fya26[ii]);
		free(fya27[ii]);
		free(fya28[ii]);
		free(fya29[ii]);
		free(fya30[ii]);
		free(fya31[ii]);

		free(fza0[ii]);
		free(fza1[ii]);
		free(fza2[ii]);
		free(fza3[ii]);
		free(fza4[ii]);
		free(fza5[ii]);
		free(fza6[ii]);
		free(fza7[ii]);
		free(fza8[ii]);
		free(fza9[ii]);
		free(fza10[ii]);
		free(fza11[ii]);
		free(fza12[ii]);
		free(fza13[ii]);
		free(fza14[ii]);
		free(fza15[ii]);
		free(fza16[ii]);
		free(fza17[ii]);
		free(fza18[ii]);
		free(fza19[ii]);
		free(fza20[ii]);
		free(fza21[ii]);
		free(fza22[ii]);
		free(fza23[ii]);
		free(fza24[ii]);
		free(fza25[ii]);
		free(fza26[ii]);
		free(fza27[ii]);
		free(fza28[ii]);
		free(fza29[ii]);
		free(fza30[ii]);
		free(fza31[ii]);

	    }
	}
    }

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
    int atomid, moleid;
    int isFirstAtom;
    int last_mole_id;

    fprintf(stderr,"Reading input file...\n");
    fprintf(fpouts,"Reading input file...\n");

    /* read input file */
    fpins = fopen(INPUT,"r");

    sscanf(fgets(buffer,datalen,fpins), "%d %d", &ij, &jk);
    sscanf(fgets(buffer,datalen,fpins), "%lf", &treq);
    sscanf(fgets(buffer,datalen,fpins), "%lf %lf %lf", &boxlx, &boxly, &boxlz);
    sscanf(fgets(buffer,datalen,fpins), "%lf %lf %lf", &rcuton, &rcutoff, &rcutoffelec);
    // sscanf(fgets(buffer,datalen,fpins), "%s", coords_file);
    sscanf(fgets(buffer,datalen,fpins), "%d %d", &nstep, &fStart_option);
    sscanf(fgets(buffer,datalen,fpins), "%d %d %d %d %d", 
	    &nstep_ave, &nstep_print, &nstep_save, &nstep_ss, &nstep_trj);
    sscanf(fgets(buffer,datalen,fpins), "%lf %d", &delt, &nstep_inner);
    sscanf(fgets(buffer,datalen,fpins), "%lf", &f0);
    sscanf(fgets(buffer,datalen,fpins), "%d %d", &isNVTnh, &whichNH);
    sscanf(fgets(buffer,datalen,fpins), "%lf %lf", &qq, &qqs);
    sscanf(fgets(buffer,datalen,fpins), "%d", &isLJswitchOn);
    sscanf(fgets(buffer,datalen,fpins), "%d %d %d", &isEwaldOn,&fEwald_BC,&fEwald_Dim);
    sscanf(fgets(buffer,datalen,fpins), "%d", &isWolfOn);
    sscanf(fgets(buffer,datalen,fpins), "%d %d %d", &KMAXX, &KMAXY, &KMAXZ);
    sscanf(fgets(buffer,datalen,fpins), "%d", &KSQMAX);
    sscanf(fgets(buffer,datalen,fpins), "%lf", &kappa);

    sscanf(fgets(buffer,datalen,fpins), "%d", &nconstraint);

    sscanf(fgets(buffer,datalen,fpins), "%d", &isSFon);
    sscanf(fgets(buffer,datalen,fpins), "%d", &sf_type);

    fclose(fpins);

    fprintf(stderr,"reading cfg file...\n");
    fprintf(fpouts,"reading cfg file...\n");
    
    /* read config file */
    fpcfg = fopen(CONFIG,"r");
    // read the first line of comments
    // fscanf(fpcfg, "%[^\n]", buffer); /* read in everything to the buffer except the return */
    fscanf(fpcfg,"%s %d\n",sysname,&nspecie);
    assert(nspecie<=nspecie_max);

    // read in atom list
    fscanf(fpcfg, "%d atoms\n", &natom);
    assert(natom<=natom_max);
    nmole = 0;
    isFirstAtom = 1;
    last_mole_id = 0;
    for (ii=0;ii<natom;ii++)
    {
	fscanf(fpcfg,"%d %d %s %lf %lf %lf %lf %d %d\n",&atomid,&moleid,atomname[ii],&aw[ii],&epsilon[ii],&sigma[ii],
		&charge[ii],&isghost[ii],&tasostype[ii]);
	// assign the molecule index to the atom
	atom_to_mole_idx[ii] = moleid;
	// assign the correct index of the first atom in a molecule
	if (moleid == last_mole_id+1)
	{
	    isFirstAtom = 1;
	    last_mole_id = moleid;
	}
	if (isFirstAtom)
	{
	    mole_first_atom_idx[nmole] = atomid;
	    nmole++;
	    assert(nmole<nmole_max+1);
	    isFirstAtom = 0;
	}
    }
    mole_first_atom_idx[nmole] = atomid + 1; // set the boundary of the last molecule

    // read in bond list
    fscanf(fpcfg, "%d bonds\n", &nbond);
    assert(nbond<=nbond_max);
    last_mole_id = 0;
    for (ii=0;ii<nbond;ii++)
    {
	fscanf(fpcfg,"%d %d %d %lf %lf %lf\n",&bond_idx[ii][0],&bond_idx[ii][1],
		&bond_type[ii],&Kb[ii],&Req[ii],&alpha[ii]);
	// check if the bond atom is belong to the molecule
	// if it is, set this bond as the 1st bond of the molecule and increase
	// the counter of the molecule by 1 to the next molecule
	if (mole_first_atom_idx[last_mole_id]<=bond_idx[ii][0] && bond_idx[ii][0]< mole_first_atom_idx[last_mole_id+1]) 
	{
	    mole_first_bond_idx[last_mole_id] = ii;
	    last_mole_id++;
	    assert(last_mole_id<nmole_max+1);
	}
    }
    // set the upper bound for the last molecule
    mole_first_bond_idx[last_mole_id] = ii;

    // read in angle list
    fscanf(fpcfg, "%d angles\n", &nangle);
    assert(nangle<=nangle_max);
    last_mole_id = 0;
    for (ii=0;ii<nangle;ii++)
    {
	fscanf(fpcfg,"%d %d %d %d %lf %lf %lf %lf %lf\n",&angle_idx[ii][0],&angle_idx[ii][1],
		&angle_idx[ii][2],&angle_type[ii],&Ktheta[ii],&Thetaeq[ii],
		&agl_para_3[ii], &agl_para_4[ii], &agl_para_5[ii]); // these 3 parameters only for TRwater

	// check if the angle atom is belong to the molecule
	// if it is, set this angle as the 1st angle of the molecule and increase
	// the counter of the molecule by 1 to the next molecule
	if (mole_first_atom_idx[last_mole_id]<=angle_idx[ii][0] && angle_idx[ii][0]< mole_first_atom_idx[last_mole_id+1]) 
	{
	    mole_first_angle_idx[last_mole_id] = ii;
	    last_mole_id++;
	    assert(last_mole_id<nmole_max+1);
	}
    }
    // set the upper bound for the last molecule
    mole_first_angle_idx[last_mole_id] = ii;

    // read in dihedral list
    fscanf(fpcfg, "%d dihedrals\n", &ndih);
    assert(ndih<=ndih_max);
    last_mole_id = 0;
    for (ii=0;ii<ndih;ii++)
    {
	fscanf(fpcfg,"%d %d %d %d %d %lf %lf %lf %lf\n",&dih_idx[ii][0],&dih_idx[ii][1],
		&dih_idx[ii][2],&dih_idx[ii][3],&dih_type[ii],&c1[ii],&c2[ii],&c3[ii],&c4[ii]);
	// check if the dihedral atom is belong to the molecule
	// if it is, set this diheral as the 1st dihedral of the molecule and increase
	// the counter of the molecule by 1 to the next molecule
	if (mole_first_atom_idx[last_mole_id]<=dih_idx[ii][0] && dih_idx[ii][0]< mole_first_atom_idx[last_mole_id+1]) 
	{
	    mole_first_dih_idx[last_mole_id] = ii;
	    last_mole_id++;
	    assert(last_mole_id<nmole_max+1);
	}
    }
    // set the upper bound for the last molecule
    mole_first_dih_idx[last_mole_id] = ii;

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
    fprintf(fpouts,"=====================================\n");
    fprintf(fpouts,"Units:\n");
    fprintf(fpouts,"  Energy(J/mol)\n  Temperature(K)\n  Velocity(m/s)\n  Distance(Angstrom)\n");
    fprintf(fpouts,"  Force(J/Angstrom/mol)\n");
    fprintf(fpouts,"=====================================\n");
    fprintf(fpouts,"%s\n",sysname);
    fprintf(fpouts,"natom=%d\n",natom);
    fprintf(fpouts,"nmole = %d\n",nmole);
    fprintf(fpouts,"nconstraint=%d\n",nconstraint);

    fprintf(fpouts,"%d atoms.\n",natom);
    fprintf(fpouts,"%d bonds.\n",nbond);
    fprintf(fpouts,"%d angles.\n",nangle);
    fprintf(fpouts,"%d dihedrals.\n",ndih);
    fprintf(fpouts,"%d impropers.\n",nimp);
    fprintf(fpouts,"%d nonbonded pairs.\n",nnbp);
    fprintf(fpouts,"%lf kappa\n",kappa);
    fprintf(fpouts,"%d %d %d %d KMAX etc.\n",KMAXX,KMAXY,KMAXZ,KSQMAX);
    // nvt
    // charge on
    // sf on
    // LJ switch on
    // copy in.mspms here?

    fprintf(fpouts,"utot=%lf\n",uinter+uintra+ukin); // utot
    fprintf(fpouts,"upot=%lf\n",uinter+uintra); // upot
    fprintf(fpouts,"ukin=%lf\n",ukin);
    fprintf(fpouts,"uinter=%lf\n",uinter);
    fprintf(fpouts,"uintra=%lf\n",uintra);
    fprintf(fpouts,"uvdw=%lf\n",uvdw);
    fprintf(fpouts,"ubond=%lf\n",ubond);
    fprintf(fpouts,"uangle=%lf\n",uangle);
    fprintf(fpouts,"udih=%lf\n",udih);
    fprintf(fpouts,"uimp=%lf\n",uimp);
    fprintf(fpouts,"uewald=%lf\n",uewald);
    fprintf(fpouts,"ureal=%lf\n",ureal);
    fprintf(fpouts,"ufourier=%lf\n",ufourier);
    fprintf(fpouts,"uself=%lf\n",uself);
    fprintf(fpouts,"uexcl=%lf\n",uexcl);
    fprintf(fpouts,"uvaccum=%lf\n",uvaccum);
    fprintf(fpouts,"uGz0=%lf\n",uGz0);

    fprintf(fpouts,"usflj=%lf\n",usflj);
    fprintf(fpouts,"uwolf=%lf\n",uwolf);
    fprintf(fpouts,"uwolf_real=%lf\n",uwolf_real);
    fprintf(fpouts,"uwolf_con=%lf\n",uwolf_con);

    fflush(fpouts);
}

int make_exclude_list()
{
    int ii, jj, nexcllist;

    fprintf(fpouts,"making exclude list...\n");

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
       printf("%d excludes are:\n",ii);
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

    fprintf(fpouts,"initializing velocities...\n");

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
    utot = upot + ukin;
    // add energy of thermostat, if nose hoover is not used, they will just be zero
    utot = utot + upot_nhts + unhts + unhtss;
    fprintf(stderr,"%10d %10.4le %10.4le %10.4le %10.4le\n",istep,utot,upot,ukin,tinst);
    fprintf(fplog,"%10d %10.4le %10.4le %10.4le %10.4le %10.4le %10.4le %10.4le %10.4le %10.4le %10.4le %10.4le %10.4le %10.4le %10.4le %10.4le\n",
	    istep,utot,upot,ukin,tinst,uinter,uintra,uvdw,ubond,uangle,udih,uimp,uewald,usflj,unhts,unhtss);
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
	fprintf(fpss,"%s  %lf  %lf  %lf\n",atomname[ii],xx[ii],yy[ii],zz[ii]);
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

    fwrite(&qq,sizeof(double),1,fpsave);
    fwrite(&ps,sizeof(double),1,fpsave);
    fwrite(&gg,sizeof(double),1,fpsave);
    fwrite(&ss,sizeof(double),1,fpsave);
    fwrite(&qqs,sizeof(double),1,fpsave);
    fwrite(&pss,sizeof(double),1,fpsave);
    fwrite(&ggs,sizeof(double),1,fpsave);
    fwrite(&sss,sizeof(double),1,fpsave);

    fwrite(&istep,sizeof(int),1,fpsave);
    fwrite(icounter,sizeof(int),num_counter_max,fpsave);
    fwrite(accumulator,sizeof(double),num_counter_max,fpsave);

    fclose(fpsave);
}

int loadit()
{
    int ii;

    fprintf(stderr,"loading from saved file...\n");
    fprintf(fpouts,"loading from saved file...\n");

    fpload = fopen(LOADFILE,"rb");

    fread(xx,sizeof(double),natom,fpload);
    fread(yy,sizeof(double),natom,fpload);
    fread(zz,sizeof(double),natom,fpload);
    fread(vx,sizeof(double),natom,fpload);
    fread(vy,sizeof(double),natom,fpload);
    fread(vz,sizeof(double),natom,fpload);

    fread(&qq,sizeof(double),1,fpload);
    fread(&ps,sizeof(double),1,fpload);
    fread(&gg,sizeof(double),1,fpload);
    fread(&ss,sizeof(double),1,fpload);
    fread(&qqs,sizeof(double),1,fpload);
    fread(&pss,sizeof(double),1,fpload);
    fread(&ggs,sizeof(double),1,fpload);
    fread(&sss,sizeof(double),1,fpload);

    fread(&istep,sizeof(int),1,fpload);
    fread(icounter,sizeof(int),num_counter_max,fpload);
    fread(accumulator,sizeof(double),num_counter_max,fpload);

    fclose(fpload);
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

int opening()
{
    // open output file at the very beginning to keep log of the run
    fpouts = fopen(OUTPUT,"w");

    fprintf(fpouts,"WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW\n");
    fprintf(fpouts,"WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW\n");
    fprintf(fpouts,"WW......WWW.....WW.....WWW......WWW.....WWW\n");
    fprintf(fpouts,"WW..W.W..W..WWWWWW..WW..WW..W.W..W..WWWWWWW\n");
    fprintf(fpouts,"WW..W.W..W..WWWWWW..WW..WW..W.W..W..WWWWWWW\n");
    fprintf(fpouts,"WW..W.W..WW....WWW..WW..WW..W.W..WW....WWWW\n");
    fprintf(fpouts,"WW..W.W..WWWWW..WW..WW..WW..W.W..WWWWW..WWW\n");
    fprintf(fpouts,"WW..W.W..WWWWW..WW..WW..WW..W.W..WWWWW..WWW\n");
    fprintf(fpouts,"WW..WWW..W.....WWW.....WWW..WWW..W.....WW-2\n");
    fprintf(fpouts,"WWWWWWWWWWWWWWWWWW..WWWWWWWWWWWWWWWWWWWWWWW\n");
    fprintf(fpouts,"WWWWWWWWWWWWWWWWWW..WWWWWWWWWWWWWWWWWWWWWWW\n");
    fprintf(fpouts,"WWWWWWWWWWWWWWWWWW..WWWWWWWWWWWWWWWWWWWWWWW\n");
    fprintf(fpouts,"WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWYANG\n");
}

extern int get_values_from_grid(double, double, double, int, double*, double*);
int main (int argc, char *argv[])
{
    opening();
    readins();
    init_vars();
    make_exclude_list();

    // start of the MD simulation
    // initialize velocities
    velinit(); 

    if (fStart_option!=new_run) loadit(); // if not new run, load from old file

    // calculate total energies
    erfrc(); 
    rafrc();
    // print out initial values
    echo();
    // print initial properties
    printit();
    // make snapshots & movies
    trajectory();

    for (istep=nstep_start;istep<=nstep;istep++) // NOTE: start from 1 and <=
    {
	if (isNVTnh) // velocity verlet with nose hoover
	{
	    switch (whichNH)
	    {
		case 1:
		    vver_nh_1();
		    break;
		case 2:
		    vver_nh_2();
		    break;
		case 3:
		    vver_nh_3();
		    break;
	    }	
	}
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
	accumulator[16][0] += uvaccum;
	accumulator[16][1] += uvaccum*uvaccum;
	accumulator[17][0] += uwolf;
	accumulator[17][1] += uwolf*uwolf;

	if (istep%nstep_print == 0) printit();

	if (nstep_ss && istep%nstep_ss == 0) snapshot();

	if (nstep_trj && istep%nstep_trj==0) trajectory();

	if (istep%nstep_ave==0) averages();

	if (istep%nstep_save==0) saveit();

    }

    snapshot();

    calres();

    fprintf(stderr,"%d frames in the trajectory file.\n",nframe);
    fprintf(fpouts,"%d frames in the trajectory file.\n",nframe);

    // write out results and clean up
    ending();

}

