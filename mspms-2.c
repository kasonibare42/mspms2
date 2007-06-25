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
#include <ctype.h>
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
extern int init_npt_respa();
extern int npt_respa();


int calculate_ljlrc();

/* Initiate variables */
int init_vars()
{
    int ii, jj;
    FILE *fpcoords;
    char buffer[200];
    const int datalen = 200;

    fprintf(stderr,"Initializing variables...\n");
    fprintf(fpouts,"Initializing variables...\n");

    // initiate files
    // output file is initialized already at the very beginning of the run
    fplog = fopen(LOGFILE,"w");
    fptrj = fopen(MOVIE,"wb");

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
    istep = 0; // istep is used for printit, the first print should be at step zero
    nstep_start = 1;

    // initiate counters and accumulators
    for (ii=0;ii<num_counter_max;ii++)
    {
	icounter[ii] = 0;
	for (jj=0;jj<5;jj++)
	    accumulator[ii][jj] = 0.0;
    }
    // set the counter for equilibrium
    // it will decrease during the run
    icounter[11] = nstep_eq;

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

    // volume
    boxv = boxlx*boxly*boxlz;

    // delt
    deltby2 = delt/2.0;
    delts = delt/nstep_inner;
    deltsby2 = delts/2.0;

    dt_outer2 = deltby2;
    dt_outer4 = delt/4.0;
    dt_outer8 = delt/8.0;

    // initialize thermostat/baron stat input data
    if (what_ensemble == npt_run)
	init_npt_respa();

    // nose hoover
    // following for Dr. Maginn's nose hoover
    // Gts = 0.0;
    // vts = 0.0;
    // rts = 0.0;
    /*
     * Qts = Rgas*treq*nfree/Omega
     * where Omega is a parameter related to the mass of the thermostat
     * for this program, we read in the Qts directly.
     */
    
    // variables for nose-hoover NVT, see frenkel and smit
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
    fprintf(stderr,"reading initial coordinates of the system...\n");
    fprintf(fpouts,"reading initial coordinates of the system...\n");
    // fscanf(fpcoords, "%[^\n]", buffer);
    fgets(buffer,datalen,fpcoords);
    fgets(buffer,datalen,fpcoords);
    for (ii=0;ii<natom;ii++)
    {
	fscanf(fpcoords,"%s %lf %lf %lf\n",buffer,&xx[ii],&yy[ii],&zz[ii]);
    }
    fclose(fpcoords);

    double uljlrc_term1, uljlrc_term2;
    double pljlrc_term1, pljlrc_term2;
    int mm, nn;
    int atomid_1, atomid_2;
    double sigmaij, epsilonij;
    int isSameSpecie;
    double temp1, temp2, temp3;

    fprintf(stderr,"calculating LJ long range corrections...\n");
    fprintf(fpouts,"calculating LJ long range corrections...\n");

    uljlrc_term1 = (8.0/9.0)*pi*pow(rcutoff,-9.0);
    uljlrc_term2 = -(8.0/3.0)*pi*pow(rcutoff,-3.0); 
    pljlrc_term1 = (32.0/9.0)*pi*pow(rcutoff,-9.0);
    pljlrc_term2 = -(16.0/3.0)*pi*pow(rcutoff,-3.0);
    // calculate long range correction terms for single molecules
    for (mm=0;mm<nspecie;mm++)
    {
	for (nn=0;nn<nspecie;nn++)
	{
	    uljlrc_term[mm][nn] = 0.0;
	    pljlrc_term[mm][nn] = 0.0;
	    // loop through the atoms in one molecule of specie mm
	    for (ii=0;ii<natom_per_mole[mm];ii++)
	    {
		// use the first molecule of one specie to do the calculation
		atomid_1 = specie_first_atom_idx[mm] + ii;
		// loop through all the atoms in one molecule of specie nn
		for (jj=0;jj<natom_per_mole[nn];jj++)
		{
		    atomid_2 = specie_first_atom_idx[nn] + ii;
		    sigmaij = 0.5*(sigma[atomid_1]+sigma[atomid_2]);
		    epsilonij = sqrt(epsilon[atomid_1]*epsilon[atomid_2]);
		    temp1 = pow(sigmaij,9.0)*uljlrc_term1 + pow(sigmaij,3.0)*uljlrc_term2;
		    temp2 = epsilonij*pow(sigmaij,3.0);
		    temp3 = pow(sigmaij,9.0)*pljlrc_term1 + pow(sigmaij,3.0)*pljlrc_term2;
		    uljlrc_term[mm][nn] += temp1*temp2;
		    pljlrc_term[mm][nn] += temp3*temp2;

		} // through all atom in one molecule of specie 2
	    } // through all atom in one molecule of specie 1

	} // loop through specie 2
    } // loop through speice 1

    // calculate the total lj lrc
    calculate_ljlrc();

}

int calculate_ljlrc()
{
    int mm, nn;

    // calculate the lrc for the real system
    uljlrc = 0.0;
    pljlrc = 0.0;
    for (mm=0;mm<nspecie;mm++)
    {
	for (nn=0;nn<nspecie;nn++)
	{
	    uljlrc += (uljlrc_term[mm][nn]*nmole_per_specie[mm]*nmole_per_specie[nn]/boxv);
	    pljlrc += 
		(pljlrc_term[mm][nn]*nmole_per_specie[mm]*nmole_per_specie[nn]/boxv/boxv);
	}
    }
    pljlrc = pljlrc*J_mol_A3_to_Pascal; // convert J/mol/A^3 to Pascal

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

    fprintf(fpouts,"Pressure                    %15.6le\n",accumulator[18][5]);
    fprintf(fpouts,"   std. dev.                %15.4lf\n",accumulator[18][6]);
    fprintf(fpouts,"   fluctuation              %15.4lf\n",accumulator[18][7]);

    fprintf(fpouts,"Box volume                  %15.6le\n",accumulator[19][5]);
    fprintf(fpouts,"   std. dev.                %15.4lf\n",accumulator[19][6]);
    fprintf(fpouts,"   fluctuation              %15.4lf\n",accumulator[19][7]);

    fprintf(fpouts,"Ideal pressure              %15.6le\n",accumulator[20][5]);
    fprintf(fpouts,"   std. dev.                %15.4lf\n",accumulator[20][6]);
    fprintf(fpouts,"   fluctuation              %15.4lf\n",accumulator[20][7]);

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
    char keyword[100];
    int isDataFound;
    int atomid, moleid;
    int isFirstAtom;
    int last_mole_id;

    int position_counter;
    int starting_position;

    fprintf(stderr,"Reading input file...\n");
    fprintf(fpouts,"Reading input file...\n");

    /* read input file */
    fpins = fopen(INPUT,"r");

    sscanf(fgets(buffer,datalen,fpins), "%d %d", &ij, &jk);
    sscanf(fgets(buffer,datalen,fpins), "%lf", &treq);
    sscanf(fgets(buffer,datalen,fpins), "%lf %lf %lf", &boxlx, &boxly, &boxlz);
    sscanf(fgets(buffer,datalen,fpins), "%lf %lf %lf", &rcuton, &rcutoff, &rcutoffelec);
    // sscanf(fgets(buffer,datalen,fpins), "%s", coords_file);
    sscanf(fgets(buffer,datalen,fpins), "%d %d %d", &nstep, &fStart_option, &nstep_eq);
    sscanf(fgets(buffer,datalen,fpins), "%d %d %d %d %d", 
	    &nstep_ave, &nstep_print, &nstep_save, &nstep_ss, &nstep_trj);
    sscanf(fgets(buffer,datalen,fpins), "%lf %d", &delt, &nstep_inner);
    sscanf(fgets(buffer,datalen,fpins), "%lf", &f0);
    sscanf(fgets(buffer,datalen,fpins), "%d %d", &what_ensemble, &whichNH);
    sscanf(fgets(buffer,datalen,fpins), "%lf %lf", &qq, &qqs);
    sscanf(fgets(buffer,datalen,fpins), "%d", &isLJswitchOn);
    sscanf(fgets(buffer,datalen,fpins), "%d", &isChargeOn);
    // sscanf(fgets(buffer,datalen,fpins), "%d", &isWolfOn);
    // sscanf(fgets(buffer,datalen,fpins), "%d %d %d", &KMAXX, &KMAXY, &KMAXZ);
    // sscanf(fgets(buffer,datalen,fpins), "%d", &KSQMAX);
    // sscanf(fgets(buffer,datalen,fpins), "%lf", &kappa);

    sscanf(fgets(buffer,datalen,fpins), "%d", &nconstraint);

    sscanf(fgets(buffer,datalen,fpins), "%d", &isSFon);
    sscanf(fgets(buffer,datalen,fpins), "%d", &sf_type);


    // read in electrostatic parametes if needed
    int electype; // type of electrostatic
    if (isChargeOn)
    {
	// rewind the input file stream
	rewind(fpins);
	isDataFound = false;
       	while (fgets(buffer,datalen,fpins)!=NULL)
       	{
	    sscanf(buffer,"%s",keyword);
	    for (ii=0;ii<strlen(keyword);ii++)
	       	keyword[ii] = toupper(keyword[ii]);
	    if (!strcmp(keyword,"ELECTROSTATIC"))
	    {
	       	fprintf(stderr,"Data section for electrostatics found...\n");
	       	fprintf(fpouts,"Data section for electrostatics found...\n");
		isDataFound = true;
		// preset all electrostatic methods to false
		// then turn on them according to the input
		isEwaldOn = isWolfOn = isSimpleCoulomb = 0;
		sscanf(fgets(buffer,datalen,fpins),"%d",&electype);
		if (electype==elec_ewald)
		{
		    isEwaldOn = 1;
		    sscanf(fgets(buffer,datalen,fpins), "%lf", &kappa);
		    sscanf(fgets(buffer,datalen,fpins), "%d %d %d %d", &KMAXX, &KMAXY, &KMAXZ, &KSQMAX);
		    sscanf(fgets(buffer,datalen,fpins), "%d %d",&fEwald_BC,&fEwald_Dim);
		    fprintf(fpouts,"%dD Ewald summation requested...\n",fEwald_Dim);
		    fprintf(fpouts,"kappa=%lf\n",kappa);
		    fprintf(fpouts,"kmaxx=%d kmaxy=%d kmaxz=%d kmaxsq=%d\n",KMAXX, KMAXY, KMAXZ, KSQMAX);
		    fprintf(fpouts,"Boundary condition is %d\n",fEwald_BC);
		}
		else if (electype==elec_wolf)
		{
		    isWolfOn = 1;
		    sscanf(fgets(buffer,datalen,fpins), "%lf", &kappa);
		    fprintf(fpouts,"Wolf method requested...\n");
		    fprintf(fpouts,"kappa=%lf\n",kappa);
		}
		else if (electype==elec_simple_coulomb)
		{
		    isSimpleCoulomb = 1;
		    fprintf(fpouts,"Simple coulomb interaction requested...\n");
		} 
		break;
	    }
	} // read through lines
	if (isDataFound==false)
	{ 
	    fprintf(stderr,"Error: data for electrostatics not found.\n");
	    fprintf(fpouts,"Error: data for electrostatics not found.\n");
	    exit(1);
	}
    } // if charge needed


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
    fgets(buffer,datalen,fpcfg);
    sscanf(buffer,"%d atoms%n",&natom,&position_counter);
    // printf("natom=%d  counter=%d\n",natom,position_counter);
    // printf("%s",&buffer[position_counter]);
    // fscanf(fpcfg, "%d atoms%n", &natom,&position_counter);
    assert(natom<=natom_max);
    // read in the information of how many molecules in one specie and
    // how many atoms in one molecule for a certain specie
    specie_first_atom_idx[0] = 0;
    starting_position = position_counter;
    for (ii=0;ii<nspecie;ii++)
    {
	sscanf(&buffer[starting_position],"%d %d%n",&nmole_per_specie[ii],&natom_per_mole[ii],&position_counter);
	// advance the reading position in the string to next data section
	starting_position += position_counter;
	// printf("%d %d ",nmole_per_specie[ii],natom_per_mole[ii]);
	specie_first_atom_idx[ii+1] = specie_first_atom_idx[ii] + nmole_per_specie[ii]*natom_per_mole[ii];
    }
    // readin detailed atom information
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
    // fscanf(fpcfg, "%d bonds\n", &nbond);
    sscanf(fgets(buffer,datalen,fpcfg),"%d bonds%n",&nbond,&position_counter);
    assert(nbond<=nbond_max);
    // readin in how many bonds in a molecule of a certain specie
    starting_position = position_counter;
    for (ii=0;ii<nspecie;ii++)
    {
	sscanf(&buffer[starting_position],"%d%n",&nbond_per_mole[ii],&position_counter);
	starting_position += position_counter;
    }
    // readin detailed bond information
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
    // fscanf(fpcfg, "%d angles\n", &nangle);
    sscanf(fgets(buffer,datalen,fpcfg),"%d angles%n",&nangle,&position_counter);
    assert(nangle<=nangle_max);
    // readin in how many angles in a molecule of a certain specie
    starting_position = position_counter;
    for (ii=0;ii<nspecie;ii++)
    {
	sscanf(&buffer[starting_position],"%d%n",&nangle_per_mole[ii],&position_counter);
	starting_position += position_counter;
    }
    // read in detailed angle information
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
    // fscanf(fpcfg, "%d dihedrals\n", &ndih);
    sscanf(fgets(buffer,datalen,fpcfg),"%d dihedrals%n",&ndih,&position_counter);
    assert(ndih<=ndih_max);
    // read in how many angles in a molecule of a certain specie
    starting_position = position_counter;
    for (ii=0;ii<nspecie;ii++)
    {
	sscanf(&buffer[starting_position],"%d%n",&nangle_per_mole[ii],&position_counter);
	starting_position += position_counter;
    }
    // readin detailed dihedral information
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
    // fscanf(fpcfg, "%d impropers\n", &nimp);
    sscanf(fgets(buffer,datalen,fpcfg),"%d impropers%n",&nimp,&position_counter);
    assert(nimp<=nimp_max);
    // readin how many impropers in a molecule of a certain specie
    starting_position = position_counter;
    for (ii=0;ii<nspecie;ii++)
    {
	sscanf(&buffer[starting_position],"%d%n",&nimp_per_mole[ii],&position_counter);
	starting_position += position_counter;
    }

    // read in nonbonded pair list
    // fscanf(fpcfg, "%d nonbonded\n", &nnbp);
    sscanf(fgets(buffer,datalen,fpcfg),"%d nonbonded%n",&nnbp,&position_counter);
    assert(nnbp<=nnbp_max);
    // read in how many nonbonded pairs in a molecule of a certain specie
    starting_position = position_counter;
    for (ii=0;ii<nspecie;ii++)
    {
	sscanf(&buffer[starting_position],"%d%n",&nnbp_per_mole[ii],&position_counter);
	starting_position += position_counter;
    }
    // readin detailed nonbonded pair information
    last_mole_id = 0;
    for (ii=0;ii<nnbp;ii++)
    {
	fscanf(fpcfg,"%d %d\n",&nbp_idx[ii][0],&nbp_idx[ii][1]);
	// check if nonbonded atom is belong to the molecule
	if (mole_first_atom_idx[last_mole_id]<nbp_idx[ii][0] && nbp_idx[ii][0]< mole_first_atom_idx[last_mole_id+1])
	{
	    mole_first_nbp_idx[last_mole_id] = ii;
	    last_mole_id++;
	    assert(last_mole_id<nmole_max+1);
	}
    }

    fclose(fpcfg);
}

int echo()
{
    fprintf(fpouts,"==================================================\n");
    fprintf(fpouts,"Energies are for all molecules.\n");
    fprintf(fpouts,"Only total energy contains long range corrections.\n");
    fprintf(fpouts,"Pressure is with long range corrections.\n");
    fprintf(fpouts,"==================================================\n");
    fprintf(fpouts,"Units:\n");
    fprintf(fpouts,"  Energy(J/mol)\n  Temperature(K)\n  Velocity(m/s)\n  Distance(Angstrom)\n");
    fprintf(fpouts,"  Force(J/Angstrom/mol)\n");
    fprintf(fpouts,"  virial(J/mol)\n");
    fprintf(fpouts,"  Pressure(Pascal)\n");
    fprintf(fpouts,"==================================================\n");
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
    fprintf(fpouts,"preq=%lf\nQts=%lf\nQbs=%lf\n",preq,Qts,Qbs);
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
    fprintf(fpouts,"uer_vdw=%lf\n",uvdw-unbp_vdw);
    fprintf(fpouts,"unbp_vdw=%lf\n",unbp_vdw);
    fprintf(fpouts,"ubond=%lf\n",ubond);
    fprintf(fpouts,"uangle=%lf\n",uangle);
    fprintf(fpouts,"udih=%lf\n",udih);
    fprintf(fpouts,"uimp=%lf\n",uimp);
    fprintf(fpouts,"uewald=%lf\n",uewald);
    fprintf(fpouts,"ureal=%lf\n",ureal);
    fprintf(fpouts,"ufourier=%lf\n",ufourier);
    fprintf(fpouts,"uself=%lf\n",uself);
    fprintf(fpouts,"uexcl=%lf\n",uexcl);
    fprintf(fpouts,"uvacuum=%lf\n",uvacuum);
    fprintf(fpouts,"uGz0=%lf\n",uGz0);

    fprintf(fpouts,"usflj=%lf\n",usflj);
    fprintf(fpouts,"uwolf=%lf\n",uwolf);
    fprintf(fpouts,"uwolf_real=%lf\n",uwolf_real);
    fprintf(fpouts,"uwolf_con=%lf\n",uwolf_con);
    fprintf(fpouts,"ucoulomb=%lf\n",ucoulomb);

    fprintf(fpouts,"tinst=%lf\n",tinst);

    fprintf(fpouts,"virial=%lf\n",virial_inter+virial_intra);
    fprintf(fpouts,"virial_inter=%lf\n",virial_inter);
    fprintf(fpouts,"virial_intra=%lf\n",virial_intra);
    fprintf(fpouts,"pideal=%lf\n",pideal=natom/(boxlx*boxly*boxlz)*tinst*kb_1e30);
    // pljlrc already calculated in initialization part
    pinst = pideal
	+ (virial_inter+virial_intra)*virial_to_pressure/(boxlx*boxly*boxlz)
	+ pljlrc;
    fprintf(fpouts,"pressure=%lf\n",pinst);

    fprintf(fpouts,"uljlrc=%lf\npljlrc=%lf\n",uljlrc,pljlrc);

    fflush(fpouts);
}

int make_exclude_list()
{
    int ii, jj, nexcllist;

    fprintf(fpouts,"making exclude list... ");

    nexcllist = 0;
    for (ii=0;ii<natom;ii++)
    {
	pointexcl[ii] = nexcllist;
	// exclude bonded atoms
	// In order to improve it and save the memory usage
	// Only set atom jj to be an exclude atom for atom ii
	// if jj>ii. It should avoid double count of the exclulde pairs.
	// Since the calculations loop ij is from i=0,n-1 and j=i+1,n
	// i.e. if 0 and 1 are exclude pair, 1 will be in 0's exclude list
	// but 0 will not be in 1's exclude list
	// Be careful of it while coding other interaction functions
	for (jj=0;jj<nbond;jj++)
	{
	    if (bond_idx[jj][0] == ii && bond_idx[jj][1] > ii)
	    {
		excllist[nexcllist] = bond_idx[jj][1];
		nexcllist++;
	    }
	    else if (bond_idx[jj][1] == ii && bond_idx[jj][0] > ii)
	    {
		excllist[nexcllist] = bond_idx[jj][0];
		nexcllist++;
	    }
	    assert(nexcllist<exclude_max);
	}
	// exclude angled atoms
	for (jj=0;jj<nangle;jj++)
	{
	    if (angle_idx[jj][0] == ii && angle_idx[jj][2] > ii)
	    {
		excllist[nexcllist] = angle_idx[jj][2];
		nexcllist++;
	    }
	    else if (angle_idx[jj][2] == ii && angle_idx[jj][0] > ii)
	    {
		excllist[nexcllist] = angle_idx[jj][0];
		nexcllist++;
	    }
	    assert(nexcllist<exclude_max);
	}
	// exclude dihedraled atoms
	for (jj=0;jj<ndih;jj++)
	{
	    if (dih_idx[jj][0] == ii && dih_idx[jj][3] > ii)
	    {
		excllist[nexcllist] = dih_idx[jj][3];
		nexcllist++;
	    }
	    else if (dih_idx[jj][3] == ii && dih_idx[jj][0] > ii)
	    {
		excllist[nexcllist] = dih_idx[jj][0];
		nexcllist++;
	    }
	    assert(nexcllist<exclude_max);
	}
	// exclude non-bonded atoms
	for (jj=0;jj<nnbp;jj++)
	{
	    if (nbp_idx[jj][0] == ii && nbp_idx[jj][1] > ii)
	    {
		excllist[nexcllist] = nbp_idx[jj][1];
		nexcllist++;
	    }
	    else if (nbp_idx[jj][1] == ii && nbp_idx[jj][0] > ii)
	    {
		excllist[nexcllist] = nbp_idx[jj][0];
		nexcllist++;
	    }
	    assert(nexcllist<exclude_max);
	}
    }
    pointexcl[natom] = nexcllist;

    fprintf(fpouts,"%d excluding pairs\n",nexcllist);

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
    virial = virial_inter + virial_intra;
    // add energy of thermostat, if nose hoover is not used, they will just be zero
    utot = utot + unhts + unhtss + utsbs;
    // add long range corrections into total energy
    utot = utot + uljlrc;
    // calculate ideal pressure part
    pideal=natom/(boxlx*boxly*boxlz)*tinst*kb_1e30;
    // do not need to recalculate lrc here, it should be calculated
    // elsewhere when variables changed
    pinst = pideal
	+ (virial_inter+virial_intra)*virial_to_pressure/(boxlx*boxly*boxlz)
	+ pljlrc;
    fprintf(stderr,"%10d %10.4le %10.4le %10.4le %10.4le %10.4le %10.4le\n",
	    istep,utot,upot,ukin,tinst,pinst,boxv);

    fprintf(fplog,"%10d %10.4le %10.4le %10.4le %10.4le %10.4le %10.4le %10.4le %10.4le %10.4le ",
	    istep,utot,upot,ukin,tinst,uinter,uintra,uvdw,ubond,uangle);

    fprintf(fplog,"%10.4le %10.4le %10.4le %10.4le %10.4le %10.4le %10.4le %10.4le %10.4le ",
	    uangle,udih,uimp,uewald,usflj,unhts,unhtss,virial,virial_inter,virial_intra,utsbs);

    fprintf(fplog,"%10.4le %10.4le %10.4le\n",
	    utsbs,pinst,boxv);
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

    fwrite(&vts,sizeof(double),1,fpsave);
    fwrite(&rts,sizeof(double),1,fpsave);
    fwrite(&vbs,sizeof(double),1,fpsave);

    fwrite(&boxlx,sizeof(double),1,fpsave);
    fwrite(&boxly,sizeof(double),1,fpsave);
    fwrite(&boxlz,sizeof(double),1,fpsave);
    fwrite(&boxv,sizeof(double),1,fpsave);

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

    // recaculate kinetic energy and temperature using the loaded velocities
    ukin = 0.0;
    for (ii=0;ii<natom;ii++)
	ukin += aw[ii]*(vx[ii]*vx[ii]+vy[ii]*vy[ii]+vz[ii]*vz[ii]);
    ukin = 0.5*ukin;
    tinst = 2.0*ukin/(Rgas*nfree);

    // calculate instant temperature
    fread(&qq,sizeof(double),1,fpload);
    fread(&ps,sizeof(double),1,fpload);
    fread(&gg,sizeof(double),1,fpload);
    fread(&ss,sizeof(double),1,fpload);
    fread(&qqs,sizeof(double),1,fpload);
    fread(&pss,sizeof(double),1,fpload);
    fread(&ggs,sizeof(double),1,fpload);
    fread(&sss,sizeof(double),1,fpload);

    fread(&vts,sizeof(double),1,fpload);
    fread(&rts,sizeof(double),1,fpload);
    fread(&vbs,sizeof(double),1,fpload);

    fread(&boxlx,sizeof(double),1,fpload);
    fread(&boxly,sizeof(double),1,fpload);
    fread(&boxlz,sizeof(double),1,fpload);
    fread(&boxv,sizeof(double),1,fpload);

    // recalculate the thermostat energy
    utsbs = 0.5*Qbs*vbs*vbs + 0.5*Qts*vts*vts + (nfree+1)*Rgas*treq*rts + preq*boxv*PascalA3_to_J_mol;
    // recalculate the long range corrections since the box size may be changed
    calculate_ljlrc();

    // read counters and accumulators only if its a continue run
    if (fStart_option==continue_run)
    {
       	fread(&istep,sizeof(int),1,fpload);
	nstep_start = istep + 1;
       	fread(icounter,sizeof(int),num_counter_max,fpload);
       	fread(accumulator,sizeof(double),num_counter_max,fpload);
    }

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
	if (what_ensemble == nvt_run) // velocity verlet with nose hoover for NVT
	{
	    vver_nh_3();
	    /*
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
	    */
	}
	else if (what_ensemble == npt_run)
	{
	    npt_respa();
	}
	else 
	    vver(); // velocity verlet

	// print out, snapshot, trajectory, save
	if (istep%nstep_print == 0) printit();

	if (nstep_ss && istep%nstep_ss == 0) snapshot();

	if (nstep_trj && istep%nstep_trj==0) trajectory();

	if (istep%nstep_save==0) saveit();

	icounter[11]--;
	// if still in equilibrium run
	// do not do averages
	if (icounter[11]>=0) continue;

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
	accumulator[16][0] += uvacuum;
	accumulator[16][1] += uvacuum*uvacuum;
	accumulator[17][0] += uwolf;
	accumulator[17][1] += uwolf*uwolf;
	accumulator[18][0] += pinst;
	accumulator[18][1] += pinst*pinst;
	accumulator[19][0] += boxv;
	accumulator[19][1] += boxv*boxv;
	accumulator[20][0] += pideal;
	accumulator[20][1] += pideal*pideal;

	if ((istep-nstep_eq)%nstep_ave==0) averages();

    } // loop through all simulation time

    snapshot();

    calres();

    fprintf(stderr,"%d frames in the trajectory file.\n",nframe);
    fprintf(fpouts,"%d frames in the trajectory file.\n",nframe);

    // write out results and clean up
    ending();

}

