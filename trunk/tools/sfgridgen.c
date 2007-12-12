/*
 * Generate grids of Solid-Fluid potentials and forces.
 * Used for interpolations
 * Calculate the potentials, forces (1st derivatives),
 * 2nd derivatives and save them as binary format
 *
 * Input and output have different units
 * Input: energy is K
 * Output: energy is J/mol. force is J/mol/Angstrom
 *
 * To compile: gcc sfgridgen.c -o sfgridgen.x -lm
 *
 * Currently only for nanotubes which are peridocial in z direction
 *
 * Written by Yang Wang 2007
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

// maximum number of grids in x,y,z directions
#define max_grids	999999
#define solid_atom_max	1000   //number of atoms in one solid (absorbent) unit cell 
#define Rgas		8.314472 // J/mol/K 

// f(x,y,z) potentials
double Ene[max_grids];
// f'(x,y,z) 1st derivatives, forces
double dEx[max_grids], dEy[max_grids], dEz[max_grids];
// f''(x,y,z) 2nd derivatives
double dFxx[max_grids], dFxy[max_grids], dFxz[max_grids];
double dFyx[max_grids], dFyy[max_grids], dFyz[max_grids];
double dFzx[max_grids], dFzy[max_grids], dFzz[max_grids];
// x,y,z coordinates of solid atoms
double sxx[solid_atom_max];
double syy[solid_atom_max];
double szz[solid_atom_max];
// x,y,z of fluid atom
double fxx, fyy, fzz;

int s_natom; // number of atoms in solid
double uclx, ucly, uclz; // the length in x,y,z directions of a unit cell
double grid_itvl; // grid interval
int ngrid_x, ngrid_y, ngrid_z; // number of grids in x,y,z directions
int ngrid_total;
double rcuton, rcutoff; // cuton and cutoff
double rcutonsq, rcutoffsq;
double roff2_minus_ron2_cube;
int isLJswitch; // switch energy for LJ interactions
double factorLJswitch;
int nuc_z; // number of unit cell in z direction

char f_atomname[100]; // name of the fluid atom
char s_atomname[100]; // solid atom name
double f_sigma, f_epsilon; // fluid sigma, epsilon
double s_sigma, s_epsilon; // solid sigma, epsilon
double sigmasf, epsilonsf; // mixing terms

FILE *fpins; // input file
FILE *fpgridfile; // file pointer of output grid file

// lower and higher bound of x,y,z
double xcenter, ycenter, zcenter;
double xmin, xmax;
double ymin, ymax;
double zmin, zmax;

int main (int argc, char *argv[])
{
    int ii, jj, kk;
    int pp, nn;
    double rxsf, rysf, rzsf;
    double rsfsq, rsf;
    double r_rsfsq;
    double r_r6, r_r12, r_r12_minus_r_r6;
    double usf_vdw_temp;
    double usf_vdw;
    double fsf, fxsf, fysf, fzsf; // 1st
    double Dusfx, Dusfy, Dusfz;

    double ddf; // 2nd 
    double dfxx, dfxy, dfxz;
    double dfyx, dfyy, dfyz;
    double dfzx, dfzy, dfzz;

    // accumulator
    double Dfxx, Dfxy, Dfxz;
    double Dfyx, Dfyy, Dfyz;
    double Dfzx, Dfzy, Dfzz;

    double temp1, temp2, temp3, temp4;
    double temp1x, temp1y, temp1z;
    double temp2x, temp2y, temp2z;
    double temp3x, temp3y, temp3z;
    double temp4x, temp4y, temp4z;
    double AA;
    int index;

    // open input file
    fpins = fopen(argv[1],"r");

    // read in input data
    // fluid atom properties, name, sigma, epsilon (Angstrom and K)
    fscanf(fpins,"%s %lf %lf\n",f_atomname,&f_sigma,&f_epsilon);
    // cut off and if use LJ switch potential
    fscanf(fpins,"%lf %lf %lf %d\n",&grid_itvl,&rcuton,&rcutoff,&isLJswitch);
    // solid properties, number of atoms, length in x,y,z directions
    fscanf(fpins,"%d %lf %lf %lf\n",&s_natom,&uclx,&ucly,&uclz);
    // more solid properties, sigma and epsilon (Angstrom and K)
    fscanf(fpins,"%lf %lf\n",&s_sigma,&s_epsilon);
    // read in solid coordinates
    for (ii=0;ii<s_natom;ii++)
	fscanf(fpins,"%s %lf %lf %lf\n",s_atomname,&sxx[ii],&syy[ii],&szz[ii]);

    fclose(fpins);

    // initialize
    rcutonsq = rcuton*rcuton;
    rcutoffsq = rcutoff*rcutoff;
    roff2_minus_ron2_cube = (rcutoffsq-rcutonsq)*(rcutoffsq-rcutonsq)*(rcutoffsq-rcutonsq);

    sigmasf = (f_sigma + s_sigma)/2.0;
    epsilonsf = sqrt(f_epsilon*s_epsilon)*Rgas; // mixing and change it to J/mol

    // number of grids in x,y,z directions
    ngrid_x = (int)floor(uclx/grid_itvl)+2;
    ngrid_y = (int)floor(ucly/grid_itvl)+2;
    ngrid_z = (int)floor(uclz/grid_itvl)+2;
    printf("number of grids in x,y,z directions: %d %d %d\n",ngrid_x,ngrid_y,ngrid_z);
    ngrid_total = ngrid_x*ngrid_y*ngrid_z;
    printf("number of total grids = %d\n",ngrid_total);
    assert(ngrid_total<max_grids);

    // number of repeat unit cells in z direction based on the cut off
    nuc_z = 1;
    while (nuc_z*uclz<rcutoff)
	nuc_z++;
    printf("number of unit cells in z direction: %d\n",nuc_z);

    // set up the boundary
    xcenter = uclx/2.0;
    ycenter = ucly/2.0;
    zcenter = 0.0;
    xmax = xcenter + uclx/2.0;
    xmin = xcenter - uclx/2.0;
    ymax = ycenter + ucly/2.0;
    ymin = ycenter - ucly/2.0;
    zmin = 0.0;
    zmax = uclz;
    printf("X goes from %lf to %lf\n",xmin,xmax);
    printf("Y goes from %lf to %lf\n",ymin,ymax);
    printf("Z goes from %lf to %lf\n",zmin,zmax);

    for (ii=0;ii<ngrid_x;ii++)
    {
	fxx = xmin + grid_itvl*ii;
	for (jj=0;jj<ngrid_y;jj++)
	{
	    fyy = ymin + grid_itvl*jj;
	    for (kk=0;kk<ngrid_z;kk++)
	    {
		fzz = zmin + grid_itvl*kk;

		// reset s-f energy for particle on position of ii,jj,kk
		usf_vdw = 0.0;
		Dusfx = 0.0;
		Dusfy = 0.0;
		Dusfz = 0.0;
		Dfxx = Dfxy = Dfxz = 0.0;
		Dfyx = Dfyy = Dfyz = 0.0;
		Dfzx = Dfzy = Dfzz = 0.0;
		double xx,yy,zz;
		// xx=8.959794; yy=5.002037; zz=0.257057;
		// xx=9.045419; yy=3.828860; zz=1.608055;
		// xx=4.212336; yy=6.596563; zz=2.424706;
		// xx=5.184985; yy=4.861523; zz=2.065608;
		// xx=9.915239; yy=10.662925; zz=0.076776;
		// xx=4.284484; yy=5.045358; zz=3.193456;
		// xx=4.455993; yy=4.990464; zz=0.727417;
		// xx=9.915239; yy=10.662925; zz=2.535997;
		// fxx=xx;  fyy=yy;  fzz=zz;
		for (pp=-nuc_z;pp<=nuc_z;pp++)
		{
		    for (nn=0;nn<s_natom;nn++)
		    {
			rxsf = fxx - sxx[nn];
			rysf = fyy - syy[nn];
			rzsf = fzz - szz[nn];
			// add peridoic z length
			rzsf += pp*uclz;
			rsfsq = rxsf*rxsf + rysf*rysf + rzsf*rzsf;
			if (rsfsq<rcutoffsq)
			{
			    rsf = sqrt(rsfsq);
			    if (isLJswitch)
			    {
				if (rsfsq<rcutonsq)
				    factorLJswitch = 1.0;
				else
				    factorLJswitch = (rcutoffsq-rsfsq)*(rcutoffsq-rsfsq)
					*(rcutoffsq+2.0*rsfsq-3.0*rcutonsq)/roff2_minus_ron2_cube;
			    } // if LJ switch is used
			    r_rsfsq = sigmasf*sigmasf/rsfsq;
			    r_r6 = r_rsfsq*r_rsfsq*r_rsfsq;
			    r_r12 = r_r6*r_r6;
			    r_r12_minus_r_r6 = r_r12 - r_r6;
			    usf_vdw_temp = epsilonsf*r_r12_minus_r_r6; // still need 4.0
			    if (isLJswitch)
				usf_vdw += usf_vdw_temp*factorLJswitch; // still need 4.0
			    else
				usf_vdw += usf_vdw_temp; // still need 4.0

			    // forces (1st derivatives)
			    if (isLJswitch)
			    {
				if (rsfsq<rcutonsq)
				    fsf = 24.0*epsilonsf*(r_r12_minus_r_r6+r_r12)/rsfsq;
				else
				    fsf = 24.0*epsilonsf*(r_r12_minus_r_r6+r_r12)/rsfsq*factorLJswitch
					-4.0*usf_vdw_temp*12.0*(rcutoffsq-rsfsq)*(rcutonsq-rsfsq)/roff2_minus_ron2_cube;
			    }
			    else
				fsf = 24.0*epsilonsf*(r_r12_minus_r_r6+r_r12)/rsfsq;
			    fxsf = fsf*rxsf;
			    fysf = fsf*rysf;
			    fzsf = fsf*rzsf;
			    // the derivative is the opposite of the force since force=-dU/dr
			    Dusfx -= fxsf;
			    Dusfy -= fysf;
			    Dusfz -= fzsf;

			    // 2nd derivatives
			    if (isLJswitch)
			    {
				if (rsfsq<rcutonsq)
				{
				    ddf = 96.0*epsilonsf*(r_r12_minus_r_r6+6.0*r_r12-r_r6)/rsfsq/rsfsq;
				    // dfx 2nd 
				    dfxx = -rxsf*ddf*rxsf + fsf;
				    dfxy = -rxsf*ddf*rysf;
				    dfxz = -rxsf*ddf*rzsf;
				    // dfy 2nd
				    dfyx = -rysf*ddf*rxsf;
				    dfyy = -rysf*ddf*rysf + fsf;
				    dfyz = -rysf*ddf*rzsf;
				    // dfz 2nd
				    dfzx = -rzsf*ddf*rxsf;
				    dfzy = -rzsf*ddf*rysf;
				    dfzz = -rzsf*ddf*rzsf + fsf;
				}
				else
				{
				}
			    }
			    else // no switch 
			    {
				ddf = 96.0*epsilonsf*(r_r12_minus_r_r6+6.0*r_r12-r_r6)/rsfsq/rsfsq;
				// dfx 2nd 
				dfxx = -rxsf*ddf*rxsf + fsf;
				dfxy = -rxsf*ddf*rysf;
				dfxz = -rxsf*ddf*rzsf;
				// dfy 2nd
				dfyx = -rysf*ddf*rxsf;
				dfyy = -rysf*ddf*rysf + fsf;
				dfyz = -rysf*ddf*rzsf;
				// dfz 2nd
				dfzx = -rzsf*ddf*rxsf;
				dfzy = -rzsf*ddf*rysf;
				dfzz = -rzsf*ddf*rzsf + fsf;

			    }
			    Dfxx += dfxx;
			    Dfxy += dfxy;
			    Dfxz += dfxz;
			    Dfyx += dfyx;
			    Dfyy += dfyy;
			    Dfyz += dfyz;
			    Dfzx += dfzx;
			    Dfzy += dfzy;
			    Dfzz += dfzz;
			} // cut off check
		    } // nn loop, solid atoms
		} // pp loop, periodical in z direction
		usf_vdw *= 4.0;
		// calculate index
		index = ii*ngrid_y*ngrid_z + jj*ngrid_z + kk;
		Ene[index] = usf_vdw;
		// 1st 
		dEx[index] = Dusfx;
		dEy[index] = Dusfy;
		dEz[index] = Dusfz;
		// 2nd
		dFxx[index] = Dfxx;
		dFxy[index] = Dfxy;
		dFxz[index] = Dfxz;
		dFyx[index] = Dfyx;
		dFyy[index] = Dfyy;
		dFyz[index] = Dfyz;
		dFzx[index] = Dfzx;
		dFzy[index] = Dfzy;
		dFzz[index] = Dfzz;

		/*
		   printf("fxx=%lf  fyy=%lf   fzz=%lf\n",fxx,fyy,fzz);
		   printf("sigmasf=%lf epsilonsf=%lf\n",sigmasf,epsilonsf);
		   printf("usf_vdw=%lf\n",usf_vdw);
		   printf("Udx=%lf Udy=%lf Udz=%lf\n",Dusfx, Dusfy, Dusfz);
		   printf("Dfxx=%lf Dfxy=%lf Dfxz=%lf\n",Dfxx,Dfxy,Dfxz);
		   printf("Dfyx=%lf Dfyy=%lf Dfyz=%lf\n",Dfyx,Dfyy,Dfyz);
		   printf("Dfzx=%lf Dfzy=%lf Dfzz=%lf\n",Dfzx,Dfzy,Dfzz);
		   exit(1);
		   */

	    } // kk loop
	} // jj loop
    } // ii loop


    // output the grid data
    fpgridfile = fopen("out.sfgridgen","wb");
    fwrite(&uclx,sizeof(double),1,fpgridfile);
    fwrite(&ucly,sizeof(double),1,fpgridfile);
    fwrite(&uclz,sizeof(double),1,fpgridfile);
    fwrite(&xcenter,sizeof(double),1,fpgridfile);
    fwrite(&ycenter,sizeof(double),1,fpgridfile);
    fwrite(&zcenter,sizeof(double),1,fpgridfile);
    fwrite(&ngrid_total,sizeof(int),1,fpgridfile);
    fwrite(&ngrid_x,sizeof(int),1,fpgridfile);
    fwrite(&ngrid_y,sizeof(int),1,fpgridfile);
    fwrite(&ngrid_z,sizeof(int),1,fpgridfile);
    fwrite(&grid_itvl,sizeof(double),1,fpgridfile);
    fwrite(&grid_itvl,sizeof(double),1,fpgridfile);
    fwrite(&grid_itvl,sizeof(double),1,fpgridfile);
    fwrite(Ene,sizeof(double),ngrid_total,fpgridfile);
    fwrite(dEx,sizeof(double),ngrid_total,fpgridfile);
    fwrite(dEy,sizeof(double),ngrid_total,fpgridfile);
    fwrite(dEz,sizeof(double),ngrid_total,fpgridfile);
    fwrite(dFxx,sizeof(double),ngrid_total,fpgridfile);
    fwrite(dFxy,sizeof(double),ngrid_total,fpgridfile);
    fwrite(dFxz,sizeof(double),ngrid_total,fpgridfile); 
    fwrite(dFyx,sizeof(double),ngrid_total,fpgridfile);
    fwrite(dFyy,sizeof(double),ngrid_total,fpgridfile);
    fwrite(dFyz,sizeof(double),ngrid_total,fpgridfile); 
    fwrite(dFzx,sizeof(double),ngrid_total,fpgridfile);
    fwrite(dFzy,sizeof(double),ngrid_total,fpgridfile);
    fwrite(dFzz,sizeof(double),ngrid_total,fpgridfile); 

    fclose(fpgridfile);
}


