#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <ctype.h>
#include <gsl/gsl_linalg.h>
#include "vars.h"

extern void cforce_atom_(int*, double*, double*, double*, double*, double*);
extern void initpotentialgrid_(int*, double*, double*, double*, int*, int*, int*, double*);
extern void pass_grid_file_name_(char*, int*);
extern void read_grids_(int* nspecies_yang);

double *amatrix, *amatrix_backup;
double *interp_vector;
double *bvector;
gsl_vector *work; // = gsl_vector_alloc(32);
gsl_vector *Svector; //  = gsl_vector_alloc(32);
gsl_vector *xvector; // = gsl_vector_alloc(32);
gsl_matrix *Vmatrix; // = gsl_matrix_alloc(32,32);
gsl_matrix *Xmatrix;

// initialize and reading the tasos grids
int init_tasos_grid()
{
    int	ii, jj;
    const int datalen = 200;
    char buffer[200];
    double auc, buc, cuc, nanotuberadius;
    int	na, nb, nc;
    int	nspecies_yang;
    int tnuatoms;
    char szgrid[12][200]; // assuming only 12 types of grids
    char keyword[100];

    fprintf(stderr,"Reading input data for tasos grids...\n");
    fprintf(fpouts,"Reading input data for tasos grids...\n");

    // re-open input file to read extra data section
    fpins = fopen(INPUT,"r");

    while (fgets(buffer,datalen,fpins)!=NULL)
    {
	sscanf(buffer,"%s",keyword);
	for (ii=0;ii<strlen(keyword);ii++)
	    keyword[ii] = toupper(keyword[ii]);
	if (!strcmp(keyword,"TASOS"))
	{
	    fprintf(stderr,"Data section for tasos grids found...\n");
	    fprintf(fpouts,"Data section for tasos grids found...\n");

	    sscanf(fgets(buffer,datalen,fpins), "%lf %lf %lf %lf", &auc, &buc, &cuc, &nanotuberadius);
	    sscanf(fgets(buffer,datalen,fpins), "%d %d %d", &na, &nb, &nc);
	    sscanf(fgets(buffer,datalen,fpins), "%d", &nspecies_yang);
	    sscanf(fgets(buffer,datalen,fpins), "%d", &tnuatoms);
	    for (ii=0;ii<tnuatoms;ii++)
		sscanf(fgets(buffer,datalen,fpins), "%s", szgrid+ii);

	    initpotentialgrid_(&tnuatoms, &auc, &buc, &cuc, &na, &nb, &nc, &nanotuberadius);
	    fprintf(stderr,"init potential ok\n");
	    fprintf(fpouts,"init potential ok\n");
	    for (ii=0;ii<tnuatoms;ii++)
	    {
		jj = ii + 1;
		pass_grid_file_name_(szgrid[ii], &jj);
		fprintf(stderr,"%s\n",szgrid[ii]);
	    }
	    fprintf(stderr,"pass grid file name ok\n");
	    fprintf(fpouts,"pass grid file name ok\n");
	    read_grids_(&nspecies_yang);
	    fprintf(stderr,"read grids ok\n");
	    fprintf(fpouts,"read grids ok\n");

	    fclose(fpins);
	    return 0;
	} // if keyword found
    } // read through the lines
    fprintf(stderr,"Error: data for tasos grids not found.\n");
    fprintf(fpouts,"Error: data for tasos grids not found.\n");
    fclose(fpins);
    exit(1);
}


int init_my_interp()
{
    int ii;
    FILE *fpgridfile;
    const int datalen = 200;
    char buffer[200];
    int nunique_atom;
    char szgrid[200];
    char keyword[100];
    double uclx_chk, ucly_chk, uclz_chk;
    double temp;

    fprintf(stderr,"Warning: this solid-fluid interpolation method is not fully tested.\n");
    fprintf(fpouts,"Warning: this solid-fluid interpolation method is not fully tested.\n");
    fprintf(stderr,"Warning: the algorithm is not efficient.\n");
    fprintf(fpouts,"Warning: the algorithm is not efficient.\n");
    fprintf(stderr,"Reading input data for myinterp...\n");
    fprintf(fpouts,"Reading input data for myinterp...\n");

    // re-open input file to read extra data section
    fpins = fopen(INPUT,"r");

    // change the atom type from tasos type to myinterp type
    // myinterp type = tasos type - 1
    for (ii=0;ii<natom;ii++)
	tasostype[ii] -= 1;

    // set up the variables needed for solving matrix
    work = gsl_vector_alloc(32);
    Svector = gsl_vector_alloc(32);
    xvector = gsl_vector_alloc(32);
    Vmatrix = gsl_matrix_alloc(32,32);
    Xmatrix = gsl_matrix_alloc(32,32);
    amatrix = calloc(1024,sizeof(double));
    amatrix_backup = calloc(1024,sizeof(double));
    interp_vector = calloc(32,sizeof(double));
    bvector = calloc(32,sizeof(double));

    // initiate pointers for s-f grids
    for (ii=0;ii<nunique_atom_max;ii++)
    {
	Ene[ii] = NULL;
	dEx[ii] = dEy[ii] = dEz[ii] = NULL;
	dFxx[ii] = dFxy[ii] = dFxz[ii] = NULL;
	dFyx[ii] = dFyy[ii] = dFyz[ii] = NULL;
	dFzx[ii] = dFzy[ii] = dFzz[ii] = NULL;
    }

    while (fgets(buffer,datalen,fpins)!=NULL)
    { 
	sscanf(buffer,"%s",keyword);
	for (ii=0;ii<strlen(keyword);ii++)
	    keyword[ii] = toupper(keyword[ii]);
	if (!strcmp(keyword,"MYINTERP"))
	{
	    fprintf(stderr,"Data section for myinterp found...\n");
	    fprintf(fpouts,"Data section for myinterp found...\n");

	    sscanf(fgets(buffer,datalen,fpins),"%lf %lf %lf %lf %lf %lf",
		    &uclx,&ucly,&uclz,&xcenter,&ycenter,&zcenter);
	    xmax = xcenter + uclx/2.0; 
	    xmin = xcenter - uclx/2.0;
	    ymax = ycenter + ucly/2.0; 
	    ymin = ycenter - ucly/2.0;
	    // note: zmin and zmax are different from x,y
	    zmin = 0.0; 
	    zmax = uclz; 
	    fprintf(stderr,"xmin=%lf  xmax=%lf\n",xmin,xmax);
	    fprintf(stderr,"ymin=%lf  ymax=%lf\n",ymin,ymax);
	    fprintf(stderr,"zmin=%lf  zmax=%lf\n",zmin,zmax);
	    fprintf(fpouts,"xmin=%lf  xmax=%lf\n",xmin,xmax);
	    fprintf(fpouts,"ymin=%lf  ymax=%lf\n",ymin,ymax);
	    fprintf(fpouts,"zmin=%lf  zmax=%lf\n",zmin,zmax);
	    sscanf(fgets(buffer,datalen,fpins),"%d",&nunique_atom);
	    for (ii=0;ii<nunique_atom;ii++)
	    {
		sscanf(fgets(buffer,datalen,fpins),"%s",szgrid);
		fpgridfile = fopen(szgrid,"rb");
		fread(&uclx_chk,sizeof(double),1,fpgridfile);
		fread(&ucly_chk,sizeof(double),1,fpgridfile);
		fread(&uclz_chk,sizeof(double),1,fpgridfile);
		if (uclx==uclx_chk && ucly==ucly_chk && uclz==uclz_chk)
		{
		    fprintf(stderr,"unit cell size matches for grid: %s\n",szgrid);
		    fprintf(fpouts,"unit cell size matches for grid: %s\n",szgrid);
		}
		else
		{
		    fprintf(stderr,
			    "Error: unit cell mismatch: %s\nuclx=%lf/%lf  ucly=%lf/%lf  uclz=%lf/%lf\n",
			    szgrid, uclx,uclx_chk,ucly,ucly_chk,uclz,uclz_chk);
		    fprintf(fpouts,
			    "Error: unit cell mismatch: %s\nuclx=%lf/%lf  ucly=%lf/%lf  uclz=%lf/%lf\n",
			    szgrid, uclx,uclx_chk,ucly,ucly_chk,uclz,uclz_chk);
		    exit(1);
		}
		// read in the center coordinates
		// not really useful here
		fread(&temp,sizeof(double),1,fpgridfile);
		fread(&temp,sizeof(double),1,fpgridfile);
		fread(&temp,sizeof(double),1,fpgridfile);
		fread(&ngrid_total,sizeof(int),1,fpgridfile);
		fread(&ngrid_x,sizeof(int),1,fpgridfile);
		fread(&ngrid_y,sizeof(int),1,fpgridfile);
		fread(&ngrid_z,sizeof(int),1,fpgridfile);
		fprintf(stderr,"%d=%d*%d*%d grids in grid file %d.\n",
			ngrid_total,ngrid_x,ngrid_y,ngrid_z,ii);
		fprintf(fpouts,"%d=%d*%d*%d grids in grid file %d.\n",
			ngrid_total,ngrid_x,ngrid_y,ngrid_z,ii);
		fread(&grid_itvl_x,sizeof(double),1,fpgridfile);
		fread(&grid_itvl_y,sizeof(double),1,fpgridfile);
		fread(&grid_itvl_z,sizeof(double),1,fpgridfile);
		fprintf(stderr,"grid size x,y,z = %lf,%lf,%lf\n",grid_itvl_x,grid_itvl_y,grid_itvl_z);
		fprintf(fpouts,"grid size x,y,z = %lf,%lf,%lf\n",grid_itvl_x,grid_itvl_y,grid_itvl_z);

		Ene[ii] = calloc(ngrid_total,sizeof(double));
		dEx[ii] = calloc(ngrid_total,sizeof(double));
		dEy[ii] = calloc(ngrid_total,sizeof(double));
		dEz[ii] = calloc(ngrid_total,sizeof(double));
		dFxx[ii] = calloc(ngrid_total,sizeof(double));
		dFxy[ii] = calloc(ngrid_total,sizeof(double));
		dFxz[ii] = calloc(ngrid_total,sizeof(double));
		dFyx[ii] = calloc(ngrid_total,sizeof(double));
		dFyy[ii] = calloc(ngrid_total,sizeof(double));
		dFyz[ii] = calloc(ngrid_total,sizeof(double));
		dFzx[ii] = calloc(ngrid_total,sizeof(double));
		dFzy[ii] = calloc(ngrid_total,sizeof(double));
		dFzz[ii] = calloc(ngrid_total,sizeof(double));

		assert(Ene[ii]!=NULL);
		assert(dEx[ii]!=NULL);
		assert(dEy[ii]!=NULL);
		assert(dEz[ii]!=NULL);
		assert(dFxx[ii]!=NULL);
		assert(dFxy[ii]!=NULL);
		assert(dFxz[ii]!=NULL);
		assert(dFyx[ii]!=NULL);
		assert(dFyy[ii]!=NULL);
		assert(dFyz[ii]!=NULL);
		assert(dFzx[ii]!=NULL);
		assert(dFzy[ii]!=NULL);
		assert(dFzz[ii]!=NULL);

		fread(Ene[ii],sizeof(double),ngrid_total,fpgridfile);
		fread(dEx[ii],sizeof(double),ngrid_total,fpgridfile);
		fread(dEy[ii],sizeof(double),ngrid_total,fpgridfile);
		fread(dEz[ii],sizeof(double),ngrid_total,fpgridfile);
		fread(dFxx[ii],sizeof(double),ngrid_total,fpgridfile);
		fread(dFxy[ii],sizeof(double),ngrid_total,fpgridfile);
		fread(dFxz[ii],sizeof(double),ngrid_total,fpgridfile);
		fread(dFyx[ii],sizeof(double),ngrid_total,fpgridfile);
		fread(dFyy[ii],sizeof(double),ngrid_total,fpgridfile);
		fread(dFyz[ii],sizeof(double),ngrid_total,fpgridfile);
		fread(dFzx[ii],sizeof(double),ngrid_total,fpgridfile);
		fread(dFzy[ii],sizeof(double),ngrid_total,fpgridfile);
		fread(dFzz[ii],sizeof(double),ngrid_total,fpgridfile);
		fclose(fpgridfile); 
	    } // loop through unique atoms grid files
	    fclose(fpins);
	    return 0;
	} // if key word found
    } // read through lines
    fprintf(stderr,"Error: data for myinterp not found.\n");
    fprintf(fpouts,"Error: data for myinterp not found.\n");
    fclose(fpins);
    exit(1);
}

int end_my_interp()
{
    gsl_vector_free(work);
    gsl_vector_free(Svector);
    gsl_vector_free(xvector);
    gsl_matrix_free(Vmatrix);
    gsl_matrix_free(Xmatrix);
    free(amatrix);
    free(amatrix_backup);
    free(interp_vector);
    free(bvector);
}

int Amatrix_ele_assign_value(double *line, double x, double y, double z)
{
    double x2, y2, z2;
    double x3, y3, z3;

    x2=x*x; x3=x2*x;
    y2=y*y; y3=y2*y;
    z2=z*z; z3=z2*z;

    // 32 elements in a line
    line[0] = 1.0;
    line[1] = x;
    line[2] = y;
    line[3] = z;
    line[4] = x2;
    line[5] = x*y;
    line[6] = x*z;
    line[7] = y2;
    line[8] = y*z;
    line[9] = z2;
    line[10] = x3;
    line[11] = x2*y;
    line[12] = x*y2;
    line[13] = x2*z;
    line[14] = x*z2;
    line[15] = y3;
    line[16] = y2*z;
    line[17] = y*z2;
    line[18] = z3;
    line[19] = x*y*z;
    line[20] = x3*y;
    line[21] = x2*y2;
    line[22] = x*y3;
    line[23] = x3*z;
    line[24] = x2*z2;
    line[25] = x*z3;
    line[26] = y3*z;
    line[27] = y2*z2;
    line[28] = y*z3;
    line[29] = x2*y*z;
    line[30] = x*y2*z;
    line[31] = x*y*z2;
}

int Amatrix_ele_assign_dx(double *line, double x, double y, double z)
{
    double x2, y2, z2;
    double x3, y3, z3;

    x2=x*x; x3=x2*x;
    y2=y*y; y3=y2*y;
    z2=z*z; z3=z2*z;

    // 32 elements in a line
    line[0] = 0;
    line[1] = 1;
    line[2] = 0;
    line[3] = 0;
    line[4] = 2*x;
    line[5] = y;
    line[6] = z;
    line[7] = 0;
    line[8] = 0;
    line[9] = 0;
    line[10] = 3*x2;
    line[11] = 2*x*y;
    line[12] = y2;
    line[13] = 2*x*z;
    line[14] = z2;
    line[15] = 0;
    line[16] = 0;
    line[17] = 0;
    line[18] = 0;
    line[19] = y*z;
    line[20] = 3*x2*y;
    line[21] = 2*x*y2;
    line[22] = y3;
    line[23] = 3*x2*z;
    line[24] = 2*x*z2;
    line[25] = z3;
    line[26] = 0;
    line[27] = 0;
    line[28] = 0;
    line[29] = 2*x*y*z;
    line[30] = y2*z;
    line[31] = y*z2;
}


int Amatrix_ele_assign_dy(double *line, double x, double y, double z)
{
    double x2, y2, z2;
    double x3, y3, z3;

    x2=x*x; x3=x2*x;
    y2=y*y; y3=y2*y;
    z2=z*z; z3=z2*z;

    // 32 elements in a line
    line[0] = 0;
    line[1] = 0;
    line[2] = 1;
    line[3] = 0;
    line[4] = 0;
    line[5] = x;
    line[6] = 0;
    line[7] = 2*y;
    line[8] = z;
    line[9] = 0;
    line[10] = 0;
    line[11] = x2;
    line[12] = 2*x*y;
    line[13] = 0;
    line[14] = 0;
    line[15] = 3*y2;
    line[16] = 2*y*z;
    line[17] = z2;
    line[18] = 0;
    line[19] = x*z;
    line[20] = x3;
    line[21] = 2*x2*y;
    line[22] = 3*x*y2;
    line[23] = 0;
    line[24] = 0;
    line[25] = 0;
    line[26] = 3*y2*z;
    line[27] = 2*y*z2;
    line[28] = z3;
    line[29] = x2*z;
    line[30] = 2*x*y*z;
    line[31] = x*z2;
}

int Amatrix_ele_assign_dz(double *line, double x, double y, double z)
{
    double x2, y2, z2;
    double x3, y3, z3;

    x2=x*x; x3=x2*x;
    y2=y*y; y3=y2*y;
    z2=z*z; z3=z2*z;

    // 32 elements in a line
    line[0] = 0;
    line[1] = 0;
    line[2] = 0;
    line[3] = 1;
    line[4] = 0;
    line[5] = 0;
    line[6] = x;
    line[7] = 0;
    line[8] = y;
    line[9] = 2*z;
    line[10] = 0;
    line[11] = 0;
    line[12] = 0;
    line[13] = x2;
    line[14] = 2*x*z;
    line[15] = 0;
    line[16] = y2;
    line[17] = 2*y*z;
    line[18] = 3*z2;
    line[19] = x*y;
    line[20] = 0;
    line[21] = 0;
    line[22] = 0;
    line[23] = x3;
    line[24] = 2*x2*z;
    line[25] = 3*x*z2;
    line[26] = y3;
    line[27] = 2*y2*z;
    line[28] = 3*y*z2;
    line[29] = x2*y;
    line[30] = x*y2;
    line[31] = 2*x*y*z;
}



int get_values_from_grid(double fxx,double fyy, double fzz, int type, double *usf, double *fsf)
{
    // use cubic hermite interpolation for energy and force calculations
    /*
     *         H---------G
     *        /         /|
     *       E---------F |
     *       | |       | |
     *  z(v) | D-      | C
     *       |/        |/   y(u)
     *       A---------B
     *
     *          x(t)
     */
    /*
    // values at the 8 corners
    double Avalue, Bvalue, Cvalue, Dvalue;
    double Evalue, Fvalue, Gvalue, Hvalue;
    // 1st and 2nd order derivaties at 8 corners
    double Adx, Ady, Adz;
    double Bdx, Bdy, Bdz;
    double Cdx, Cdy, Cdz;
    double Ddx, Ddy, Ddz;
    double Edx, Edy, Edz;
    double Fdx, Fdy, Fdz;
    double Gdx, Gdy, Gdz;
    double Hdx, Hdy, Hdz;
    // 2nd order derivaties
    double Adxx, Adxy, Adxz, Adyx, Adyy, Adyz, Adzx, Adzy, Adzz;
    double Bdxx, Bdxy, Bdxz, Bdyx, Bdyy, Bdyz, Bdzx, Bdzy, Bdzz;
    double Cdxx, Cdxy, Cdxz, Cdyx, Cdyy, Cdyz, Cdzx, Cdzy, Cdzz;;
    double Ddxx, Ddxy, Ddxz, Ddyx, Ddyy, Ddyz, Ddzx, Ddzy, Ddzz;;
    double Edxx, Edxy, Edxz, Edyx, Edyy, Edyz, Edzx, Edzy, Edzz;
    double Fdxx, Fdxy, Fdxz, Fdyx, Fdyy, Fdyz, Fdzx, Fdzy, Fdzz;;
    double Gdxx, Gdxy, Gdxz, Gdyx, Gdyy, Gdyz, Gdzx, Gdzy, Gdzz;
    double Hdxx, Hdxy, Hdxz, Hdyx, Hdyy, Hdyz, Hdzx, Hdzy, Hdzz;
    // hermite functions on x,y,z directions
    double h00x, h10x, h01x, h11x;
    double h00y, h10y, h01y, h11y;
    double h00z, h10z, h01z, h11z;
    // interpolated values at 12 boarders
    double H1, // A-B
    H2, // A-D
    H3, // A-E
    H4, // D-C
    H5, // D-H
    H6, // B-C
    H7, // B-F
    H8, // C-G
    H9, // E-F
    H10, // E-H
    H11, // H-G
    H12; // F-G
    // factors for the 12 interpolated values
    double f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12;
    // final interpolate value
    double PP;
    // index of the lowest corner A
    int iAx, iAy, iAz;
    // index of the 8 corners in the grid
    int Aidx, Bidx, Cidx, Didx, Eidx, Fidx, Gidx, Hidx;
    // the fractional increament in x,y,z
    double tx, uy, vz;
    double tx2, tx3, uy2, uy3, vz2, vz3;
    double tempidx;

    // fxx = 9.8;
    // fyy = 10.6;
    // fzz = 0.0;



    // calculate the index of corner A
    tx = modf((fxx-xmin)/grid_itvl_x,&tempidx);
    iAx = (int)tempidx;
    uy = modf((fyy-ymin)/grid_itvl_y,&tempidx);
    iAy = (int)tempidx;
    // periodical boundary
    fzz = fzz - zmax*rint(fzz/zmax);
    // make sure its in the primal box
    if (fzz<0)
    fzz = zmax + fzz;
    vz = modf((fzz-zmin)/grid_itvl_z,&tempidx);
    iAz = (int)tempidx;
    tx2 = tx*tx;
    tx3 = tx2*tx;
    uy2 = uy*uy;
    uy3 = uy2*uy;
    vz2 = vz*vz;
    vz3 = vz2*vz;

    // coordinates of 8 corners
    double Ax, Ay, Az;
    double Bx, By, Bz;
    double Cx, Cy, Cz;
    double Dx, Dy, Dz;
    double Ex, Ey, Ez;
    double Fx, Fy, Fz;
    double Gx, Gy, Gz;
    double Hx, Hy, Hz;
    // calculate the corner coordinates
    Ax = iAx*grid_itvl_x + xmin;
    Ay = iAy*grid_itvl_y + ymin;
    Az = iAz*grid_itvl_z + zmin;
    Bx = Ax + grid_itvl_x; By = Ay; Bz = Az;
    Cx = Ax + grid_itvl_x; Cy = Ay + grid_itvl_y; Cz = Az;
    Dx = Ax; Dy = Ay + grid_itvl_y; Dz = Az;
    Ex = Ax; Ey = Ay; Ez = Az + grid_itvl_z;
    Fx = Ax + grid_itvl_x; Fy = Ay; Fz = Az + grid_itvl_z;
    Gx = Ax + grid_itvl_x; Gy = Ay + grid_itvl_y; Gz = Az + grid_itvl_z;
    Hx = Ax; Hy = Ay + grid_itvl_y; Hz = Az + grid_itvl_z;

    // calculate the index of the 8 corners in the grid
    Aidx = iAx*ngrid_y*ngrid_z + iAy*ngrid_z + iAz;
    Bidx = (iAx+1)*ngrid_y*ngrid_z + iAy*ngrid_z + iAz;
    Cidx = (iAx+1)*ngrid_y*ngrid_z + (iAy+1)*ngrid_z + iAz;
    Didx = iAx*ngrid_y*ngrid_z + (iAy+1)*ngrid_z + iAz;
    Eidx = iAx*ngrid_y*ngrid_z + iAy*ngrid_z + (iAz+1);
    Fidx = (iAx+1)*ngrid_y*ngrid_z + iAy*ngrid_z + (iAz+1);
    Gidx = (iAx+1)*ngrid_y*ngrid_z + (iAy+1)*ngrid_z + (iAz+1);
    Hidx = iAx*ngrid_y*ngrid_z + (iAy+1)*ngrid_z + (iAz+1);

    // make sure its within the boundary
    assert(Hidx<ngrid_total);

    // get the values
    Avalue = Ene[type][Aidx];
    Bvalue = Ene[type][Bidx];
    Cvalue = Ene[type][Cidx];
    Dvalue = Ene[type][Didx];
    Evalue = Ene[type][Eidx];
    Fvalue = Ene[type][Fidx];
    Gvalue = Ene[type][Gidx];
    Hvalue = Ene[type][Hidx];

    // get the 1st order derivatives
    Adx = dEx[type][Aidx];
    Ady = dEy[type][Aidx];
    Adz = dEz[type][Aidx];
    Bdx = dEx[type][Bidx];
    Bdy = dEy[type][Bidx];
    Bdz = dEz[type][Bidx];
    Cdx = dEx[type][Cidx];
    Cdy = dEy[type][Cidx];
    Cdz = dEz[type][Cidx];
    Ddx = dEx[type][Didx];
    Ddy = dEy[type][Didx];
    Ddz = dEz[type][Didx];
    Edx = dEx[type][Eidx];
    Edy = dEy[type][Eidx];
    Edz = dEz[type][Eidx];
    Fdx = dEx[type][Fidx];
    Fdy = dEy[type][Fidx];
    Fdz = dEz[type][Fidx];
    Gdx = dEx[type][Gidx];
    Gdy = dEy[type][Gidx];
    Gdz = dEz[type][Gidx];
    Hdx = dEx[type][Hidx];
    Hdy = dEy[type][Hidx];
    Hdz = dEz[type][Hidx];

    // get the 2nd order derivaties
    Adxx = dFxx[type][Aidx];
    Adxy = dFxy[type][Aidx];
    Adxz = dFxz[type][Aidx];
    Adyx = dFyx[type][Aidx];
    Adyy = dFyy[type][Aidx];
    Adyz = dFyz[type][Aidx];
    Adzx = dFzx[type][Aidx];
    Adzy = dFzy[type][Aidx];
    Adzz = dFzz[type][Aidx];

    Bdxx = dFxx[type][Bidx];
    Bdxy = dFxy[type][Bidx];
    Bdxz = dFxz[type][Bidx];
    Bdyx = dFyx[type][Bidx];
    Bdyy = dFyy[type][Bidx];
    Bdyz = dFyz[type][Bidx];
    Bdzx = dFzx[type][Bidx];
    Bdzy = dFzy[type][Bidx];
    Bdzz = dFzz[type][Bidx];

    Cdxx = dFxx[type][Cidx];
    Cdxy = dFxy[type][Cidx];
    Cdxz = dFxz[type][Cidx];
    Cdyx = dFyx[type][Cidx];
    Cdyy = dFyy[type][Cidx];
    Cdyz = dFyz[type][Cidx];
    Cdzx = dFzx[type][Cidx];
    Cdzy = dFzy[type][Cidx];
    Cdzz = dFzz[type][Cidx];

    Ddxx = dFxx[type][Didx];
    Ddxy = dFxy[type][Didx];
    Ddxz = dFxz[type][Didx];
    Ddyx = dFyx[type][Didx];
    Ddyy = dFyy[type][Didx];
    Ddyz = dFyz[type][Didx];
    Ddzx = dFzx[type][Didx];
    Ddzy = dFzy[type][Didx];
    Ddzz = dFzz[type][Didx];

    Edxx = dFxx[type][Eidx];
    Edxy = dFxy[type][Eidx];
    Edxz = dFxz[type][Eidx];
    Edyx = dFyx[type][Eidx];
    Edyy = dFyy[type][Eidx];
    Edyz = dFyz[type][Eidx];
    Edzx = dFzx[type][Eidx];
    Edzy = dFzy[type][Eidx];
    Edzz = dFzz[type][Eidx];

    Fdxx = dFxx[type][Fidx];
    Fdxy = dFxy[type][Fidx];
    Fdxz = dFxz[type][Fidx];
    Fdyx = dFyx[type][Fidx];
    Fdyy = dFyy[type][Fidx];
    Fdyz = dFyz[type][Fidx];
    Fdzx = dFzx[type][Fidx];
    Fdzy = dFzy[type][Fidx];
    Fdzz = dFzz[type][Fidx];

    Gdxx = dFxx[type][Gidx];
    Gdxy = dFxy[type][Gidx];
    Gdxz = dFxz[type][Gidx];
    Gdyx = dFyx[type][Gidx];
    Gdyy = dFyy[type][Gidx];
    Gdyz = dFyz[type][Gidx];
    Gdzx = dFzx[type][Gidx];
    Gdzy = dFzy[type][Gidx];
    Gdzz = dFzz[type][Gidx];

    Hdxx = dFxx[type][Hidx];
    Hdxy = dFxy[type][Hidx];
    Hdxz = dFxz[type][Hidx];
    Hdyx = dFyx[type][Hidx];
    Hdyy = dFyy[type][Hidx];
    Hdyz = dFyz[type][Hidx];
    Hdzx = dFzx[type][Hidx];
    Hdzy = dFzy[type][Hidx];
    Hdzz = dFzz[type][Hidx];

    // calculate hermite functions in X direction
    h00x = 2.0*tx3 - 3.0*tx2 + 1.0;
    h10x = (tx3 - 2.0*tx2 + tx)*grid_itvl_x;
    h01x = -2.0*tx3 + 3.0*tx2;
    h11x = (tx3 - tx2)*grid_itvl_x;
    // calculate hermite functions in Y direction
    h00y = 2.0*uy3 - 3.0*uy2 + 1.0;
    h10y = (uy3 - 2.0*uy2 + uy)*grid_itvl_y;
    h01y = -2.0*uy3 + 3.0*uy2;
    h11y = (uy3 - uy2)*grid_itvl_y;
    // calculate hermite functions in X direction
    h00z = 2.0*vz3 - 3.0*vz2 + 1.0;
    h10z = (vz3 - 2.0*vz2 + vz)*grid_itvl_z;
    h01z = -2.0*vz3 + 3.0*vz2;
    h11z = (vz3 - vz2)*grid_itvl_z;

    // calculate the factors
    f1 = (1-uy)*(1-vz);
    f2 = (1-tx)*(1-vz);
    f3 = (1-uy)*(1-tx);
    f4 = uy*(1-vz);
    f5 = uy*(1-tx);
    f6 = tx*(1-vz);
    f7 = tx*(1-uy);
    f8 = tx*uy;
    f9 = vz*(1-uy);
    f10 = vz*(1-tx);
    f11 = uy*vz;
    f12 = vz*tx;

    // Energy
    // calculate the 12 intermediate interpolations
    H1 = h00x*Avalue + h10x*Adx + h01x*Bvalue + h11x*Bdx; // A-B
    H2 = h00y*Avalue + h10y*Ady + h01y*Dvalue + h11y*Ddy; // A-D
    H3 = h00z*Avalue + h10z*Adz + h01z*Evalue + h11z*Edz; // A-E
    H4 = h00x*Dvalue + h10x*Ddx + h01x*Cvalue + h11x*Cdx; // D-C
    H5 = h00z*Dvalue + h10z*Ddz + h01z*Hvalue + h11z*Hdz; // D-H
    H6 = h00y*Bvalue + h10y*Bdy + h01y*Cvalue + h11y*Cdy; // B-C
    H7 = h00z*Bvalue + h10z*Bdz + h01z*Fvalue + h11z*Fdz; // B-F
    H8 = h00z*Cvalue + h10z*Cdz + h01z*Gvalue + h11z*Gdz; // C-G
    H9 = h00x*Evalue + h10x*Edx + h01x*Fvalue + h11x*Fdx; // E-F
    H10 = h00y*Evalue + h10y*Edy + h01y*Hvalue + h11y*Hdy; // E-H
    H11 = h00x*Hvalue + h10x*Hdx + h01x*Gvalue + h11x*Gdx; // H-G
    H12 = h00y*Fvalue + h10y*Fdy + h01y*Gvalue + h11y*Gdy; // F-G
    // final value 
    PP = (H1*f1+H2*f2+H3*f3+H4*f4+H5*f5+H6*f6+H7*f7+H8*f8+H9*f9+H10*f10+H11*f11+H12*f12)/3.0;
    *usf = PP;

    // fx
    // calculate the 12 intermediate interpolations
    H1 = h00x*(-Adx) + h10x*Adxx + h01x*(-Bdx) + h11x*Bdxx; // A-B
    H2 = h00y*(-Adx) + h10y*Adxy + h01y*(-Ddx) + h11y*Ddxy; // A-D
    H3 = h00z*(-Adx) + h10z*Adxz + h01z*(-Edx) + h11z*Edxz; // A-E
    H4 = h00x*(-Ddx) + h10x*Ddxx + h01x*(-Cdx) + h11x*Cdxx; // D-C
    H5 = h00z*(-Ddx) + h10z*Ddxz + h01z*(-Hdx) + h11z*Hdxz; // D-H
    H6 = h00y*(-Bdx) + h10y*Bdxy + h01y*(-Cdx) + h11y*Cdxy; // B-C
    H7 = h00z*(-Bdx) + h10z*Bdxz + h01z*(-Fdx) + h11z*Fdxz; // B-F
    H8 = h00z*(-Cdx) + h10z*Cdxz + h01z*(-Gdx) + h11z*Gdxz; // C-G
    H9 = h00x*(-Edx) + h10x*Edxx + h01x*(-Fdx) + h11x*Fdxx; // E-F
    H10 = h00y*(-Edx) + h10y*Edxy + h01y*(-Hdx) + h11y*Hdxy; // E-H
    H11 = h00x*(-Hdx) + h10x*Hdxx + h01x*(-Gdx) + h11x*Gdxx; // H-G
    H12 = h00y*(-Fdx) + h10y*Fdxy + h01y*(-Gdx) + h11y*Gdxy; // F-G
    // final value 
    PP = (H1*f1+H2*f2+H3*f3+H4*f4+H5*f5+H6*f6+H7*f7+H8*f8+H9*f9+H10*f10+H11*f11+H12*f12)/3.0;
    fsf[0] = PP;

    // fy
    // calculate the 12 intermediate interpolations
    H1 = h00x*(-Ady) + h10x*Adyx + h01x*(-Bdy) + h11x*Bdyx; // A-B
    H2 = h00y*(-Ady) + h10y*Adyy + h01y*(-Ddy) + h11y*Ddyy; // A-D
    H3 = h00z*(-Ady) + h10z*Adyz + h01z*(-Edy) + h11z*Edyz; // A-E
    H4 = h00x*(-Ddy) + h10x*Ddyx + h01x*(-Cdy) + h11x*Cdyx; // D-C
    H5 = h00z*(-Ddy) + h10z*Ddyz + h01z*(-Hdy) + h11z*Hdyz; // D-H
    H6 = h00y*(-Bdy) + h10y*Bdyy + h01y*(-Cdy) + h11y*Cdyy; // B-C
    H7 = h00z*(-Bdy) + h10z*Bdyz + h01z*(-Fdy) + h11z*Fdyz; // B-F
    H8 = h00z*(-Cdy) + h10z*Cdyz + h01z*(-Gdy) + h11z*Gdyz; // C-G
    H9 = h00x*(-Edy) + h10x*Edyx + h01x*(-Fdy) + h11x*Fdyx; // E-F
    H10 = h00y*(-Edy) + h10y*Edyy + h01y*(-Hdy) + h11y*Hdyy; // E-H
    H11 = h00x*(-Hdy) + h10x*Hdyx + h01x*(-Gdy) + h11x*Gdyx; // H-G
    H12 = h00y*(-Fdy) + h10y*Fdyy + h01y*(-Gdy) + h11y*Gdyy; // F-G
    // final value 
    PP = (H1*f1+H2*f2+H3*f3+H4*f4+H5*f5+H6*f6+H7*f7+H8*f8+H9*f9+H10*f10+H11*f11+H12*f12)/3.0;
    fsf[1] = PP;

    // fz
    // calculate the 12 intermediate interpolations
    H1 = h00x*(-Adz) + h10x*Adzx + h01x*(-Bdz) + h11x*Bdzx; // A-B
    H2 = h00y*(-Adz) + h10y*Adzy + h01y*(-Ddz) + h11y*Ddzy; // A-D
    H3 = h00z*(-Adz) + h10z*Adzz + h01z*(-Edz) + h11z*Edzz; // A-E
    H4 = h00x*(-Ddz) + h10x*Ddzx + h01x*(-Cdz) + h11x*Cdzx; // D-C
    H5 = h00z*(-Ddz) + h10z*Ddzz + h01z*(-Hdz) + h11z*Hdzz; // D-H
    H6 = h00y*(-Bdz) + h10y*Bdzy + h01y*(-Cdz) + h11y*Cdzy; // B-C
    H7 = h00z*(-Bdz) + h10z*Bdzz + h01z*(-Fdz) + h11z*Fdzz; // B-F
    H8 = h00z*(-Cdz) + h10z*Cdzz + h01z*(-Gdz) + h11z*Gdzz; // C-G
    H9 = h00x*(-Edz) + h10x*Edzx + h01x*(-Fdz) + h11x*Fdzx; // E-F
    H10 = h00y*(-Edz) + h10y*Edzy + h01y*(-Hdz) + h11y*Hdzy; // E-H
    H11 = h00x*(-Hdz) + h10x*Hdzx + h01x*(-Gdz) + h11x*Gdzx; // H-G
    H12 = h00y*(-Fdz) + h10y*Fdzy + h01y*(-Gdz) + h11y*Gdzy; // F-G
    // final value 
    PP = (H1*f1+H2*f2+H3*f3+H4*f4+H5*f5+H6*f6+H7*f7+H8*f8+H9*f9+H10*f10+H11*f11+H12*f12)/3.0;
    fsf[2] = PP;

    */

	/*
	   printf("xx=%lf  yy=%lf  zz=%lf\n",fxx,fyy,fzz);
	   printf("iAx=%d  iAy=%d  iAz=%d\n",iAx,iAy,iAz);
	   printf("Aidx=%d  Eidx=%d  Hidx=%d\n",Aidx,Eidx,Hidx);

	   printf("Ax=%lf  Ay=%lf  Az=%lf\n",Ax,Ay,Az);
	   printf("Bx=%lf  By=%lf  Bz=%lf\n",Bx,By,Bz);
	   printf("Cx=%lf  Cy=%lf  Cz=%lf\n",Cx,Cy,Cz);
	   printf("Dx=%lf  Dy=%lf  Dz=%lf\n",Dx,Dy,Dz);
	   printf("Ex=%lf  Ey=%lf  Ez=%lf\n",Ex,Ey,Ez);
	   printf("Fx=%lf  Fy=%lf  Fz=%lf\n",Fx,Fy,Fz);
	   printf("Gx=%lf  Gy=%lf  Gz=%lf\n",Gx,Gy,Gz);
	   printf("Hx=%lf  Hy=%lf  Hz=%lf\n",Hx,Hy,Hz);

	   printf("A=%lf B=%lf C=%lf D=%lf E=%lf F=%lf G=%lf H=%lf\n",
	   Avalue,Bvalue,Cvalue,Dvalue,Evalue,Fvalue,Gvalue,Hvalue);

	   printf("Adx=%lf Ady=%lf Adz=%lf\n",Adx,Ady,Adz);
	   printf("Bdx=%lf Bdy=%lf Bdz=%lf\n",Bdx,Bdy,Bdz);
	   printf("Cdx=%lf Cdy=%lf Cdz=%lf\n",Cdx,Cdy,Cdz);
	   printf("Ddx=%lf Ddy=%lf Ddz=%lf\n",Ddx,Ddy,Ddz);
	   printf("Edx=%lf Edy=%lf Edz=%lf\n",Edx,Edy,Edz);
	   printf("Fdx=%lf Fdy=%lf Fdz=%lf\n",Fdx,Fdy,Fdz);
	   printf("Gdx=%lf Gdy=%lf Gdz=%lf\n",Gdx,Gdy,Gdz);
	   printf("Hdx=%lf Hdy=%lf Hdz=%lf\n",Hdx,Hdy,Hdz);

	   printf("\n");
	   printf("Adxx=%lf Adxy=%lf Adxz=%lf\n",Adxx,Adxy,Adxz);
	   printf("Adyx=%lf Adyy=%lf Adyz=%lf\n",Adyx,Adyy,Adyz);
	   printf("Adzx=%lf Adzy=%lf Adzz=%lf\n",Adzx,Adzy,Adzz);
	   printf("Bdxx=%lf Bdxy=%lf Bdxz=%lf\n",Bdxx,Bdxy,Bdxz);
	   printf("Bdyx=%lf Bdyy=%lf Bdyz=%lf\n",Bdyx,Bdyy,Bdyz);
	   printf("Bdzx=%lf Bdzy=%lf Bdzz=%lf\n",Bdzx,Bdzy,Bdzz);
	   printf("Cdxx=%lf Cdxy=%lf Cdxz=%lf\n",Cdxx,Cdxy,Cdxz);
	   printf("Cdyx=%lf Cdyy=%lf Cdyz=%lf\n",Cdyx,Cdyy,Cdyz);
	   printf("Cdzx=%lf Cdzy=%lf Cdzz=%lf\n",Cdzx,Cdzy,Cdzz);
	   printf("Ddxx=%lf Ddxy=%lf Ddxz=%lf\n",Ddxx,Ddxy,Ddxz);
	   printf("Ddyx=%lf Ddyy=%lf Ddyz=%lf\n",Ddyx,Ddyy,Ddyz);
	   printf("Ddzx=%lf Ddzy=%lf Ddzz=%lf\n",Ddzx,Ddzy,Ddzz);
	   printf("Edxx=%lf Edxy=%lf Edxz=%lf\n",Edxx,Edxy,Edxz);
	   printf("Edyx=%lf Edyy=%lf Edyz=%lf\n",Edyx,Edyy,Edyz);
	   printf("Edzx=%lf Edzy=%lf Edzz=%lf\n",Edzx,Edzy,Edzz);
	   printf("Fdxx=%lf Fdxy=%lf Fdxz=%lf\n",Fdxx,Fdxy,Fdxz);
	   printf("Fdyx=%lf Fdyy=%lf Fdyz=%lf\n",Fdyx,Fdyy,Fdyz);
	   printf("Fdzx=%lf Fdzy=%lf Fdzz=%lf\n",Fdzx,Fdzy,Fdzz);
	   printf("Gdxx=%lf Gdxy=%lf Gdxz=%lf\n",Gdxx,Gdxy,Gdxz);
	   printf("Gdyx=%lf Gdyy=%lf Gdyz=%lf\n",Gdyx,Gdyy,Gdyz);
	   printf("Gdzx=%lf Gdzy=%lf Gdzz=%lf\n",Gdzx,Gdzy,Gdzz);
	   printf("Hdxx=%lf Hdxy=%lf Hdxz=%lf\n",Hdxx,Hdxy,Hdxz);
	   printf("Hdyx=%lf Hdyy=%lf Hdyz=%lf\n",Hdyx,Hdyy,Hdyz);
	   printf("Hdzx=%lf Hdzy=%lf Hdzz=%lf\n",Hdzx,Hdzy,Hdzz);

	   printf("\nH1=%lf H2=%lf H3=%lf H4=%lf H5=%lf H6=%lf H7=%lf H8=%lf H9=%lf H10=%lf H11=%lf H12=%lf\n",
	   H1,H2,H3,H4,H5,H6,H7,H8,H9,H10,H11,H12);
	   printf("f1=%lf f2=%lf f3=%lf f4=%lf f5=%lf f6=%lf f7=%lf f8=%lf f9=%lf f10=%lf f11=%lf f12=%lf\n",
	   f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12);

	   printf("usf=%lf\n",*usf);
	   printf("fx=%lf   fy=%lf   fz=%lf\n",fsf[0],fsf[1],fsf[2]);

	   exit(1);
	   */










    // values at the 8 corners
    double Avalue, Bvalue, Cvalue, Dvalue;
    double Evalue, Fvalue, Gvalue, Hvalue;
    // 1st order derivaties at 8 corners
    double Adx, Ady, Adz; 
    double Bdx, Bdy, Bdz; 
    double Cdx, Cdy, Cdz; 
    double Ddx, Ddy, Ddz; 
    double Edx, Edy, Edz; 
    double Fdx, Fdy, Fdz; 
    double Gdx, Gdy, Gdz; 
    double Hdx, Hdy, Hdz; 
    // 2nd order derivaties
    double Adxx, Adxy, Adxz, Adyx, Adyy, Adyz, Adzx, Adzy, Adzz;
    double Bdxx, Bdxy, Bdxz, Bdyx, Bdyy, Bdyz, Bdzx, Bdzy, Bdzz;
    double Cdxx, Cdxy, Cdxz, Cdyx, Cdyy, Cdyz, Cdzx, Cdzy, Cdzz;;
    double Ddxx, Ddxy, Ddxz, Ddyx, Ddyy, Ddyz, Ddzx, Ddzy, Ddzz;;
    double Edxx, Edxy, Edxz, Edyx, Edyy, Edyz, Edzx, Edzy, Edzz;
    double Fdxx, Fdxy, Fdxz, Fdyx, Fdyy, Fdyz, Fdzx, Fdzy, Fdzz;;
    double Gdxx, Gdxy, Gdxz, Gdyx, Gdyy, Gdyz, Gdzx, Gdzy, Gdzz;
    double Hdxx, Hdxy, Hdxz, Hdyx, Hdyy, Hdyz, Hdzx, Hdzy, Hdzz;
    // final interpolate value
    double PP;
    // index of the lowest corner A
    int iAx, iAy, iAz;
    // index of the 8 corners in the grid
    int Aidx, Bidx, Cidx, Didx, Eidx, Fidx, Gidx, Hidx;
    // coordinates of 8 corners
    double Ax, Ay, Az;
    double Bx, By, Bz;
    double Cx, Cy, Cz;
    double Dx, Dy, Dz;
    double Ex, Ey, Ez;
    double Fx, Fy, Fz;
    double Gx, Gy, Gz;
    double Hx, Hy, Hz;

    int ii;

    // double amatrix[1024], amatrix_backup[1024];
    // gsl_vector *work = gsl_vector_alloc(32);
    // gsl_vector *Svector = gsl_vector_alloc(32);
    // gsl_vector *xvector = gsl_vector_alloc(32);
    // gsl_matrix *Vmatrix = gsl_matrix_alloc(32,32);
    // double interp_vector[32];
    // double bvector[32];

    // printf("xx=%lf  yy=%lf  zz=%lf\n",fxx,fyy,fzz);

    // calculate the index of corner A
    iAx = (int)((fxx-xmin)/grid_itvl_x);
    iAy = (int)((fyy-ymin)/grid_itvl_y);
    // periodical for z direction
    fzz = fzz - zmax*rint(fzz/zmax);
    // make sure its in the primal box
    if (fzz<0) fzz=zmax+fzz;
    iAz = (int)((fzz-zmin)/grid_itvl_z);


    // calculate the corner coordinates
    Ax = iAx*grid_itvl_x + xmin;
    Ay = iAy*grid_itvl_y + ymin;
    Az = iAz*grid_itvl_z + zmin;
    Bx = Ax + grid_itvl_x; By = Ay; Bz = Az;
    Cx = Ax + grid_itvl_x; Cy = Ay + grid_itvl_y; Cz = Az;
    Dx = Ax; Dy = Ay + grid_itvl_y; Dz = Az;
    Ex = Ax; Ey = Ay; Ez = Az + grid_itvl_z;
    Fx = Ax + grid_itvl_x; Fy = Ay; Fz = Az + grid_itvl_z;
    Gx = Ax + grid_itvl_x; Gy = Ay + grid_itvl_y; Gz = Az + grid_itvl_z;
    Hx = Ax; Hy = Ay + grid_itvl_y; Hz = Az + grid_itvl_z;

    // calculate the index of the 8 corners in the grid
    Aidx = iAx*ngrid_y*ngrid_z + iAy*ngrid_z + iAz;
    Bidx = (iAx+1)*ngrid_y*ngrid_z + iAy*ngrid_z + iAz;
    Cidx = (iAx+1)*ngrid_y*ngrid_z + (iAy+1)*ngrid_z + iAz;
    Didx = iAx*ngrid_y*ngrid_z + (iAy+1)*ngrid_z + iAz;
    Eidx = iAx*ngrid_y*ngrid_z + iAy*ngrid_z + (iAz+1);
    Fidx = (iAx+1)*ngrid_y*ngrid_z + iAy*ngrid_z + (iAz+1);
    Gidx = (iAx+1)*ngrid_y*ngrid_z + (iAy+1)*ngrid_z + (iAz+1);
    Hidx = iAx*ngrid_y*ngrid_z + (iAy+1)*ngrid_z + (iAz+1);

    // printf("xx=%lf  yy=%lf  zz=%lf\n",fxx,fyy,fzz);
    // printf("iAx=%d  iAy=%d  iAz=%d\n",iAx,iAy,iAz);
    // printf("Aidx=%d  Eidx=%d  Hidx=%d\n",Aidx,Eidx,Hidx);

    // make sure its within the boundary
    assert(Hidx<ngrid_total);

    // get the values
    Avalue = Ene[type][Aidx];
    Bvalue = Ene[type][Bidx];
    Cvalue = Ene[type][Cidx];
    Dvalue = Ene[type][Didx];
    Evalue = Ene[type][Eidx];
    Fvalue = Ene[type][Fidx];
    Gvalue = Ene[type][Gidx];
    Hvalue = Ene[type][Hidx];

    // get the 1st order derivatives
    Adx = dEx[type][Aidx];
    Ady = dEy[type][Aidx];
    Adz = dEz[type][Aidx];
    Bdx = dEx[type][Bidx];
    Bdy = dEy[type][Bidx];
    Bdz = dEz[type][Bidx];
    Cdx = dEx[type][Cidx];
    Cdy = dEy[type][Cidx];
    Cdz = dEz[type][Cidx];
    Ddx = dEx[type][Didx];
    Ddy = dEy[type][Didx];
    Ddz = dEz[type][Didx];
    Edx = dEx[type][Eidx];
    Edy = dEy[type][Eidx];
    Edz = dEz[type][Eidx];
    Fdx = dEx[type][Fidx];
    Fdy = dEy[type][Fidx];
    Fdz = dEz[type][Fidx];
    Gdx = dEx[type][Gidx];
    Gdy = dEy[type][Gidx];
    Gdz = dEz[type][Gidx];
    Hdx = dEx[type][Hidx];
    Hdy = dEy[type][Hidx];
    Hdz = dEz[type][Hidx];

    // get 2nd order derivaties
    Adxx = dFxx[type][Aidx];
    Adxy = dFxy[type][Aidx];
    Adxz = dFxz[type][Aidx];
    Adyx = dFyx[type][Aidx];
    Adyy = dFyy[type][Aidx];
    Adyz = dFyz[type][Aidx];
    Adzx = dFzx[type][Aidx];
    Adzy = dFzy[type][Aidx];
    Adzz = dFzz[type][Aidx];

    Bdxx = dFxx[type][Bidx];
    Bdxy = dFxy[type][Bidx];
    Bdxz = dFxz[type][Bidx];
    Bdyx = dFyx[type][Bidx];
    Bdyy = dFyy[type][Bidx];
    Bdyz = dFyz[type][Bidx];
    Bdzx = dFzx[type][Bidx];
    Bdzy = dFzy[type][Bidx];
    Bdzz = dFzz[type][Bidx];

    Cdxx = dFxx[type][Cidx];
    Cdxy = dFxy[type][Cidx];
    Cdxz = dFxz[type][Cidx];
    Cdyx = dFyx[type][Cidx];
    Cdyy = dFyy[type][Cidx];
    Cdyz = dFyz[type][Cidx];
    Cdzx = dFzx[type][Cidx];
    Cdzy = dFzy[type][Cidx];
    Cdzz = dFzz[type][Cidx];

    Ddxx = dFxx[type][Didx];
    Ddxy = dFxy[type][Didx];
    Ddxz = dFxz[type][Didx];
    Ddyx = dFyx[type][Didx];
    Ddyy = dFyy[type][Didx];
    Ddyz = dFyz[type][Didx];
    Ddzx = dFzx[type][Didx];
    Ddzy = dFzy[type][Didx];
    Ddzz = dFzz[type][Didx];

    Edxx = dFxx[type][Eidx];
    Edxy = dFxy[type][Eidx];
    Edxz = dFxz[type][Eidx];
    Edyx = dFyx[type][Eidx];
    Edyy = dFyy[type][Eidx];
    Edyz = dFyz[type][Eidx];
    Edzx = dFzx[type][Eidx];
    Edzy = dFzy[type][Eidx];
    Edzz = dFzz[type][Eidx];

    Fdxx = dFxx[type][Fidx];
    Fdxy = dFxy[type][Fidx];
    Fdxz = dFxz[type][Fidx];
    Fdyx = dFyx[type][Fidx];
    Fdyy = dFyy[type][Fidx];
    Fdyz = dFyz[type][Fidx];
    Fdzx = dFzx[type][Fidx];
    Fdzy = dFzy[type][Fidx];
    Fdzz = dFzz[type][Fidx];

    Gdxx = dFxx[type][Gidx];
    Gdxy = dFxy[type][Gidx];
    Gdxz = dFxz[type][Gidx];
    Gdyx = dFyx[type][Gidx];
    Gdyy = dFyy[type][Gidx];
    Gdyz = dFyz[type][Gidx];
    Gdzx = dFzx[type][Gidx];
    Gdzy = dFzy[type][Gidx];
    Gdzz = dFzz[type][Gidx];

    Hdxx = dFxx[type][Hidx];
    Hdxy = dFxy[type][Hidx];
    Hdxz = dFxz[type][Hidx];
    Hdyx = dFyx[type][Hidx];
    Hdyy = dFyy[type][Hidx];
    Hdyz = dFyz[type][Hidx];
    Hdzx = dFzx[type][Hidx];
    Hdzy = dFzy[type][Hidx];
    Hdzz = dFzz[type][Hidx];


    /*
    Ax=3; Ay=3; Az=3;
    Bx=3.2; By=3; Bz=3;
    Cx=3.2; Cy=3.2; Cz=3;
    Dx=3; Dy=3.2; Dz=3;
    Ex=3; Ey=3; Ez=3.2; 
    Fx=3.2; Fy=3; Fz=3.2;
    Gx=3.2; Gy=3.2; Gz=3.2;
    Hx=3; Hy=3.2; Hz=3.2;
    fxx=3.1; fyy=3.1; fzz=3.1;
    */

    /*
       Avalue = 27.000000;
       Bvalue = 28.240000;
       Cvalue = 29.480000;
       Dvalue = 28.240000;
       Evalue = 28.240000;
       Fvalue = 29.480000;
       Gvalue = 30.720000;
       Hvalue = 29.480000;

       Adx=6.00000000; Ady=6.00000000; Adz=6.00000000;
       Bdx=6.40000000; Bdy=6.00000000; Bdz=6.00000000;
       Cdx=6.40000000; Cdy=6.40000000; Cdz=6.00000000;
       Ddx=6.00000000; Ddy=6.40000000; Ddz=6.00000000;
       Edx=6.00000000; Edy=6.00000000; Edz=6.40000000;
       Fdx=6.40000000; Fdy=6.00000000; Fdz=6.40000000;
       Gdx=6.40000000; Gdy=6.40000000; Gdz=6.40000000;
       Hdx=6.00000000; Hdy=6.40000000; Hdz=6.40000000;
       */



    /*
       Avalue = -165.1804;
       Bvalue = -146.5216;
       Cvalue = -130.3909;
       Dvalue = -146.5216;
       Evalue = -146.5216;
       Fvalue = -130.3909;
       Gvalue = -116.4187;
       Hvalue = -130.3909;

       Adx=97.0547788; Ady=97.0547788; Adz=97.0547788;
       Bdx=103.525097; Bdy=97.0547788; Bdz=97.0547788;
       Cdx=103.525097; Cdy=103.525097; Cdz=97.0547788;
       Ddx=97.0547788; Ddy=103.525097; Ddz=97.0547788;
       Edx=97.0547788; Edy=97.0547788; Edz=103.525097;
       Fdx=103.525097; Fdy=97.0547788; Fdz=103.525097;
       Gdx=103.525097; Gdy=103.525097; Gdz=103.525097;
       Hdx=97.0547788; Hdy=103.525097; Hdz=103.525097;
       */



    /*
       printf("Ax=%lf  Ay=%lf  Az=%lf\n",Ax,Ay,Az);

       printf("A=%lf B=%lf C=%lf D=%lf E=%lf F=%lf G=%lf H=%lf\n",
       Avalue,Bvalue,Cvalue,Dvalue,Evalue,Fvalue,Gvalue,Hvalue);

       printf("Adx=%lf Ady=%lf Adz=%lf\n",Adx,Ady,Adz);
       printf("Bdx=%lf Bdy=%lf Bdz=%lf\n",Bdx,Bdy,Bdz);
       printf("Cdx=%lf Cdy=%lf Cdz=%lf\n",Cdx,Cdy,Cdz);
       printf("Ddx=%lf Ddy=%lf Ddz=%lf\n",Ddx,Ddy,Ddz);
       printf("Edx=%lf Edy=%lf Edz=%lf\n",Edx,Edy,Edz);
       printf("Fdx=%lf Fdy=%lf Fdz=%lf\n",Fdx,Fdy,Fdz);
       printf("Gdx=%lf Gdy=%lf Gdz=%lf\n",Gdx,Gdy,Gdz);
       printf("Hdx=%lf Hdy=%lf Hdz=%lf\n",Hdx,Hdy,Hdz);

	   printf("\n");
	   printf("Adxx=%lf Adxy=%lf Adxz=%lf\n",Adxx,Adxy,Adxz);
	   printf("Adyx=%lf Adyy=%lf Adyz=%lf\n",Adyx,Adyy,Adyz);
	   printf("Adzx=%lf Adzy=%lf Adzz=%lf\n",Adzx,Adzy,Adzz);
	   printf("Bdxx=%lf Bdxy=%lf Bdxz=%lf\n",Bdxx,Bdxy,Bdxz);
	   printf("Bdyx=%lf Bdyy=%lf Bdyz=%lf\n",Bdyx,Bdyy,Bdyz);
	   printf("Bdzx=%lf Bdzy=%lf Bdzz=%lf\n",Bdzx,Bdzy,Bdzz);
	   printf("Cdxx=%lf Cdxy=%lf Cdxz=%lf\n",Cdxx,Cdxy,Cdxz);
	   printf("Cdyx=%lf Cdyy=%lf Cdyz=%lf\n",Cdyx,Cdyy,Cdyz);
	   printf("Cdzx=%lf Cdzy=%lf Cdzz=%lf\n",Cdzx,Cdzy,Cdzz);
	   printf("Ddxx=%lf Ddxy=%lf Ddxz=%lf\n",Ddxx,Ddxy,Ddxz);
	   printf("Ddyx=%lf Ddyy=%lf Ddyz=%lf\n",Ddyx,Ddyy,Ddyz);
	   printf("Ddzx=%lf Ddzy=%lf Ddzz=%lf\n",Ddzx,Ddzy,Ddzz);
	   printf("Edxx=%lf Edxy=%lf Edxz=%lf\n",Edxx,Edxy,Edxz);
	   printf("Edyx=%lf Edyy=%lf Edyz=%lf\n",Edyx,Edyy,Edyz);
	   printf("Edzx=%lf Edzy=%lf Edzz=%lf\n",Edzx,Edzy,Edzz);
	   printf("Fdxx=%lf Fdxy=%lf Fdxz=%lf\n",Fdxx,Fdxy,Fdxz);
	   printf("Fdyx=%lf Fdyy=%lf Fdyz=%lf\n",Fdyx,Fdyy,Fdyz);
	   printf("Fdzx=%lf Fdzy=%lf Fdzz=%lf\n",Fdzx,Fdzy,Fdzz);
	   printf("Gdxx=%lf Gdxy=%lf Gdxz=%lf\n",Gdxx,Gdxy,Gdxz);
	   printf("Gdyx=%lf Gdyy=%lf Gdyz=%lf\n",Gdyx,Gdyy,Gdyz);
	   printf("Gdzx=%lf Gdzy=%lf Gdzz=%lf\n",Gdzx,Gdzy,Gdzz);
	   printf("Hdxx=%lf Hdxy=%lf Hdxz=%lf\n",Hdxx,Hdxy,Hdxz);
	   printf("Hdyx=%lf Hdyy=%lf Hdyz=%lf\n",Hdyx,Hdyy,Hdyz);
	   printf("Hdzx=%lf Hdzy=%lf Hdzz=%lf\n",Hdzx,Hdzy,Hdzz);
	   */


    // calcualte the element values in matrix A
    // Ax=b
    // Matrix A is ordered as
    // 1-8 	are A-H values
    // 9-16 	are A-H dx 
    // 17-24	are A-H dy 
    // 25-32	are A-H dz 
    Amatrix_ele_assign_value(amatrix,Ax,Ay,Az);
    Amatrix_ele_assign_value(amatrix+32,Bx,By,Bz);
    Amatrix_ele_assign_value(amatrix+64,Cx,Cy,Cz);
    Amatrix_ele_assign_value(amatrix+96,Dx,Dy,Dz);
    Amatrix_ele_assign_value(amatrix+128,Ex,Ey,Ez);
    Amatrix_ele_assign_value(amatrix+160,Fx,Fy,Fz);
    Amatrix_ele_assign_value(amatrix+192,Gx,Gy,Gz);
    Amatrix_ele_assign_value(amatrix+224,Hx,Hy,Hz);
    Amatrix_ele_assign_dx(amatrix+256,Ax,Ay,Az);
    Amatrix_ele_assign_dx(amatrix+288,Bx,By,Bz);
    Amatrix_ele_assign_dx(amatrix+320,Cx,Cy,Cz);
    Amatrix_ele_assign_dx(amatrix+352,Dx,Dy,Dz);
    Amatrix_ele_assign_dx(amatrix+384,Ex,Ey,Ez);
    Amatrix_ele_assign_dx(amatrix+416,Fx,Fy,Fz);
    Amatrix_ele_assign_dx(amatrix+448,Gx,Gy,Gz);
    Amatrix_ele_assign_dx(amatrix+480,Hx,Hy,Hz);
    Amatrix_ele_assign_dy(amatrix+512,Ax,Ay,Az);
    Amatrix_ele_assign_dy(amatrix+544,Bx,By,Bz);
    Amatrix_ele_assign_dy(amatrix+576,Cx,Cy,Cz);
    Amatrix_ele_assign_dy(amatrix+608,Dx,Dy,Dz);
    Amatrix_ele_assign_dy(amatrix+640,Ex,Ey,Ez);
    Amatrix_ele_assign_dy(amatrix+672,Fx,Fy,Fz);
    Amatrix_ele_assign_dy(amatrix+704,Gx,Gy,Gz);
    Amatrix_ele_assign_dy(amatrix+736,Hx,Hy,Hz);
    Amatrix_ele_assign_dz(amatrix+768,Ax,Ay,Az);
    Amatrix_ele_assign_dz(amatrix+800,Bx,By,Bz);
    Amatrix_ele_assign_dz(amatrix+832,Cx,Cy,Cz);
    Amatrix_ele_assign_dz(amatrix+864,Dx,Dy,Dz);
    Amatrix_ele_assign_dz(amatrix+896,Ex,Ey,Ez);
    Amatrix_ele_assign_dz(amatrix+928,Fx,Fy,Fz);
    Amatrix_ele_assign_dz(amatrix+960,Gx,Gy,Gz);
    Amatrix_ele_assign_dz(amatrix+992,Hx,Hy,Hz);

    // set the backup matrix since the SVD function will change A's value
    for (ii=0;ii<1024;ii++) amatrix_backup[ii] = amatrix[ii];

    // assign the value of the interpolation vector at x,y,z
    Amatrix_ele_assign_value(interp_vector,fxx,fyy,fzz);

    // printf("interp =\n"); for (ii=0;ii<32;ii++) printf("%lf\n",interp_vector[ii]);


    /*
       int jj;
       for (ii=0;ii<32;ii++)
       {
       for (jj=0;jj<32;jj++)
       {
       printf("%7.4lf ",amatrix[ii*32+jj]);
       }
       printf("\n");
       }
       */

    // assign data to b vector for energy
    // ordered as A-H, value, dx, dy, dz
    // This is the only thing that needs to be changed
    // for energy and force interpolations
    bvector[0] = Avalue;
    bvector[1] = Bvalue;
    bvector[2] = Cvalue;
    bvector[3] = Dvalue;
    bvector[4] = Evalue;
    bvector[5] = Fvalue;
    bvector[6] = Gvalue;
    bvector[7] = Hvalue;
    bvector[8] = Adx;
    bvector[9] = Bdx;
    bvector[10] = Cdx;
    bvector[11] = Ddx;
    bvector[12] = Edx;
    bvector[13] = Fdx;
    bvector[14] = Gdx;
    bvector[15] = Hdx;
    bvector[16] = Ady;
    bvector[17] = Bdy;
    bvector[18] = Cdy;
    bvector[19] = Ddy;
    bvector[20] = Edy;
    bvector[21] = Fdy;
    bvector[22] = Gdy;
    bvector[23] = Hdy;
    bvector[24] = Adz;
    bvector[25] = Bdz;
    bvector[26] = Cdz;
    bvector[27] = Ddz;
    bvector[28] = Edz;
    bvector[29] = Fdz;
    bvector[30] = Gdz;
    bvector[31] = Hdz;

    // printf("b = \n"); for (ii=0;ii<32;ii++) printf("%lf\n",bvector[ii]);

    // create matrix view object for amatrix
    gsl_matrix_view Umatrix
	= gsl_matrix_view_array(amatrix, 32, 32);

    // create vector view object for bvector
    gsl_vector_view b
	= gsl_vector_view_array(bvector, 32);

    // singular value decomposition
    // gsl_linalg_SV_decomp(&Umatrix.matrix, Vmatrix, Svector, work);
    gsl_linalg_SV_decomp_mod(&Umatrix.matrix, Xmatrix, Vmatrix, Svector, work);
    // modified Golub-Reinsch algorithm, faster for M>>N

    // printf("s = \n"); gsl_vector_fprintf(stdout,Svector,"%lf");

    // set the trivial values of svector to exactly 0.0 to avoild 
    // numerical accuracy problem. very important!
    for (ii=0;ii<32;ii++)
    {
	if (Svector->data[ii] <1.0e-8)
	    Svector->data[ii] = 0.0;
    }

    // printf("s = \n"); gsl_vector_fprintf(stdout,Svector,"%lf");

    // solve the equation of Ax=b using singular vectors
    // minimizing least squares
    gsl_linalg_SV_solve(&Umatrix.matrix, Vmatrix, Svector, &b.vector, xvector);

    // calculate the value at the interpolatin position
    PP = 0.0;
    for (ii=0;ii<32;ii++)
	PP += interp_vector[ii]*xvector->data[ii];
    // assign the value
    *usf = PP;

    // assign the value of the interpolation vector at x,y,z
    Amatrix_ele_assign_dx(interp_vector,fxx,fyy,fzz);
    PP = 0.0;
    for (ii=0;ii<32;ii++)
	PP += interp_vector[ii]*xvector->data[ii];
    fsf[0] = -PP;
    // printf("dx=%lf\n",-PP);
    Amatrix_ele_assign_dy(interp_vector,fxx,fyy,fzz);
    PP = 0.0;
    for (ii=0;ii<32;ii++)
	PP += interp_vector[ii]*xvector->data[ii];
    fsf[1] = -PP;
    // printf("dy=%lf\n",-PP);
    Amatrix_ele_assign_dz(interp_vector,fxx,fyy,fzz);
    PP = 0.0;
    for (ii=0;ii<32;ii++)
	PP += interp_vector[ii]*xvector->data[ii];
    fsf[2] = -PP;
    // printf("dz=%lf\n",-PP);
    // Amatrix_ele_assign_value(interp_vector,fxx,fyy,fzz);


    /*
    // ---------------- fx -------------------
    // assign values to bvector for fx
    // ordered as A-H, value, dx, dy, dz
    // This is the only thing that needs to be changed
    // for energy and force interpolations
    bvector[0] = -Adx;
    bvector[1] = -Bdx;
    bvector[2] = -Cdx;
    bvector[3] = -Ddx;
    bvector[4] = -Edx;
    bvector[5] = -Fdx;
    bvector[6] = -Gdx;
    bvector[7] = -Hdx;
    bvector[8] = Adxx;
    bvector[9] = Bdxx;
    bvector[10] = Cdxx;
    bvector[11] = Ddxx;
    bvector[12] = Edxx;
    bvector[13] = Fdxx;
    bvector[14] = Gdxx;
    bvector[15] = Hdxx;
    bvector[16] = Adxy;
    bvector[17] = Bdxy;
    bvector[18] = Cdxy;
    bvector[19] = Ddxy;
    bvector[20] = Edxy;
    bvector[21] = Fdxy;
    bvector[22] = Gdxy;
    bvector[23] = Hdxy;
    bvector[24] = Adxz;
    bvector[25] = Bdxz;
    bvector[26] = Cdxz;
    bvector[27] = Ddxz;
    bvector[28] = Edxz;
    bvector[29] = Fdxz;
    bvector[30] = Gdxz;
    bvector[31] = Hdxz;

    // printf("b = \n"); for (ii=0;ii<32;ii++) printf("%lf\n",bvector[ii]);

    // set the amatrix to the old values for next SVD 
    for (ii=0;ii<1024;ii++) amatrix[ii] = amatrix_backup[ii];

    // singular value decomposition
    // gsl_linalg_SV_decomp(&Umatrix.matrix, Vmatrix, Svector, work);
    gsl_linalg_SV_decomp_mod(&Umatrix.matrix, Xmatrix, Vmatrix, Svector, work);

    // printf("s = \n"); gsl_vector_fprintf(stdout,Svector,"%lf");

    // set the trivial values of svector to exactly 0.0 to avoild 
    // numerical accuracy problem. very important!
    for (ii=0;ii<32;ii++)
    {
	if (Svector->data[ii] <1.0e-8)
	    Svector->data[ii] = 0.0;
    }

    // solve the equation of Ax=b using singular vectors
    // minimizing least squares
    gsl_linalg_SV_solve(&Umatrix.matrix, Vmatrix, Svector, &b.vector, xvector);
    
    // calculate the value at the interpolatin position
    PP = 0.0;
    for (ii=0;ii<32;ii++)
	PP += interp_vector[ii]*xvector->data[ii];
    // assign the value to fx
    fsf[0] = PP;
    // printf("fx = %lf\n",PP);

    // --------------- fy ----------------
    // assign values to bvector for fx
    // ordered as A-H, value, dx, dy, dz
    // This is the only thing that needs to be changed
    // for energy and force interpolations
    bvector[0] = -Ady;
    bvector[1] = -Bdy;
    bvector[2] = -Cdy;
    bvector[3] = -Ddy;
    bvector[4] = -Edy;
    bvector[5] = -Fdy;
    bvector[6] = -Gdy;
    bvector[7] = -Hdy;
    bvector[8] = Adyx;
    bvector[9] = Bdyx;
    bvector[10] = Cdyx;
    bvector[11] = Ddyx;
    bvector[12] = Edyx;
    bvector[13] = Fdyx;
    bvector[14] = Gdyx;
    bvector[15] = Hdyx;
    bvector[16] = Adyy;
    bvector[17] = Bdyy;
    bvector[18] = Cdyy;
    bvector[19] = Ddyy;
    bvector[20] = Edyy;
    bvector[21] = Fdyy;
    bvector[22] = Gdyy;
    bvector[23] = Hdyy;
    bvector[24] = Adyz;
    bvector[25] = Bdyz;
    bvector[26] = Cdyz;
    bvector[27] = Ddyz;
    bvector[28] = Edyz;
    bvector[29] = Fdyz;
    bvector[30] = Gdyz;
    bvector[31] = Hdyz;

    // set the amatrix to the old values for next SVD 
    for (ii=0;ii<1024;ii++) amatrix[ii] = amatrix_backup[ii];

    // singular value decomposition
    // gsl_linalg_SV_decomp(&Umatrix.matrix, Vmatrix, Svector, work);
    gsl_linalg_SV_decomp_mod(&Umatrix.matrix, Xmatrix, Vmatrix, Svector, work);

    // set the trivial values of svector to exactly 0.0 to avoild 
    // numerical accuracy problem. very important!
    for (ii=0;ii<32;ii++)
    {
	if (Svector->data[ii] <1.0e-8)
	    Svector->data[ii] = 0.0;
    }

    // solve the equation of Ax=b using singular vectors
    // minimizing least squares
    gsl_linalg_SV_solve(&Umatrix.matrix, Vmatrix, Svector, &b.vector, xvector);
    
    // calculate the value at the interpolatin position
    PP = 0.0;
    for (ii=0;ii<32;ii++)
	PP += interp_vector[ii]*xvector->data[ii];
    // assign the value to fy
    fsf[1] = PP;
    // printf("fy = %lf\n",PP);

    // --------------- fz ----------------------
    // assign values to bvector for fx
    // ordered as A-H, value, dx, dy, dz
    // This is the only thing that needs to be changed
    // for energy and force interpolations
    bvector[0] = -Adz;
    bvector[1] = -Bdz;
    bvector[2] = -Cdz;
    bvector[3] = -Ddz;
    bvector[4] = -Edz;
    bvector[5] = -Fdz;
    bvector[6] = -Gdz;
    bvector[7] = -Hdz;
    bvector[8] = Adzx;
    bvector[9] = Bdzx;
    bvector[10] = Cdzx;
    bvector[11] = Ddzx;
    bvector[12] = Edzx;
    bvector[13] = Fdzx;
    bvector[14] = Gdzx;
    bvector[15] = Hdzx;
    bvector[16] = Adzy;
    bvector[17] = Bdzy;
    bvector[18] = Cdzy;
    bvector[19] = Ddzy;
    bvector[20] = Edzy;
    bvector[21] = Fdzy;
    bvector[22] = Gdzy;
    bvector[23] = Hdzy;
    bvector[24] = Adzz;
    bvector[25] = Bdzz;
    bvector[26] = Cdzz;
    bvector[27] = Ddzz;
    bvector[28] = Edzz;
    bvector[29] = Fdzz;
    bvector[30] = Gdzz;
    bvector[31] = Hdzz;

    // set the amatrix to the old values for next SVD 
    for (ii=0;ii<1024;ii++) amatrix[ii] = amatrix_backup[ii];

    // singular value decomposition
    // gsl_linalg_SV_decomp(&Umatrix.matrix, Vmatrix, Svector, work);
    gsl_linalg_SV_decomp_mod(&Umatrix.matrix, Xmatrix, Vmatrix, Svector, work);

    // set the trivial values of svector to exactly 0.0 to avoild 
    // numerical accuracy problem. very important!
    for (ii=0;ii<32;ii++)
    {
	if (Svector->data[ii] <1.0e-8)
	    Svector->data[ii] = 0.0;
    }

    // solve the equation of Ax=b using singular vectors
    // minimizing least squares
    gsl_linalg_SV_solve(&Umatrix.matrix, Vmatrix, Svector, &b.vector, xvector);
    
    // calculate the value at the interpolatin position
    PP = 0.0;
    for (ii=0;ii<32;ii++)
	PP += interp_vector[ii]*xvector->data[ii];
    // assign the value to fz
    fsf[2] = PP;
    // printf("fz = %lf\n",PP);



    printf("fx=%lf   fy=%lf   fz=%lf\n",fsf[0],fsf[1],fsf[2]);
    exit(1);
    */




    /*
       printf("x = \n");
       gsl_vector_fprintf(stdout,xvector,"%8.1f");

       printf("PP=%lf\n",PP);
       */








}

int init_sf_hypergeo()
{
    int ii;
    const int datalen = 200;
    char buffer[200];
    char keyword[100];

    fprintf(stderr,"Reading input data for hypergeometric nanotubes...\n");
    fprintf(fpouts,"Reading input data for hypergeometric nanotubes...\n");

    // re-open input file to read extra data section
    fpins = fopen(INPUT,"r");

    while (fgets(buffer,datalen,fpins)!=NULL)
    {
	sscanf(buffer,"%s",keyword);
	for (ii=0;ii<strlen(keyword);ii++)
	    keyword[ii] = toupper(keyword[ii]);
	if (!strcmp(keyword,"HYPERGEO"))
	{
	    fprintf(stderr,"Data section for hypergeometric nanotoubes found...\n");
	    fprintf(fpouts,"Data section for hypergeometric nanotoubes found...\n");
	    // use solid sigma and epsilon for hypergeometric parameters
	    // assume all the tubes have the same parameters
	    solid_sigma = malloc(sizeof(double));
	    solid_epsilon = malloc(sizeof(double));
	    sscanf(fgets(buffer,datalen,fpins),"%d %lf %lf",&ntube,solid_sigma,solid_epsilon);
	    fprintf(stderr,"ntube=%d  sigma=%lf  epsilon=%lf\n",ntube,*solid_sigma,*solid_epsilon);
	    fprintf(fpouts,"ntube=%d  sigma=%lf  epsilon=%lf\n",ntube,*solid_sigma,*solid_epsilon);
	    // allocate memories
	    hgntc_xx = calloc(ntube,sizeof(double));
	    hgntc_yy = calloc(ntube,sizeof(double));
	    hgnt_radius = calloc(ntube,sizeof(double));
	    for (ii=0;ii<ntube;ii++)
	    {
		sscanf(fgets(buffer,datalen,fpins),"%lf %lf %lf",&hgntc_xx[ii],&hgntc_yy[ii],&hgnt_radius[ii]);
		fprintf(stderr,"tube %d: xx=%lf  y=%lf  radius=%lf\n",ii,hgntc_xx[ii],hgntc_yy[ii],hgnt_radius[ii]);
		fprintf(fpouts,"tube %d: xx=%lf  y=%lf  radius=%lf\n",ii,hgntc_xx[ii],hgntc_yy[ii],hgnt_radius[ii]);
	    }
	    fclose(fpins);
	    return 0;
	} // if keyword found
    } // read through the lines
    fprintf(stderr,"Error: data for heypergeometric nanotubes not found.\n");
    fprintf(fpouts,"Error: data for heypergeometric nanotubes not found.\n");
    fclose(fpins);
    exit(1);
}

// calculate the hypergeo series
int hypergeo(double a, double b, double c,double z, double rr, double *sumV, double *sumF)
{
    double zn, an1, an;
    int n, i;

    *sumV = 1.0;
    zn = 1.0;
    an1 = 1.0;
    n = 10;
    *sumF = 0.0;

    for (i=0;i<=n;i++)
    {
	an = (a+i)*(b+i)/((c+i)*(1+i));
	an1 = an1*an;
	zn *= z;
	*sumV = *sumV + an1*zn;
	// sumF is for outside nanotube
	// for inside nanotube, sumF should be -sumF
	*sumF = *sumF - 2.0*(i+1.0)*an1*zn/rr;
    }

}

// calculate the hypergeometric nanotube potentials
// modified from Matt LaBrosse's code
int cal_sf_hypergeo()
{
    int ii, jj;
    double xxi, yyi;
    double fxi, fyi;
    double sigmaij, epsilonij;
    double rxij, ryij;
    double rijsq, rij;

    double tubeR, tubeRsq;
    double kfac, kfacsq, kfac4, kfac10;
    const double afac1 = -4.5, afac2 = -1.5; // -9.0/2.0 and -3.0/2.0
    const double bfac1 = -4.5, bfac2 = -1.5; // -9.0/2.0 and -3.0/2.0
    const double cfac = 1.0;
    double zfac;
    double hgrep, hgfrep;
    double rep1, rep2;
    const double c21_by_32 = 0.65625; // 21/32
    double hgattr, hgfattr;
    double attr1, attr2;
    double ljc;
    double uhyper;
    double fij, fxij, fyij;

    for (ii=0;ii<natom;ii++)
    {
	if (isghost[ii] == all_ghost) // ghost atom check
	    continue;
	xxi = xx[ii];
	yyi = yy[ii];
	fxi = fxl[ii];
	fyi = fyl[ii];
	sigmaij = 0.5*(sigma[ii]+(*solid_sigma));
	epsilonij = sqrt(epsilon[ii]*(*solid_epsilon));
	for (jj=0;jj<ntube;jj++)
	{
	    rxij = xxi - hgntc_xx[jj];
	    ryij = yyi - hgntc_yy[jj];
	    // minimum image convention
	    rxij = rxij - boxlx*rint(rxij/boxlx);
	    ryij = ryij - boxly*rint(ryij/boxly);
	    rijsq = rxij*rxij + ryij*ryij;
	    rij = sqrt(rijsq);
	    tubeR = hgnt_radius[jj];
	    tubeRsq = tubeR*tubeR;

	    if (rij<tubeR) // if inside the tube
	    {
		ljc = c3_pisq_theta*epsilonij*sigmaij*sigmaij;
		kfac = sigmaij*tubeR/(tubeRsq-rijsq);
		zfac = rijsq/tubeRsq;
		hypergeo(afac1,bfac1,cfac,zfac,rij,&hgrep,&hgfrep); // hypergeo series
		kfacsq = kfac*kfac; kfac4 = kfacsq*kfacsq; kfac10 = kfac4*kfac4*kfacsq;
		rep1 = c21_by_32*kfac10;
		rep2 = -hgfrep + hgrep*(20.0*rij/(tubeRsq-rijsq));
		hypergeo(afac2,bfac2,cfac,zfac,rij,&hgattr,&hgfattr); // calculate hypergeo series
		attr1 = -kfac4;
		attr2 = -hgfattr + hgattr*(8.0*rij/(tubeRsq-rijsq));
		uhyper = ljc*(rep1*hgrep + attr1*hgattr);
		fij = -ljc*(rep1*rep2 + attr1*attr2);
		fxij = fij*rxij/rij;
		fyij = fij*ryij/rij;
	    }
	    else // outside the tube
	    {
		ljc = c3_pisq_theta*epsilonij*sigmaij*sigmaij*tubeR;
		kfac = sigmaij*rij/(rijsq-tubeRsq);
		zfac = tubeRsq/rijsq;
		hypergeo(afac1,bfac1,cfac,zfac,rij,&hgrep,&hgfrep); // calculate the hypergeo series
		kfacsq = kfac*kfac; kfac4 = kfacsq*kfacsq; kfac10 = kfac4*kfac4*kfacsq;
		rep1 = c21_by_32*kfac10;
		rep2 = hgfrep/rij + hgrep*(-20.0/(rijsq-tubeRsq)+9.0/rijsq);
		hypergeo(afac2,bfac2,cfac,zfac,rij,&hgattr,&hgfattr); // calculate 
		attr1 = -kfac4;
		attr2 = hgfattr/rij + hgattr*(-8.0/(rijsq-tubeRsq)+3.0/rijsq);
		uhyper = ljc*(rep1*hgrep + attr1*hgattr)/rij;
		fij = -ljc*(rep1*rep2 + attr1*attr2);
		fxij = fij*rxij/rij;
		fyij = fij*ryij/rij;
	    } // outside tube
	    // energy and forces
	    usflj += uhyper;
	    fxi += fxij;
	    fyi += fyij;
	} // loop through all tubes
	fxl[ii] = fxi;
	fyl[ii] = fyi;
    } // loop through all fluid atoms
}

int init_sf_atom_explicit()
{
    int ii;
    const int datalen = 200;
    char buffer[200];
    char keyword[100];
    int itmp;

    fprintf(stderr,"Reading input data for atom explicit sorbents...\n");
    fprintf(fpouts,"Reading input data for atom explicit sorbents...\n");

    // re-open input file to read extra data section
    fpins = fopen(INPUT,"r");

    while (fgets(buffer,datalen,fpins)!=NULL)
    {
	sscanf(buffer,"%s",keyword);
	for (ii=0;ii<strlen(keyword);ii++)
	    keyword[ii] = toupper(keyword[ii]);
	if (!strcmp(keyword,"ATOMEXPLICIT"))
	{
	    fprintf(stderr,"Data section for atom explicit sorbents found...\n");
	    fprintf(fpouts,"Data section for atom explicit sorbents found...\n");
	    sscanf(fgets(buffer,datalen,fpins),"%d",&solid_natom);
	    sscanf(fgets(buffer,datalen,fpins),"%d",&fSolid_type);
	    // fSolid_type not yet in used
	    // if (fSolid_type==solid_uniform)
	    {
		solid_sigma = malloc(sizeof(double));
		solid_epsilon = malloc(sizeof(double));
		solid_charge = malloc(sizeof(double));
		// epsilon must be J/mol!! important
		sscanf(buffer,"%d %lf %lf %lf",&itmp,solid_sigma,solid_epsilon,solid_charge);
	    }
	    solid_xx = calloc(solid_natom,sizeof(double));
	    solid_yy = calloc(solid_natom,sizeof(double));
	    solid_zz = calloc(solid_natom,sizeof(double));
	    assert(solid_xx!=NULL);
	    assert(solid_yy!=NULL);
	    assert(solid_zz!=NULL);
	    // readin solid coordinates
	    for (ii=0;ii<solid_natom;ii++)
		fscanf(fpins,"%s %lf %lf %lf\n",buffer,&solid_xx[ii],&solid_yy[ii],&solid_zz[ii]);
	    fclose(fpins);
	    return 0;
	} // if keyword found
    } // read through the lines
    fprintf(stderr,"Error: data for atom explicit sorbents not found.\n");
    fprintf(fpouts,"Error: data for atom explicit sorbents not found.\n");
    fclose(fpins);
    exit(1);
}

// calcuate the interactions between solid and fluid
// using atom explicit model
int cal_sf_atom_explicit()
{
    int ii, jj;
    double xxi, yyi, zzi;
    double fxi, fyi, fzi;
    double rxij, ryij, rzij;
    double rijsq, rij, r_rijsq;
    double r_r6, r_r12, r_r12_minus_r_r6;
    double epsilonij, sigmaij;
    double usf_vdw_temp, usf_vdw;
    double fij, fxij, fyij, fzij;

    usf_vdw = 0.0;
    for (ii=0;ii<natom;ii++)
    {
	if (isghost[ii] == all_ghost) // do not calculate ghost atoms
	    continue;
	xxi = xx[ii];
	yyi = yy[ii];
	zzi = zz[ii];
	fxi = fxl[ii];
	fyi = fyl[ii];
	fzi = fzl[ii];
	// only uniform solid is considered now
	// so, calculate epsilonij and sigmaij here
	// Also assume solid has no charge
	sigmaij = 0.5*(sigma[ii]+(*solid_sigma));
	epsilonij = sqrt(epsilon[ii]*(*solid_epsilon));
	for (jj=0;jj<solid_natom;jj++)
	{
	    // assume no solid atom is ghost type
	    rxij = xxi - solid_xx[jj];
	    ryij = yyi - solid_yy[jj];
	    rzij = zzi - solid_zz[jj];
	    // minimum image convention
	    rxij = rxij - boxlx*rint(rxij/boxlx);
	    ryij = ryij - boxly*rint(ryij/boxly);
	    rzij = rzij - boxlz*rint(rzij/boxlz);
	    rijsq = rxij*rxij + ryij*ryij + rzij*rzij;
	    rij = sqrt(rijsq);

	    if (rijsq<rcutoffsq)
	    {
		if (isLJswitchOn)
		{
		    if (rijsq<rcutonsq)
			LJswitch = 1.0;
		    else
			LJswitch = (rcutoffsq-rijsq)*(rcutoffsq-rijsq)
			    *(rcutoffsq+2.0*rijsq-3.0*rcutonsq)/roff2_minus_ron2_cube;
		}
		r_rijsq = sigmaij*sigmaij/rijsq;
		r_r6 = r_rijsq*r_rijsq*r_rijsq;
		r_r12 = r_r6*r_r6;
		r_r12_minus_r_r6 = r_r12 - r_r6;
		usf_vdw_temp = epsilonij*r_r12_minus_r_r6; // still need *4.0
		if (isLJswitchOn) // if switch is used
		    usf_vdw += usf_vdw_temp*LJswitch; // still need 4.0
		else
		    usf_vdw += usf_vdw_temp; // still need *4.0
		// force calculation
		if (isLJswitchOn)
		{
		    if (rijsq<rcutonsq)
			fij = 24.0*epsilonij*(r_r12_minus_r_r6+r_r12)/rijsq;
		    else
			fij = 24.0*epsilonij*(r_r12_minus_r_r6+r_r12)/rijsq*LJswitch // 4.0 for the real energy
			    -4.0*usf_vdw_temp*12.0*(rcutoffsq-rijsq)*(rcutonsq-rijsq)/roff2_minus_ron2_cube;
		}
		else
		    fij = 24.0*epsilonij*(r_r12_minus_r_r6+r_r12)/rijsq;
		fxij = fij*rxij;
		fyij = fij*ryij;
		fzij = fij*rzij;
		// force on fluid atom ii
		fxi += fxij;
		fyi += fyij;
		fzi += fzij;
	    } // if with cutoff
	} // solid atom jj loop
	fxl[ii] = fxi;
	fyl[ii] = fyi;
	fzl[ii] = fzi;
    } // fluid atom ii loop

    // factor
    usflj += usf_vdw*4.0;
}

int sffrc()
{
    int ii;
    double usflj_tasos;
    double tasos_force[3];

    // reset energy
    usflj = 0.0;

    if (sf_type==nanotube_hypergeo)
    {
	cal_sf_hypergeo();
    }
    else if (sf_type==nanotube_atom_explicit)
    {
	cal_sf_atom_explicit();
    }
    else if (sf_type==nanotube_tasos)
    {
	for (ii=0;ii<natom;ii++)
	{
	    // call Tasos's code for force calculations
	    cforce_atom_(&tasostype[ii],&xx[ii],&yy[ii],&zz[ii],&usflj_tasos,tasos_force);

	    // the tasos energy has unit of K and the force has unit of K/Angstrom
	    // change them to J/mol and J/mol/Angstrom
	    usflj_tasos *= Rgas;
	    tasos_force[0] *= Rgas;
	    tasos_force[1] *= Rgas;
	    tasos_force[2] *= Rgas;

	    usflj += usflj_tasos;
	    fxl[ii] += tasos_force[0];
	    fyl[ii] += tasos_force[1];
	    fzl[ii] += tasos_force[2];
	} // natom loop
    }
    else if (sf_type==nanotube_my_interp)
    {
	for (ii=0;ii<natom;ii++)
	{
	    get_values_from_grid(xx[ii],yy[ii],zz[ii],tasostype[ii],&usflj_tasos,tasos_force);

	    usflj += usflj_tasos;

	    fxl[ii] += tasos_force[0];
	    fyl[ii] += tasos_force[1];
	    fzl[ii] += tasos_force[2];

	}
    }
}

