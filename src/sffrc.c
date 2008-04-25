#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <ctype.h>
#include "mspms2.h"

extern void cforce_atom_(int*, double*, double*, double*, double*, double*);
extern void initpotentialgrid_(int*, double*, double*, double*, int*, int*,
		int*, double*);
extern void pass_grid_file_name_(char*, int*);
extern void read_grids_(int* nspecies_yang);

/// Read and initialize parameters for hypergeo nanotubes.
int init_sf_hypergeo()
{
	int ii;
	char buffer[LONG_STRING_LENGTH];
	char keyword[100];

	fprintf(stderr,"Reading input data for hypergeometric nanotubes...\n");
	fprintf(fpouts, "Reading input data for hypergeometric nanotubes...\n");

	// re-open input file to read extra data section
	fpins = fopen(INPUT,"r");

	while (fgets(buffer, LONG_STRING_LENGTH, fpins)!=NULL)
	{
		sscanf(buffer, "%s", keyword);
		for (ii=0; ii<strlen(keyword); ii++)
		{
			keyword[ii] = toupper(keyword[ii]);
		}
		if (!strcmp(keyword, "HYPERGEO"))
		{
			fprintf(stderr,"Data section for hypergeometric nanotoubes found...\n");
			fprintf(fpouts, "Data section for hypergeometric nanotoubes found...\n");
			// use solid sigma and epsilon for hypergeometric parameters
			// assume all the tubes have the same parameters
			solid_sigma = malloc(sizeof(double));
			solid_epsilon = malloc(sizeof(double));
			sscanf(fgets(buffer, LONG_STRING_LENGTH, fpins), "%d", &ntube);
			sscanf(fgets(buffer, LONG_STRING_LENGTH, fpins), "%lf", solid_sigma);
			sscanf(fgets(buffer, LONG_STRING_LENGTH, fpins), "%lf", solid_epsilon);
			
			fprintf(stderr,"ntube=%d  sigma=%lf  epsilon=%lf\n",ntube,*solid_sigma,*solid_epsilon);
			fprintf(fpouts, "ntube=%d  sigma=%lf  epsilon=%lf\n", ntube,*solid_sigma, *solid_epsilon);
			
			// Reduce them
			*solid_sigma /= sigma_base;
			*solid_epsilon /= epsilon_base;
			
			// allocate memories
			hgntc_xx = calloc(ntube, sizeof(double));
			hgntc_yy = calloc(ntube, sizeof(double));
			hgnt_radius = calloc(ntube, sizeof(double));
			for (ii=0; ii<ntube; ii++)
			{
				sscanf(fgets(buffer, LONG_STRING_LENGTH, fpins), "%lf %lf %lf", &hgntc_xx[ii], &hgntc_yy[ii], &hgnt_radius[ii]);
				fprintf(stderr,"tube %d: xx=%lf  y=%lf  radius=%lf\n",ii,hgntc_xx[ii],hgntc_yy[ii],hgnt_radius[ii]);
				fprintf(fpouts, "tube %d: xx=%lf  y=%lf  radius=%lf\n", ii,hgntc_xx[ii], hgntc_yy[ii], hgnt_radius[ii]);
				// Reduce them
				hgntc_xx[ii] /= sigma_base;
				hgntc_yy[ii] /= sigma_base;
				hgnt_radius[ii] /= sigma_base;
			}
			// Reduce the constant
			const_3pisq_theta = C3_PI_SQ_THETA*sigma_base*sigma_base;
			fclose(fpins);
			return 0;
		} // if keyword found
	} // read through the lines
	fprintf(stderr,"Error: data for heypergeometric nanotubes not found.\n");
	fprintf(fpouts, "Error: data for heypergeometric nanotubes not found.\n");
	fclose(fpins);
	exit(1);
}

/// Read and initialize parameters for atom explicit adsorbents.
int init_sf_atom_explicit()
{
	int ii;
	char buffer[LONG_STRING_LENGTH];
	char keyword[100];
	int itmp;

	fprintf(stderr,"Reading input data for atom explicit sorbents...\n");
	fprintf(fpouts, "Reading input data for atom explicit sorbents...\n");

	// re-open input file to read extra data section
	fpins = fopen(INPUT,"r");

	while (fgets(buffer, LONG_STRING_LENGTH, fpins)!=NULL)
	{
		sscanf(buffer, "%s", keyword);
		for (ii=0; ii<strlen(keyword); ii++)
			keyword[ii] = toupper(keyword[ii]);
		if (!strcmp(keyword, "ATOMEXPLICIT"))
		{
			fprintf(stderr,"Data section for atom explicit sorbents found...\n");
			fprintf(fpouts,
					"Data section for atom explicit sorbents found...\n");
			sscanf(fgets(buffer, LONG_STRING_LENGTH, fpins), "%d", &solid_natom);
			sscanf(fgets(buffer, LONG_STRING_LENGTH, fpins), "%d", &fSolid_type);
			// fSolid_type not yet in used
			// if (fSolid_type==SOLID_UNIFORM)
			{
				solid_sigma = malloc(sizeof(double));
				solid_epsilon = malloc(sizeof(double));
				solid_charge = malloc(sizeof(double));
				// epsilon must be K!! important
				sscanf(buffer, "%d %lf %lf %lf", &itmp, solid_sigma,
						solid_epsilon, solid_charge);
				*solid_sigma /= sigma_base;
				*solid_epsilon /= epsilon_base;
			}
			solid_xx = calloc(solid_natom, sizeof(double));
			solid_yy = calloc(solid_natom, sizeof(double));
			solid_zz = calloc(solid_natom, sizeof(double));
			assert(solid_xx!=NULL);
			assert(solid_yy!=NULL);
			assert(solid_zz!=NULL);
			// readin solid coordinates
			for (ii=0; ii<solid_natom; ii++)
			{
				fscanf(fpins, "%s %lf %lf %lf\n", buffer, &solid_xx[ii],
						&solid_yy[ii], &solid_zz[ii]);
				// Reduce them
				solid_xx[ii] /= sigma_base;
				solid_yy[ii] /= sigma_base;
				solid_zz[ii] /= sigma_base;
			}
			fclose(fpins);
			return 0;
		} // if keyword found
	} // read through the lines
	fprintf(stderr,"Error: data for atom explicit sorbents not found.\n");
	fprintf(fpouts, "Error: data for atom explicit sorbents not found.\n");
	fclose(fpins);
	exit(1);
}

/** 
 * \brief Read and initialize and reading the tasos grids.
 * 
 * Note that this part does NOT use reduced units for historical reasons.
 */
int init_tasos_grid()
{
	int ii, jj;
	char buffer[LONG_STRING_LENGTH];
	double auc, buc, cuc, nanotuberadius;
	int na, nb, nc;
	int nspecies_yang;
	int tnuatoms;
	char szgrid[12][200]; // assuming only 12 types of grids
	char keyword[100];

	fprintf(stderr,"Reading input data for tasos grids...\n");
	fprintf(fpouts, "Reading input data for tasos grids...\n");

	// re-open input file to read extra data section
	fpins = fopen(INPUT,"r");

	while (fgets(buffer, LONG_STRING_LENGTH, fpins)!=NULL)
	{
		sscanf(buffer, "%s", keyword);
		for (ii=0; ii<strlen(keyword); ii++)
		{
			keyword[ii] = toupper(keyword[ii]);
		}
		if (!strcmp(keyword, "TASOS"))
		{
			fprintf(stderr,"Data section for tasos grids found...\n");
			fprintf(fpouts, "Data section for tasos grids found...\n");

			sscanf(fgets(buffer, LONG_STRING_LENGTH, fpins), "%lf %lf %lf %lf",
					&auc, &buc, &cuc, &nanotuberadius);
			sscanf(fgets(buffer, LONG_STRING_LENGTH, fpins), "%d %d %d", &na,
					&nb, &nc);
			sscanf(fgets(buffer, LONG_STRING_LENGTH, fpins), "%d",
					&nspecies_yang);
			sscanf(fgets(buffer, LONG_STRING_LENGTH, fpins), "%d", &tnuatoms);
			fprintf(
					fpouts,
					"The size of the unit cell is %lf long, %lf wide, %lf high, with %lf radius.\n",
					auc, buc, cuc, nanotuberadius);
			fprintf(fpouts,
					"There are %d unique species and %d unique atoms.\n",
					nspecies_yang, tnuatoms);
			for (ii=0; ii<tnuatoms; ii++)
			{
				sscanf(fgets(buffer, LONG_STRING_LENGTH, fpins), "%s", szgrid
						+ii);
			}
			initpotentialgrid_(&tnuatoms, &auc, &buc, &cuc, &na, &nb, &nc,
					&nanotuberadius);
			fprintf(stderr,"init potential ok\n");
			fprintf(fpouts, "init potential ok\n");
			for (ii=0; ii<tnuatoms; ii++)
			{
				jj = ii + 1;
				pass_grid_file_name_(szgrid[ii], &jj);
				fprintf(stderr,"%s\n",szgrid[ii]);
				fprintf(fpouts, "%s\n", szgrid[ii]);
			}
			fprintf(stderr,"pass grid file name ok\n");
			fprintf(fpouts, "pass grid file name ok\n");
			read_grids_(&nspecies_yang);
			fprintf(stderr,"read grids ok\n");
			fprintf(fpouts, "read grids ok\n");

			fclose(fpins);
			return 0;
		} // if keyword found
	} // read through the lines
	fprintf(stderr,"Error: data for tasos grids not found.\n");
	fprintf(fpouts, "Error: data for tasos grids not found.\n");
	fclose(fpins);
	exit(1);
}

/**
 * Read and initialize the parameters for my interpolation.
 * Note that this part does NOT use reduced units for historical reasons.
 */
int init_my_interp()
{
	int ii, jj;
	PSAMPLE_MOLECULE pSampleMole;
	FILE *fpgridfile;
	char buffer[LONG_STRING_LENGTH];
	int nunique_atom;
	char szgrid[200];
	char keyword[100];
	double uclx_chk, ucly_chk, uclz_chk;
	double temp;

	fprintf(stderr,"Warning: this solid-fluid interpolation method is not fully tested.\n");
	fprintf(fpouts,
			"Warning: this solid-fluid interpolation method is not fully tested.\n");
	fprintf(stderr,"Reading input data for myinterp...\n");
	fprintf(fpouts, "Reading input data for myinterp...\n");

	// re-open input file to read extra data section
	fpins = fopen(INPUT,"r");

	// Change the atom type from tasos type to myinterp type.
	// myinterp type = tasos type - 1
	for (ii=0;ii<nspecie;ii++)
	{
		pSampleMole = sample_mole + ii;
		for (jj=0;jj<pSampleMole->natom;jj++)
		{
			pSampleMole->interp_type[jj] -= 1;
		}
	}

	// set up the variable needed for interpolation
	_safealloc(interp_vector,32,sizeof(double))
	;

	while (fgets(buffer, LONG_STRING_LENGTH, fpins)!=NULL)
	{
		sscanf(buffer, "%s", keyword);
		for (ii=0; ii<strlen(keyword); ii++)
		{
			keyword[ii] = toupper(keyword[ii]);
		}
		if (!strcmp(keyword, "MYINTERP"))
		{
			fprintf(stderr,"Data section for myinterp found...\n");
			fprintf(fpouts, "Data section for myinterp found...\n");

			sscanf(fgets(buffer, LONG_STRING_LENGTH, fpins), "%lf %lf %lf", &uclx, &ucly, &uclz);
			sscanf(fgets(buffer, LONG_STRING_LENGTH, fpins), "%lf %lf %lf", &xcenter, &ycenter, &zcenter);
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
			fprintf(fpouts, "xmin=%lf  xmax=%lf\n", xmin, xmax);
			fprintf(fpouts, "ymin=%lf  ymax=%lf\n", ymin, ymax);
			fprintf(fpouts, "zmin=%lf  zmax=%lf\n", zmin, zmax);
			sscanf(fgets(buffer, LONG_STRING_LENGTH, fpins), "%d",
					&nunique_atom);
			for (ii=0; ii<nunique_atom; ii++)
			{
				sscanf(fgets(buffer, LONG_STRING_LENGTH, fpins), "%s", szgrid);
				fpgridfile = fopen(szgrid, "rb");
				fread(&uclx_chk, sizeof(double), 1, fpgridfile);
				fread(&ucly_chk, sizeof(double), 1, fpgridfile);
				fread(&uclz_chk, sizeof(double), 1, fpgridfile);
				if (uclx==uclx_chk && ucly==ucly_chk && uclz==uclz_chk)
				{
					fprintf(stderr,"unit cell size matches for grid: %s\n",szgrid);
					fprintf(fpouts, "unit cell size matches for grid: %s\n",
							szgrid);
				}
				else
				{
					fprintf(stderr,
					"Error: unit cell mismatch: %s\nuclx=%lf/%lf  ucly=%lf/%lf  uclz=%lf/%lf\n",
					szgrid, uclx,uclx_chk,ucly,ucly_chk,uclz,uclz_chk);
					fprintf(
							fpouts,
							"Error: unit cell mismatch: %s\nuclx=%lf/%lf  ucly=%lf/%lf  uclz=%lf/%lf\n",
							szgrid, uclx, uclx_chk, ucly, ucly_chk, uclz,
							uclz_chk);
					exit(1);
				}
				// read in the center coordinates
				// not really useful here
				fread(&temp, sizeof(double), 1, fpgridfile);
				fread(&temp, sizeof(double), 1, fpgridfile);
				fread(&temp, sizeof(double), 1, fpgridfile);
				fread(&ngrid_total, sizeof(int), 1, fpgridfile);
				fread(&ngrid_x, sizeof(int), 1, fpgridfile);
				fread(&ngrid_y, sizeof(int), 1, fpgridfile);
				fread(&ngrid_z, sizeof(int), 1, fpgridfile);
				fprintf(stderr,"%d=%d*%d*%d grids in grid file %d.\n",
				ngrid_total,ngrid_x,ngrid_y,ngrid_z,ii);
				fprintf(fpouts, "%d=%d*%d*%d grids in grid file %d.\n",
						ngrid_total, ngrid_x, ngrid_y, ngrid_z, ii);
				ncube_x = ngrid_x - 1;
				ncube_y = ngrid_y - 1;
				ncube_z = ngrid_z - 1;
				ncube_total = ncube_x*ncube_y*ncube_z;
				fprintf(stderr,"contains %d=%d*%d*%d cubes.\n",
				ncube_total,ncube_x,ncube_y,ncube_z);
				fprintf(fpouts, "contains %d=%d*%d*%d cubes.\n", ncube_total,
						ncube_x, ncube_y, ncube_z);
				fread(&grid_itvl_x, sizeof(double), 1, fpgridfile);
				fread(&grid_itvl_y, sizeof(double), 1, fpgridfile);
				fread(&grid_itvl_z, sizeof(double), 1, fpgridfile);
				fprintf(stderr,"grid size x,y,z = %lf,%lf,%lf\n",grid_itvl_x,grid_itvl_y,grid_itvl_z);
				fprintf(fpouts, "grid size x,y,z = %lf,%lf,%lf\n", grid_itvl_x,
						grid_itvl_y, grid_itvl_z);

				_safealloc(ene0[ii],ncube_total,sizeof(double))
				;
				_safealloc(ene1[ii],ncube_total,sizeof(double))
				;
				_safealloc(ene2[ii],ncube_total,sizeof(double))
				;
				_safealloc(ene3[ii],ncube_total,sizeof(double))
				;
				_safealloc(ene4[ii],ncube_total,sizeof(double))
				;
				_safealloc(ene5[ii],ncube_total,sizeof(double))
				;
				_safealloc(ene6[ii],ncube_total,sizeof(double))
				;
				_safealloc(ene7[ii],ncube_total,sizeof(double))
				;
				_safealloc(ene8[ii],ncube_total,sizeof(double))
				;
				_safealloc(ene9[ii],ncube_total,sizeof(double))
				;
				_safealloc(ene10[ii],ncube_total,sizeof(double))
				;
				_safealloc(ene11[ii],ncube_total,sizeof(double))
				;
				_safealloc(ene12[ii],ncube_total,sizeof(double))
				;
				_safealloc(ene13[ii],ncube_total,sizeof(double))
				;
				_safealloc(ene14[ii],ncube_total,sizeof(double))
				;
				_safealloc(ene15[ii],ncube_total,sizeof(double))
				;
				_safealloc(ene16[ii],ncube_total,sizeof(double))
				;
				_safealloc(ene17[ii],ncube_total,sizeof(double))
				;
				_safealloc(ene18[ii],ncube_total,sizeof(double))
				;
				_safealloc(ene19[ii],ncube_total,sizeof(double))
				;
				_safealloc(ene20[ii],ncube_total,sizeof(double))
				;
				_safealloc(ene21[ii],ncube_total,sizeof(double))
				;
				_safealloc(ene22[ii],ncube_total,sizeof(double))
				;
				_safealloc(ene23[ii],ncube_total,sizeof(double))
				;
				_safealloc(ene24[ii],ncube_total,sizeof(double))
				;
				_safealloc(ene25[ii],ncube_total,sizeof(double))
				;
				_safealloc(ene26[ii],ncube_total,sizeof(double))
				;
				_safealloc(ene27[ii],ncube_total,sizeof(double))
				;
				_safealloc(ene28[ii],ncube_total,sizeof(double))
				;
				_safealloc(ene29[ii],ncube_total,sizeof(double))
				;
				_safealloc(ene30[ii],ncube_total,sizeof(double))
				;
				_safealloc(ene31[ii],ncube_total,sizeof(double))
				;

				_safealloc(fxa0[ii],ncube_total,sizeof(double))
				;
				_safealloc(fxa1[ii],ncube_total,sizeof(double))
				;
				_safealloc(fxa2[ii],ncube_total,sizeof(double))
				;
				_safealloc(fxa3[ii],ncube_total,sizeof(double))
				;
				_safealloc(fxa4[ii],ncube_total,sizeof(double))
				;
				_safealloc(fxa5[ii],ncube_total,sizeof(double))
				;
				_safealloc(fxa6[ii],ncube_total,sizeof(double))
				;
				_safealloc(fxa7[ii],ncube_total,sizeof(double))
				;
				_safealloc(fxa8[ii],ncube_total,sizeof(double))
				;
				_safealloc(fxa9[ii],ncube_total,sizeof(double))
				;
				_safealloc(fxa10[ii],ncube_total,sizeof(double))
				;
				_safealloc(fxa11[ii],ncube_total,sizeof(double))
				;
				_safealloc(fxa12[ii],ncube_total,sizeof(double))
				;
				_safealloc(fxa13[ii],ncube_total,sizeof(double))
				;
				_safealloc(fxa14[ii],ncube_total,sizeof(double))
				;
				_safealloc(fxa15[ii],ncube_total,sizeof(double))
				;
				_safealloc(fxa16[ii],ncube_total,sizeof(double))
				;
				_safealloc(fxa17[ii],ncube_total,sizeof(double))
				;
				_safealloc(fxa18[ii],ncube_total,sizeof(double))
				;
				_safealloc(fxa19[ii],ncube_total,sizeof(double))
				;
				_safealloc(fxa20[ii],ncube_total,sizeof(double))
				;
				_safealloc(fxa21[ii],ncube_total,sizeof(double))
				;
				_safealloc(fxa22[ii],ncube_total,sizeof(double))
				;
				_safealloc(fxa23[ii],ncube_total,sizeof(double))
				;
				_safealloc(fxa24[ii],ncube_total,sizeof(double))
				;
				_safealloc(fxa25[ii],ncube_total,sizeof(double))
				;
				_safealloc(fxa26[ii],ncube_total,sizeof(double))
				;
				_safealloc(fxa27[ii],ncube_total,sizeof(double))
				;
				_safealloc(fxa28[ii],ncube_total,sizeof(double))
				;
				_safealloc(fxa29[ii],ncube_total,sizeof(double))
				;
				_safealloc(fxa30[ii],ncube_total,sizeof(double))
				;
				_safealloc(fxa31[ii],ncube_total,sizeof(double))
				;

				_safealloc(fya0[ii],ncube_total,sizeof(double))
				;
				_safealloc(fya1[ii],ncube_total,sizeof(double))
				;
				_safealloc(fya2[ii],ncube_total,sizeof(double))
				;
				_safealloc(fya3[ii],ncube_total,sizeof(double))
				;
				_safealloc(fya4[ii],ncube_total,sizeof(double))
				;
				_safealloc(fya5[ii],ncube_total,sizeof(double))
				;
				_safealloc(fya6[ii],ncube_total,sizeof(double))
				;
				_safealloc(fya7[ii],ncube_total,sizeof(double))
				;
				_safealloc(fya8[ii],ncube_total,sizeof(double))
				;
				_safealloc(fya9[ii],ncube_total,sizeof(double))
				;
				_safealloc(fya10[ii],ncube_total,sizeof(double))
				;
				_safealloc(fya11[ii],ncube_total,sizeof(double))
				;
				_safealloc(fya12[ii],ncube_total,sizeof(double))
				;
				_safealloc(fya13[ii],ncube_total,sizeof(double))
				;
				_safealloc(fya14[ii],ncube_total,sizeof(double))
				;
				_safealloc(fya15[ii],ncube_total,sizeof(double))
				;
				_safealloc(fya16[ii],ncube_total,sizeof(double))
				;
				_safealloc(fya17[ii],ncube_total,sizeof(double))
				;
				_safealloc(fya18[ii],ncube_total,sizeof(double))
				;
				_safealloc(fya19[ii],ncube_total,sizeof(double))
				;
				_safealloc(fya20[ii],ncube_total,sizeof(double))
				;
				_safealloc(fya21[ii],ncube_total,sizeof(double))
				;
				_safealloc(fya22[ii],ncube_total,sizeof(double))
				;
				_safealloc(fya23[ii],ncube_total,sizeof(double))
				;
				_safealloc(fya24[ii],ncube_total,sizeof(double))
				;
				_safealloc(fya25[ii],ncube_total,sizeof(double))
				;
				_safealloc(fya26[ii],ncube_total,sizeof(double))
				;
				_safealloc(fya27[ii],ncube_total,sizeof(double))
				;
				_safealloc(fya28[ii],ncube_total,sizeof(double))
				;
				_safealloc(fya29[ii],ncube_total,sizeof(double))
				;
				_safealloc(fya30[ii],ncube_total,sizeof(double))
				;
				_safealloc(fya31[ii],ncube_total,sizeof(double))
				;

				_safealloc(fza0[ii],ncube_total,sizeof(double))
				;
				_safealloc(fza1[ii],ncube_total,sizeof(double))
				;
				_safealloc(fza2[ii],ncube_total,sizeof(double))
				;
				_safealloc(fza3[ii],ncube_total,sizeof(double))
				;
				_safealloc(fza4[ii],ncube_total,sizeof(double))
				;
				_safealloc(fza5[ii],ncube_total,sizeof(double))
				;
				_safealloc(fza6[ii],ncube_total,sizeof(double))
				;
				_safealloc(fza7[ii],ncube_total,sizeof(double))
				;
				_safealloc(fza8[ii],ncube_total,sizeof(double))
				;
				_safealloc(fza9[ii],ncube_total,sizeof(double))
				;
				_safealloc(fza10[ii],ncube_total,sizeof(double))
				;
				_safealloc(fza11[ii],ncube_total,sizeof(double))
				;
				_safealloc(fza12[ii],ncube_total,sizeof(double))
				;
				_safealloc(fza13[ii],ncube_total,sizeof(double))
				;
				_safealloc(fza14[ii],ncube_total,sizeof(double))
				;
				_safealloc(fza15[ii],ncube_total,sizeof(double))
				;
				_safealloc(fza16[ii],ncube_total,sizeof(double))
				;
				_safealloc(fza17[ii],ncube_total,sizeof(double))
				;
				_safealloc(fza18[ii],ncube_total,sizeof(double))
				;
				_safealloc(fza19[ii],ncube_total,sizeof(double))
				;
				_safealloc(fza20[ii],ncube_total,sizeof(double))
				;
				_safealloc(fza21[ii],ncube_total,sizeof(double))
				;
				_safealloc(fza22[ii],ncube_total,sizeof(double))
				;
				_safealloc(fza23[ii],ncube_total,sizeof(double))
				;
				_safealloc(fza24[ii],ncube_total,sizeof(double))
				;
				_safealloc(fza25[ii],ncube_total,sizeof(double))
				;
				_safealloc(fza26[ii],ncube_total,sizeof(double))
				;
				_safealloc(fza27[ii],ncube_total,sizeof(double))
				;
				_safealloc(fza28[ii],ncube_total,sizeof(double))
				;
				_safealloc(fza29[ii],ncube_total,sizeof(double))
				;
				_safealloc(fza30[ii],ncube_total,sizeof(double))
				;
				_safealloc(fza31[ii],ncube_total,sizeof(double))
				;

				fread(ene0[ii], sizeof(double), ncube_total, fpgridfile);
				fread(ene1[ii], sizeof(double), ncube_total, fpgridfile);
				fread(ene2[ii], sizeof(double), ncube_total, fpgridfile);
				fread(ene3[ii], sizeof(double), ncube_total, fpgridfile);
				fread(ene4[ii], sizeof(double), ncube_total, fpgridfile);
				fread(ene5[ii], sizeof(double), ncube_total, fpgridfile);
				fread(ene6[ii], sizeof(double), ncube_total, fpgridfile);
				fread(ene7[ii], sizeof(double), ncube_total, fpgridfile);
				fread(ene8[ii], sizeof(double), ncube_total, fpgridfile);
				fread(ene9[ii], sizeof(double), ncube_total, fpgridfile);
				fread(ene10[ii], sizeof(double), ncube_total, fpgridfile);
				fread(ene11[ii], sizeof(double), ncube_total, fpgridfile);
				fread(ene12[ii], sizeof(double), ncube_total, fpgridfile);
				fread(ene13[ii], sizeof(double), ncube_total, fpgridfile);
				fread(ene14[ii], sizeof(double), ncube_total, fpgridfile);
				fread(ene15[ii], sizeof(double), ncube_total, fpgridfile);
				fread(ene16[ii], sizeof(double), ncube_total, fpgridfile);
				fread(ene17[ii], sizeof(double), ncube_total, fpgridfile);
				fread(ene18[ii], sizeof(double), ncube_total, fpgridfile);
				fread(ene19[ii], sizeof(double), ncube_total, fpgridfile);
				fread(ene20[ii], sizeof(double), ncube_total, fpgridfile);
				fread(ene21[ii], sizeof(double), ncube_total, fpgridfile);
				fread(ene22[ii], sizeof(double), ncube_total, fpgridfile);
				fread(ene23[ii], sizeof(double), ncube_total, fpgridfile);
				fread(ene24[ii], sizeof(double), ncube_total, fpgridfile);
				fread(ene25[ii], sizeof(double), ncube_total, fpgridfile);
				fread(ene26[ii], sizeof(double), ncube_total, fpgridfile);
				fread(ene27[ii], sizeof(double), ncube_total, fpgridfile);
				fread(ene28[ii], sizeof(double), ncube_total, fpgridfile);
				fread(ene29[ii], sizeof(double), ncube_total, fpgridfile);
				fread(ene30[ii], sizeof(double), ncube_total, fpgridfile);
				fread(ene31[ii], sizeof(double), ncube_total, fpgridfile);

				fread(fxa0[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fxa1[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fxa2[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fxa3[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fxa4[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fxa5[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fxa6[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fxa7[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fxa8[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fxa9[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fxa10[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fxa11[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fxa12[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fxa13[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fxa14[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fxa15[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fxa16[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fxa17[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fxa18[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fxa19[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fxa20[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fxa21[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fxa22[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fxa23[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fxa24[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fxa25[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fxa26[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fxa27[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fxa28[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fxa29[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fxa30[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fxa31[ii], sizeof(double), ncube_total, fpgridfile);

				fread(fya0[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fya1[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fya2[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fya3[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fya4[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fya5[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fya6[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fya7[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fya8[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fya9[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fya10[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fya11[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fya12[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fya13[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fya14[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fya15[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fya16[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fya17[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fya18[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fya19[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fya20[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fya21[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fya22[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fya23[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fya24[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fya25[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fya26[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fya27[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fya28[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fya29[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fya30[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fya31[ii], sizeof(double), ncube_total, fpgridfile);

				fread(fza0[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fza1[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fza2[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fza3[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fza4[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fza5[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fza6[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fza7[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fza8[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fza9[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fza10[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fza11[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fza12[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fza13[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fza14[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fza15[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fza16[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fza17[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fza18[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fza19[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fza20[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fza21[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fza22[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fza23[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fza24[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fza25[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fza26[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fza27[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fza28[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fza29[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fza30[ii], sizeof(double), ncube_total, fpgridfile);
				fread(fza31[ii], sizeof(double), ncube_total, fpgridfile);

				fclose(fpgridfile);
			} // loop through unique atoms grid files
			fclose(fpins);
			return 0;
		} // if key word found
	} // read through lines
	fprintf(stderr,"Error: data for myinterp not found.\n");
	fprintf(fpouts, "Error: data for myinterp not found.\n");
	fclose(fpins);
	exit(1);
}

int Amatrix_ele_assign_value(double *line, double x, double y, double z)
{
	double x2, y2, z2;
	double x3, y3, z3;

	x2=x*x;
	x3=x2*x;
	y2=y*y;
	y3=y2*y;
	z2=z*z;
	z3=z2*z;

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

	return 0;
}

int Amatrix_ele_assign_dx(double *line, double x, double y, double z)
{
	double x2, y2, z2;
	double x3, y3, z3;

	x2=x*x;
	x3=x2*x;
	y2=y*y;
	y3=y2*y;
	z2=z*z;
	z3=z2*z;

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

	return 0;
}

int Amatrix_ele_assign_dy(double *line, double x, double y, double z)
{
	double x2, y2, z2;
	double x3, y3, z3;

	x2=x*x;
	x3=x2*x;
	y2=y*y;
	y3=y2*y;
	z2=z*z;
	z3=z2*z;

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

	return 0;
}

int Amatrix_ele_assign_dz(double *line, double x, double y, double z)
{
	double x2, y2, z2;
	double x3, y3, z3;

	x2=x*x;
	x3=x2*x;
	y2=y*y;
	y3=y2*y;
	z2=z*z;
	z3=z2*z;

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

	return 0;
}

/**
 * This function accepts reduced parameters. However, they are converted
 * to real units inside the function and the output has units of J/mol
 * and J/mol/Angstrom.
 */
int get_values_from_grid(double fxx, double fyy, double fzz, int type,
		double *usf, double *fsf)
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
	// index of the lowest corner A
	int iAx, iAy, iAz;
	// index of the interpolation cube
	int cubeidx;
	
	fxx *= sigma_base;
	fyy *= sigma_base;
	fzz *= sigma_base;

	// calculate the index of corner A
	iAx = (int)((fxx-xmin)/grid_itvl_x);
	iAy = (int)((fyy-ymin)/grid_itvl_y);
	// periodical for z direction
	fzz = fzz - zmax*rint(fzz/zmax);
	// make sure its in the primal box
	if (fzz<0)
		fzz=zmax+fzz;
	iAz = (int)((fzz-zmin)/grid_itvl_z);

	// calcualte the position of the interpolation cube
	cubeidx = iAx*ncube_y*ncube_z + iAy*ncube_z + iAz;
	assert(cubeidx<ncube_total);

	// assign the value of the interpolation vector at x,y,z
	Amatrix_ele_assign_value(interp_vector, fxx, fyy, fzz);

	// energy interpolation
	*usf = ene0[type][cubeidx]*interp_vector[0] + ene1[type][cubeidx]
			*interp_vector[1] + ene2[type][cubeidx]*interp_vector[2]
			+ ene3[type][cubeidx]*interp_vector[3] + ene4[type][cubeidx]
			*interp_vector[4] + ene5[type][cubeidx]*interp_vector[5]
			+ ene6[type][cubeidx]*interp_vector[6] + ene7[type][cubeidx]
			*interp_vector[7] + ene8[type][cubeidx]*interp_vector[8]
			+ ene9[type][cubeidx]*interp_vector[9] + ene10[type][cubeidx]
			*interp_vector[10] + ene11[type][cubeidx]*interp_vector[11]
			+ ene12[type][cubeidx]*interp_vector[12] + ene13[type][cubeidx]
			*interp_vector[13] + ene14[type][cubeidx]*interp_vector[14]
			+ ene15[type][cubeidx]*interp_vector[15] + ene16[type][cubeidx]
			*interp_vector[16] + ene17[type][cubeidx]*interp_vector[17]
			+ ene18[type][cubeidx]*interp_vector[18] + ene19[type][cubeidx]
			*interp_vector[19] + ene20[type][cubeidx]*interp_vector[20]
			+ ene21[type][cubeidx]*interp_vector[21] + ene22[type][cubeidx]
			*interp_vector[22] + ene23[type][cubeidx]*interp_vector[23]
			+ ene24[type][cubeidx]*interp_vector[24] + ene25[type][cubeidx]
			*interp_vector[25] + ene26[type][cubeidx]*interp_vector[26]
			+ ene27[type][cubeidx]*interp_vector[27] + ene28[type][cubeidx]
			*interp_vector[28] + ene29[type][cubeidx]*interp_vector[29]
			+ ene30[type][cubeidx]*interp_vector[30] + ene31[type][cubeidx]
			*interp_vector[31];

	// fx
	fsf[0] = fxa0[type][cubeidx]*interp_vector[0] + fxa1[type][cubeidx]
			*interp_vector[1] + fxa2[type][cubeidx]*interp_vector[2]
			+ fxa3[type][cubeidx]*interp_vector[3] + fxa4[type][cubeidx]
			*interp_vector[4] + fxa5[type][cubeidx]*interp_vector[5]
			+ fxa6[type][cubeidx]*interp_vector[6] + fxa7[type][cubeidx]
			*interp_vector[7] + fxa8[type][cubeidx]*interp_vector[8]
			+ fxa9[type][cubeidx]*interp_vector[9] + fxa10[type][cubeidx]
			*interp_vector[10] + fxa11[type][cubeidx]*interp_vector[11]
			+ fxa12[type][cubeidx]*interp_vector[12] + fxa13[type][cubeidx]
			*interp_vector[13] + fxa14[type][cubeidx]*interp_vector[14]
			+ fxa15[type][cubeidx]*interp_vector[15] + fxa16[type][cubeidx]
			*interp_vector[16] + fxa17[type][cubeidx]*interp_vector[17]
			+ fxa18[type][cubeidx]*interp_vector[18] + fxa19[type][cubeidx]
			*interp_vector[19] + fxa20[type][cubeidx]*interp_vector[20]
			+ fxa21[type][cubeidx]*interp_vector[21] + fxa22[type][cubeidx]
			*interp_vector[22] + fxa23[type][cubeidx]*interp_vector[23]
			+ fxa24[type][cubeidx]*interp_vector[24] + fxa25[type][cubeidx]
			*interp_vector[25] + fxa26[type][cubeidx]*interp_vector[26]
			+ fxa27[type][cubeidx]*interp_vector[27] + fxa28[type][cubeidx]
			*interp_vector[28] + fxa29[type][cubeidx]*interp_vector[29]
			+ fxa30[type][cubeidx]*interp_vector[30] + fxa31[type][cubeidx]
			*interp_vector[31];

	// fy
	fsf[1] = fya0[type][cubeidx]*interp_vector[0] + fya1[type][cubeidx]
			*interp_vector[1] + fya2[type][cubeidx]*interp_vector[2]
			+ fya3[type][cubeidx]*interp_vector[3] + fya4[type][cubeidx]
			*interp_vector[4] + fya5[type][cubeidx]*interp_vector[5]
			+ fya6[type][cubeidx]*interp_vector[6] + fya7[type][cubeidx]
			*interp_vector[7] + fya8[type][cubeidx]*interp_vector[8]
			+ fya9[type][cubeidx]*interp_vector[9] + fya10[type][cubeidx]
			*interp_vector[10] + fya11[type][cubeidx]*interp_vector[11]
			+ fya12[type][cubeidx]*interp_vector[12] + fya13[type][cubeidx]
			*interp_vector[13] + fya14[type][cubeidx]*interp_vector[14]
			+ fya15[type][cubeidx]*interp_vector[15] + fya16[type][cubeidx]
			*interp_vector[16] + fya17[type][cubeidx]*interp_vector[17]
			+ fya18[type][cubeidx]*interp_vector[18] + fya19[type][cubeidx]
			*interp_vector[19] + fya20[type][cubeidx]*interp_vector[20]
			+ fya21[type][cubeidx]*interp_vector[21] + fya22[type][cubeidx]
			*interp_vector[22] + fya23[type][cubeidx]*interp_vector[23]
			+ fya24[type][cubeidx]*interp_vector[24] + fya25[type][cubeidx]
			*interp_vector[25] + fya26[type][cubeidx]*interp_vector[26]
			+ fya27[type][cubeidx]*interp_vector[27] + fya28[type][cubeidx]
			*interp_vector[28] + fya29[type][cubeidx]*interp_vector[29]
			+ fya30[type][cubeidx]*interp_vector[30] + fya31[type][cubeidx]
			*interp_vector[31];

	// fz
	fsf[2] = fza0[type][cubeidx]*interp_vector[0] + fza1[type][cubeidx]
			*interp_vector[1] + fza2[type][cubeidx]*interp_vector[2]
			+ fza3[type][cubeidx]*interp_vector[3] + fza4[type][cubeidx]
			*interp_vector[4] + fza5[type][cubeidx]*interp_vector[5]
			+ fza6[type][cubeidx]*interp_vector[6] + fza7[type][cubeidx]
			*interp_vector[7] + fza8[type][cubeidx]*interp_vector[8]
			+ fza9[type][cubeidx]*interp_vector[9] + fza10[type][cubeidx]
			*interp_vector[10] + fza11[type][cubeidx]*interp_vector[11]
			+ fza12[type][cubeidx]*interp_vector[12] + fza13[type][cubeidx]
			*interp_vector[13] + fza14[type][cubeidx]*interp_vector[14]
			+ fza15[type][cubeidx]*interp_vector[15] + fza16[type][cubeidx]
			*interp_vector[16] + fza17[type][cubeidx]*interp_vector[17]
			+ fza18[type][cubeidx]*interp_vector[18] + fza19[type][cubeidx]
			*interp_vector[19] + fza20[type][cubeidx]*interp_vector[20]
			+ fza21[type][cubeidx]*interp_vector[21] + fza22[type][cubeidx]
			*interp_vector[22] + fza23[type][cubeidx]*interp_vector[23]
			+ fza24[type][cubeidx]*interp_vector[24] + fza25[type][cubeidx]
			*interp_vector[25] + fza26[type][cubeidx]*interp_vector[26]
			+ fza27[type][cubeidx]*interp_vector[27] + fza28[type][cubeidx]
			*interp_vector[28] + fza29[type][cubeidx]*interp_vector[29]
			+ fza30[type][cubeidx]*interp_vector[30] + fza31[type][cubeidx]
			*interp_vector[31];

	return 0;
}

// calculate the hypergeo series
int hypergeo(double a, double b, double c, double z, double rr, double *sumV,
		double *sumF)
{
	double zn, an1, an;
	int n, i;

	*sumV = 1.0;
	zn = 1.0;
	an1 = 1.0;
	n = 10;
	*sumF = 0.0;

	for (i=0; i<=n; i++)
	{
		an = (a+i)*(b+i)/((c+i)*(1+i));
		an1 = an1*an;
		zn *= z;
		*sumV = *sumV + an1*zn;
		// sumF is for outside nanotube
		// for inside nanotube, sumF should be -sumF
		*sumF = *sumF - 2.0*(i+1.0)*an1*zn/rr;
	}

	return 0;
}

// Calculate the hypergeometric nanotube potentials.
// Modified from Matt LaBrosse's code
int cal_sf_hypergeo(int ii, int iSpecie, int iAtom, double *uij, double *fij)
{
	PSAMPLE_MOLECULE pSampleMole;
	int jj;
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
	double fxij, fyij;

	pSampleMole = sample_mole + iSpecie;
	xxi = xx[ii];
	yyi = yy[ii];
	fxi = fxl[ii];
	fyi = fyl[ii];
	sigmaij = 0.5*(pSampleMole->sigma[iAtom]+(*solid_sigma));
	epsilonij = sqrt(pSampleMole->epsilon[iAtom]*(*solid_epsilon));
	for (jj=0; jj<ntube; jj++)
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
			ljc = const_3pisq_theta*epsilonij*sigmaij*sigmaij;
			kfac = sigmaij*tubeR/(tubeRsq-rijsq);
			zfac = rijsq/tubeRsq;
			hypergeo(afac1, bfac1, cfac, zfac, rij, &hgrep, &hgfrep); // hypergeo series
			kfacsq = kfac*kfac;
			kfac4 = kfacsq*kfacsq;
			kfac10 = kfac4*kfac4*kfacsq;
			rep1 = c21_by_32*kfac10;
			rep2 = -hgfrep + hgrep*(20.0*rij/(tubeRsq-rijsq));
			hypergeo(afac2, bfac2, cfac, zfac, rij, &hgattr, &hgfattr); // calculate hypergeo series
			attr1 = -kfac4;
			attr2 = -hgfattr + hgattr*(8.0*rij/(tubeRsq-rijsq));
			uhyper = ljc*(rep1*hgrep + attr1*hgattr);
			*fij = -ljc*(rep1*rep2 + attr1*attr2);
			fxij = *fij*rxij/rij;
			fyij = *fij*ryij/rij;
		}
		else // outside the tube
		{
			ljc = const_3pisq_theta*epsilonij*sigmaij*sigmaij*tubeR;
			kfac = sigmaij*rij/(rijsq-tubeRsq);
			zfac = tubeRsq/rijsq;
			hypergeo(afac1, bfac1, cfac, zfac, rij, &hgrep, &hgfrep); // calculate the hypergeo series
			kfacsq = kfac*kfac;
			kfac4 = kfacsq*kfacsq;
			kfac10 = kfac4*kfac4*kfacsq;
			rep1 = c21_by_32*kfac10;
			rep2 = hgfrep/rij + hgrep*(-20.0/(rijsq-tubeRsq)+9.0/rijsq);
			hypergeo(afac2, bfac2, cfac, zfac, rij, &hgattr, &hgfattr); // calculate 
			attr1 = -kfac4;
			attr2 = hgfattr/rij + hgattr*(-8.0/(rijsq-tubeRsq)+3.0/rijsq);
			uhyper = ljc*(rep1*hgrep + attr1*hgattr)/rij;
			*fij = -ljc*(rep1*rep2 + attr1*attr2);
			fxij = *fij*rxij/rij;
			fyij = *fij*ryij/rij;
		} // outside tube
		// energy and forces
		*uij = uhyper;
		fxi += fxij;
		fyi += fyij;
	} // loop through all tubes
	fxl[ii] = fxi;
	fyl[ii] = fyi;

	return 0;
}

// calcuate the interactions between solid and fluid
// using atom explicit model
int cal_sf_atom_explicit(int ii, int iSpecie, int iAtom, double *uij,
		double *fij)
{
	int jj;
	double xxi, yyi, zzi;
	double fxi, fyi, fzi;
	double rxij, ryij, rzij;
	double rijsq, rij, r_rijsq;
	double r_r6, r_r12, r_r12_minus_r_r6;
	double epsilonij, sigmaij;
	double usf_vdw_temp, usf_vdw;
	double fxij, fyij, fzij;

	xxi = xx[ii];
	yyi = yy[ii];
	zzi = zz[ii];
	fxi = fxl[ii];
	fyi = fyl[ii];
	fzi = fzl[ii];
	// NOTE: Only uniform solid is considered now.
	// So, calculate epsilonij and sigmaij here.
	// Also assume solid has no charge.
	sigmaij = 0.5*(sample_mole[iSpecie].sigma[iAtom]+(*solid_sigma));
	epsilonij = sqrt(sample_mole[iSpecie].epsilon[iAtom]*(*solid_epsilon));
	usf_vdw = 0.0;
	for (jj=0; jj<solid_natom; jj++)
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
				{
					LJswitch = 1.0;
				}
				else
				{
					LJswitch = (rcutoffsq-rijsq)*(rcutoffsq-rijsq) *(rcutoffsq
							+2.0*rijsq-3.0*rcutonsq) /roff2_minus_ron2_cube;
				}
			}
			r_rijsq = sigmaij*sigmaij/rijsq;
			r_r6 = r_rijsq*r_rijsq*r_rijsq;
			r_r12 = r_r6*r_r6;
			r_r12_minus_r_r6 = r_r12 - r_r6;
			usf_vdw_temp = epsilonij*r_r12_minus_r_r6; // still need *4.0
			if (isLJswitchOn) // if switch is used
			{
				usf_vdw += usf_vdw_temp*LJswitch; // still need 4.0
			}
			else
			{
				usf_vdw += usf_vdw_temp; // still need *4.0
			}
			// force calculation
			if (isLJswitchOn)
			{
				if (rijsq<rcutonsq)
				{
					*fij = 24.0*epsilonij*(r_r12_minus_r_r6+r_r12)/rijsq;
				}
				else
				{
					*fij = 24.0*epsilonij*(r_r12_minus_r_r6+r_r12)/rijsq
							*LJswitch // 4.0 for the real energy
							-4.0*usf_vdw_temp*12.0*(rcutoffsq-rijsq) *(rcutonsq
									-rijsq)/roff2_minus_ron2_cube;
				}
			}
			else
			{
				*fij = 24.0*epsilonij*(r_r12_minus_r_r6+r_r12)/rijsq;
			}
			fxij = *fij*rxij;
			fyij = *fij*ryij;
			fzij = *fij*rzij;
			// force on fluid atom ii
			fxi += fxij;
			fyi += fyij;
			fzi += fzij;
		} // if with cutoff
	} // solid atom jj loop
	fxl[ii] = fxi;
	fyl[ii] = fyi;
	fzl[ii] = fzi;

	// factor
	*uij = usf_vdw*4.0;

	return 0;
}

/**
 * This function accepts reduced parameters. However, they are converted
 * to real units inside the function to call tasos Fortran subroutines.
 * The raw output has units of K and K/Angstrom. They are converted to reduced
 * units and returned.
 */
int call_tasos_forces(int itype, double fxx, double fyy, double fzz, double *usflj_tasos, double *tasos_force)
{ 
	// The tasos energy has unit of K and the force has unit of K/Angstrom. 
	// Change them to reduced units.
	fxx *= sigma_base;
	fyy *= sigma_base;
	fzz *= sigma_base;
	cforce_atom_(&itype, &fxx, &fyy, &fzz, usflj_tasos, tasos_force);
	
	*usflj_tasos /= epsilon_base;
	tasos_force[0] *= (sigma_base/epsilon_base);
	tasos_force[1] *= (sigma_base/epsilon_base);
	tasos_force[2] *= (sigma_base/epsilon_base);
	
	return 0;
}


