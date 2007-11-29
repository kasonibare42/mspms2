// This program calculates the histogram of the distannce of 
// the molecular center to the tube wall
// GNU GSL histograms are used
// Compile with
// gcc histogramtube.c -o histogramtube.x -lgsl -lgslcblas -lm

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <gsl/gsl_histogram.h>

#define natom_max 9999
#define nmole_max 1000
#define natomper_max 100 // number of atoms per molecule
#define nframe_span_max 50000 // number of frames in one span

int usage()
{}

int main(int argc, char *argv[])
{
    int ii, jj, kk, ff;
    char *trjfile; // trajectory file
    char *outputfile; // output histogram file
    int nframe_total, nframe_ignore;

    int natom, natomper;
    int nmole; // nmole = natom/natomper
    double xx[natom_max], yy[natom_max], zz[natom_max];
    double mole_xx[nmole_max], mole_yy[nmole_max], mole_zz[nmole_max];
    double atom_weight[natomper_max];
    double mole_weight;
    int atom_idx;

    double x_center, y_center, r_tube;
    double xmin, xmax;
    double xdiff, ydiff;
    int nbin;
    double mole_to_wall_dist;
    int maxval_bin_idx;
    int sum_of_bins;

    FILE *fptrj, *fpouts;

    for (ii=1;ii<argc;ii++)
    {
	if (!strcmp(argv[ii],"-in"))
	{
	    trjfile = argv[++ii];
	    printf("%s\n",trjfile);
	}
	else if (!strcmp(argv[ii],"-out"))
	    outputfile = argv[++ii];
	else if (!strcmp(argv[ii],"-for"))
	    sscanf(argv[++ii],"%d,%d",&nframe_total,&nframe_ignore);
	else if (!strcmp(argv[ii],"-natom"))
	    natom = atoi(argv[++ii]);
	else if (!strcmp(argv[ii],"-natomper"))
	    natomper = atoi(argv[++ii]);
	else if (!strcmp(argv[ii],"-xyr"))
	    sscanf(argv[++ii],"%lf,%lf,%lf",&x_center,&y_center,&r_tube);
	else if (!strcmp(argv[ii],"-nbin"))
	    nbin = atoi(argv[++ii]);
	else if (!strcmp(argv[ii],"-h"))
	{
	    usage();
	    exit(0);
	}
	else
	{
	    fprintf(stderr,"Error: Command-line argument '%s' not recognized.\n", argv[ii]);
	    exit(1);
	}
    }

    // set xmin and xmax
    xmin = 0.0;
    xmax = sqrt(x_center*x_center+y_center*y_center);
    gsl_histogram *h = gsl_histogram_alloc(nbin);
    gsl_histogram_set_ranges_uniform(h,xmin,xmax);
    printf("set xmin = %lf, xmax = %lf    with %d bins.\n",xmin,xmax,nbin);

    // calculate number of molecules
    nmole = natom/natomper;
    printf("number of molecules = %d\n",nmole);

    assert(natom<natom_max);
    assert(nmole<nmole_max);

    // read in atom weight info
    mole_weight = 0.0;
    for (ii=0;ii<natomper;ii++)
    {
	printf("weight of atom %d = ",ii+1);
	scanf("%lf",&atom_weight[ii]);
	printf("recieved input: %lf\n",atom_weight[ii]);
	mole_weight += atom_weight[ii];
    }
    printf("molecule weight = %lf\n",mole_weight);

    fptrj = fopen(trjfile,"rb");
    // ignore certain number of frames
    printf("ignoring first %d frames...\n",nframe_ignore);
    for (ii=0;ii<nframe_ignore;ii++)
    {
	fread(xx,sizeof(double),natom,fptrj);
	fread(yy,sizeof(double),natom,fptrj);
	fread(zz,sizeof(double),natom,fptrj);
    }

    // calculate the molecular to wall distance
    printf("starting histogram calculation...\n");
    for (ii=nframe_ignore;ii<nframe_total;ii++)
    {
	fread(xx,sizeof(double),natom,fptrj);
	fread(yy,sizeof(double),natom,fptrj);
	fread(zz,sizeof(double),natom,fptrj);

	for (jj=0;jj<nmole;jj++)
	{
	    // calculate the molecular coordinates
	    mole_xx[jj] = 0.0;
	    mole_yy[jj] = 0.0;
	    mole_zz[jj] = 0.0;
	    for (kk=0;kk<natomper;kk++)
	    {
		atom_idx = jj*natomper + kk;
		mole_xx[jj] += xx[atom_idx]*atom_weight[kk];
		mole_yy[jj] += yy[atom_idx]*atom_weight[kk];
		mole_zz[jj] += zz[atom_idx]*atom_weight[kk];
	    }
	    mole_xx[jj] /= mole_weight;
	    mole_yy[jj] /= mole_weight;
	    mole_zz[jj] /= mole_weight;

	    // calculate to wall distance
	    xdiff = mole_xx[jj] - x_center;
	    ydiff = mole_yy[jj] - y_center;
	    mole_to_wall_dist = r_tube - sqrt(xdiff*xdiff + ydiff*ydiff);

	    // update histogram
	    gsl_histogram_increment(h,mole_to_wall_dist);
	}

    }

    maxval_bin_idx = gsl_histogram_max_bin(h);
    sum_of_bins = gsl_histogram_sum(h);
    printf("max value of %lf in bin no.%d of range from %lf to %lf\n",
	    gsl_histogram_max_val(h),maxval_bin_idx,h->range[maxval_bin_idx],h->range[maxval_bin_idx+1]);
    printf("sum of all bin value is %lf\n",gsl_histogram_sum(h));

    fpouts = fopen(outputfile,"w");
    // gsl_histogram_fprintf(fpouts, h, "%g", "%g");

    for (ii=0;ii<h->n;ii++)
    {	
	fprintf(fpouts, "%d %lf %lf %lf %lf\n",ii,h->range[ii],h->range[ii+1],h->bin[ii],h->bin[ii]/sum_of_bins);
    }

    gsl_histogram_free(h);
    fclose(fptrj);
    fclose(fpouts);
}
