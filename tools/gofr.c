// This program calculates the g(r)
// GNU GSL histograms are used
// Compile with
// gcc gofr.c -o gofr.x -lgsl -lgslcblas -lm

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <gsl/gsl_histogram.h>

#define PI	3.14159265358979
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
    char OutputFileFullName[128];
    int nframe_start, nframe_end;

    int natom, natomper;
    int nmole; // nmole = natom/natomper
    double xx[natom_max], yy[natom_max], zz[natom_max];
    double mole_xx[nmole_max], mole_yy[nmole_max], mole_zz[nmole_max];
    double atom_weight[natomper_max];
    double mole_weight;
    int atom_idx;

    double xmin, xmax;
    double xdiff, ydiff, zdiff, rdiff;
    int nbin;
    int maxval_bin_idx;
    int sum_of_bins;

    double constant, rho=1.0, gofr;

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
	    sscanf(argv[++ii],"%d,%d",&nframe_start,&nframe_end);
	else if (!strcmp(argv[ii],"-natom"))
	    natom = atoi(argv[++ii]);
	else if (!strcmp(argv[ii],"-natomper"))
	    natomper = atoi(argv[++ii]);
	else if (!strcmp(argv[ii],"-xmax"))
	    xmax = atof(argv[++ii]);
	else if (!strcmp(argv[ii],"-nbin"))
	    nbin = atoi(argv[++ii]);
	else if (!strcmp(argv[ii],"-rho"))
	    rho = atof(argv[++ii]);
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

    xmin = 0.0;
    printf("set xmin = %lf, xmax = %lf    with %d bins.\n",xmin,xmax,nbin);
    gsl_histogram *h;

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
    printf("ignoring first %d frames...\n",nframe_start);
    for (ii=0;ii<nframe_start;ii++)
    {
	fread(xx,sizeof(double),natom,fptrj);
	fread(yy,sizeof(double),natom,fptrj);
	fread(zz,sizeof(double),natom,fptrj);
    }

    // set xmin and xmax
    h = gsl_histogram_alloc(nbin);
    gsl_histogram_set_ranges_uniform(h,xmin,xmax);

    // calculate the molecular to wall distance
    printf("starting gofr calculation...\n");
    for (ii=nframe_start;ii<nframe_end;ii++)
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
	}
	for (jj=0;jj<nmole-1;jj++)
	{
	    for (kk=jj+1;kk<nmole;kk++)
	    {

		// calculate to wall distance
		xdiff = mole_xx[jj] - mole_xx[kk];
		ydiff = mole_yy[jj] - mole_yy[kk];
		zdiff = mole_zz[jj] - mole_zz[kk];
		rdiff = sqrt(xdiff*xdiff + ydiff*ydiff + zdiff*zdiff);

		// update histogram
		gsl_histogram_increment(h,rdiff);
	    }
	}

    }

    maxval_bin_idx = gsl_histogram_max_bin(h);
    sum_of_bins = gsl_histogram_sum(h);
    printf("max value of %lf in bin no.%d of range from %lf to %lf\n",
	    gsl_histogram_max_val(h),maxval_bin_idx,h->range[maxval_bin_idx],h->range[maxval_bin_idx+1]);
    printf("sum of all bin value is %lf\n",gsl_histogram_sum(h));

    sprintf(OutputFileFullName,"%s-%d-%d",outputfile,nframe_start,nframe_end);
    fpouts = fopen(OutputFileFullName,"w");
    // gsl_histogram_fprintf(fpouts, h, "%g", "%g");

    constant = PI*rho*natom*(nframe_end-nframe_start)*4.0/3.0;
    for (ii=0;ii<h->n;ii++)
    {	
	gofr = h->bin[ii]/(pow(h->range[ii+1],3.0)-pow(h->range[ii],3.0))/constant;
	fprintf(fpouts, "%d %lf %lf %lf %lf\n",ii,h->range[ii],h->range[ii+1],gofr,h->bin[ii]);
    }

    gsl_histogram_free(h);
    fclose(fpouts);


    fclose(fptrj);
}

