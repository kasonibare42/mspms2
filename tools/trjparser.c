#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>

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
    char *outputfile; // output msd file
    int nframe_total, nframe_span, nframe_interval;

    int natom, natomper;
    int nmole; // nmole = natom/natomper
    double xx[natom_max], yy[natom_max], zz[natom_max];
    double xx0[natom_max], yy0[natom_max], zz0[natom_max];
    double mole_xx[nmole_max], mole_yy[nmole_max], mole_zz[nmole_max];
    double mole_xx0[nmole_max], mole_yy0[nmole_max], mole_zz0[nmole_max];
    double atom_weight[natomper_max];
    double mole_weight;

    int last_origin; // the frame number of last origin
    int norigin; // number of origins
    int origin_frame;
    int starting_frame, ending_frame;
    double Dsx[nframe_span_max], Dsy[nframe_span_max], Dsz[nframe_span_max];
    double D0x[nframe_span_max], D0y[nframe_span_max], D0z[nframe_span_max]; 
    double D0x_temp, D0y_temp, D0z_temp;
    int atom_idx;


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
	    sscanf(argv[++ii],"%d,%d,%d",&nframe_total,&nframe_span,&nframe_interval);
	else if (!strcmp(argv[ii],"-natom"))
	    natom = atoi(argv[++ii]);
	else if (!strcmp(argv[ii],"-natomper"))
	    natomper = atoi(argv[++ii]);
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

    // calculate number of molecules
    nmole = natom/natomper;
    printf("number of molecules = %d\n",nmole);

    assert(natom<natom_max);
    assert(nmole<nmole_max);
    assert(nframe_span<nframe_span_max);

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

    // calculate number of origins
    norigin = 0;
    origin_frame = 0;
    while (origin_frame+nframe_span<=nframe_total)
    {
	origin_frame += nframe_interval;
	norigin++;
    }
    printf("Number of origins = %d (0-%d)\n",norigin,norigin-1);

    // calculate msd
    for (ii=0;ii<nframe_span;ii++)
    {
	Dsx[ii] = 0.0;
	Dsy[ii] = 0.0;
	Dsz[ii] = 0.0;
	D0x[ii] = 0.0;
	D0y[ii] = 0.0;
	D0z[ii] = 0.0;
    }
    for (ii=0;ii<norigin;ii++)
    {
	starting_frame = ii*nframe_interval;
	ending_frame = starting_frame + nframe_span - 1;

	printf("Now calculating origin No.%d. Starting from frame No.%d to ending frame No.%d\n",ii,starting_frame,ending_frame);

	// ignore the frames before starting_frame 
	fptrj = fopen(trjfile,"rb");
	for (jj=0;jj<starting_frame;jj++)
	{
	    fread(xx,sizeof(double),natom,fptrj);
	    fread(yy,sizeof(double),natom,fptrj);
	    fread(zz,sizeof(double),natom,fptrj);
	}
	// read the starting frame
	fread(xx0,sizeof(double),natom,fptrj);
	fread(yy0,sizeof(double),natom,fptrj);
	fread(zz0,sizeof(double),natom,fptrj);
	// calculate the molecular coordinates
	for (jj=0;jj<nmole;jj++)
	{
	    mole_xx0[jj] = 0.0;
	    mole_yy0[jj] = 0.0;
	    mole_zz0[jj] = 0.0;
	    for (kk=0;kk<natomper;kk++)
	    {
		atom_idx = jj*natomper + kk;
		mole_xx0[jj] += xx0[atom_idx]*atom_weight[kk];
		mole_yy0[jj] += yy0[atom_idx]*atom_weight[kk];
		mole_zz0[jj] += zz0[atom_idx]*atom_weight[kk];
	    }
	    mole_xx0[jj] /= mole_weight;
	    mole_yy0[jj] /= mole_weight;
	    mole_zz0[jj] /= mole_weight;
	}
	// read other frames inside the span
	for (ff=1;ff<nframe_span;ff++) // starting from 1 since 0 is already read
	{
	    fread(xx,sizeof(double),natom,fptrj);
	    fread(yy,sizeof(double),natom,fptrj);
	    fread(zz,sizeof(double),natom,fptrj);
	    // calculate the molecular coordinates
	    for (jj=0;jj<nmole;jj++)
	    {
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

	    D0x_temp = 0.0;
	    D0y_temp = 0.0;
	    D0z_temp = 0.0;
	    for (jj=0;jj<nmole;jj++)
	    {
		// at time (frame) ff, calculate the Ds
		Dsx[ff] += (mole_xx[jj]-mole_xx0[jj])*(mole_xx[jj]-mole_xx0[jj]);
		Dsy[ff] += (mole_yy[jj]-mole_yy0[jj])*(mole_yy[jj]-mole_yy0[jj]);
		Dsz[ff] += (mole_zz[jj]-mole_zz0[jj])*(mole_zz[jj]-mole_zz0[jj]);

		// at time (frame) ff, calculate D0
		D0x_temp += (mole_xx[jj]-mole_xx0[jj]);
		D0y_temp += (mole_yy[jj]-mole_yy0[jj]);
		D0z_temp += (mole_zz[jj]-mole_zz0[jj]);
	    }
	    D0x[ff] += D0x_temp*D0x_temp;
	    D0y[ff] += D0y_temp*D0y_temp;
	    D0z[ff] += D0z_temp*D0z_temp;
	} // loop through nframe_span inside one span
	fclose(fptrj);

    } // loop through all origins

    // average over multiple origins and number of molecules
    fpouts = fopen(outputfile,"w");
    for (ff=0;ff<nframe_span;ff++) 
    {
	Dsx[ff] = Dsx[ff]/norigin/nmole;
	Dsy[ff] = Dsy[ff]/norigin/nmole;
	Dsz[ff] = Dsz[ff]/norigin/nmole;

	D0x[ff] = D0x[ff]/norigin/nmole;
	D0y[ff] = D0y[ff]/norigin/nmole;
	D0z[ff] = D0z[ff]/norigin/nmole;

	fprintf(fpouts,"%d   %lf   %lf   %lf   %lf   %lf   %lf\n",ff,Dsx[ff],Dsy[ff],Dsz[ff],D0x[ff],D0y[ff],D0z[ff]);

    }
    fclose(fpouts);




}





