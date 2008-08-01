#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#define NATOM_MAX 9999
#define NMOLE_MAX 1000
#define NSPECIE_MAX 3
#define SAMPLE_NATOM_MAX	200
#define nframe_span_max 10000 // number of frames in one span

int usage()
{}

int main(int argc, char *argv[])
{
	int ii, jj, kk, ff, mm;
	char *trjfile; // trajectory file
	char *outputfile; // output msd file
	int nframe_total, nframe_span, nframe_interval;

	int nspecie;
	int natom;
	int nmole; 
	int iAtom, iMole, awidstart;
	int nmole_per_specie[NSPECIE_MAX];
	int sample_natom_per_mole[NSPECIE_MAX];
	char sample_atomname[SAMPLE_NATOM_MAX][5];
	double sample_aw[SAMPLE_NATOM_MAX];
	double sample_mw[NSPECIE_MAX];
	double xx[NATOM_MAX], yy[NATOM_MAX], zz[NATOM_MAX];
	double xx0[NATOM_MAX], yy0[NATOM_MAX], zz0[NATOM_MAX];
	double mole_xx[NMOLE_MAX], mole_yy[NMOLE_MAX], mole_zz[NMOLE_MAX];
	double mole_xx0[NMOLE_MAX], mole_yy0[NMOLE_MAX], mole_zz0[NMOLE_MAX];

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

	assert(nframe_span<nframe_span_max);

	// calculate number of origins
	norigin = 0;
	origin_frame = 0;
	while (origin_frame+nframe_span<=nframe_total)
	{
		origin_frame += nframe_interval;
		norigin++;
	}
	printf("Number of origins = %d (0-%d)\n",norigin,norigin-1);

	// Zeros
	for (ii=0;ii<nframe_span;ii++)
	{
		Dsx[ii] = 0.0;
		Dsy[ii] = 0.0;
		Dsz[ii] = 0.0;
		D0x[ii] = 0.0;
		D0y[ii] = 0.0;
		D0z[ii] = 0.0;
	}

	// Calculate MSD
	for (ii=0;ii<norigin;ii++)
	{
		starting_frame = ii*nframe_interval;
		ending_frame = starting_frame + nframe_span - 1;

		// ignore the frames before starting_frame 
		fptrj = fopen(trjfile,"rb");
		fread(&nspecie, sizeof(int), 1, fptrj);
		fread(&nmole, sizeof(int), 1, fptrj);
		fread(&natom, sizeof(int), 1, fptrj);
		assert(natom<NATOM_MAX);
		assert(nmole<NMOLE_MAX);
		assert(nspecie<NSPECIE_MAX);
		iAtom = 0;
		for (kk=0;kk<nspecie;kk++)
		{
			fread(&nmole_per_specie[kk], sizeof(int), 1, fptrj);
			fread(&sample_natom_per_mole[kk], sizeof(int), 1, fptrj);
			printf("sample_natom_per_mole[%d]=%d\n",kk,sample_natom_per_mole[kk]);
			sample_mw[kk] = 0.0;
			for (jj=iAtom;jj<iAtom+sample_natom_per_mole[kk];jj++)
			{
				fread(sample_atomname[jj], sizeof(char), 5, fptrj);
				fread(&sample_aw[jj], sizeof(double), 1, fptrj);
				sample_mw[kk] += sample_aw[jj];
			}
			iAtom += sample_natom_per_mole[kk];
		}

		if (ii==0)
		{
			printf("nspecie=%d, nmole=%d, natom=%d\n", nspecie, nmole, natom);
			iAtom = 0;
			for (kk=0;kk<nspecie;kk++)
			{
				for (mm=0;mm<sample_natom_per_mole[kk];mm++)
				{
					printf("specie %d: Atom %d: AW=%lf\n",kk,mm,sample_aw[iAtom+mm]);
				}
				printf("specie %d: MW=%lf\n",kk,sample_mw[kk]);
				iAtom += sample_natom_per_mole[mm];
			}
		}

		printf("Now calculating origin No.%d. Starting from frame No.%d to ending frame No.%d\n",ii,starting_frame,ending_frame);
		for (jj=0;jj<starting_frame;jj++)
		{
			fread(nmole_per_specie, sizeof(int), nspecie, fptrj); 
			natom = 0; 
			nmole = 0;
			for (kk=0;kk<nspecie;kk++)
			{
				natom += nmole_per_specie[kk]*sample_natom_per_mole[kk]; 
				nmole += nmole_per_specie[kk];
			} 
			assert(natom<NATOM_MAX);
			assert(nmole<NMOLE_MAX);
			fread(xx,sizeof(double),natom,fptrj);
			fread(yy,sizeof(double),natom,fptrj);
			fread(zz,sizeof(double),natom,fptrj);
		}

		fread(nmole_per_specie, sizeof(int), nspecie, fptrj); 
		natom = 0;
		nmole = 0;
		for (kk=0;kk<nspecie;kk++)
		{
			natom += nmole_per_specie[kk]*sample_natom_per_mole[kk];
			nmole += nmole_per_specie[kk];
		}
		assert(natom<NATOM_MAX);
		assert(nmole<NMOLE_MAX);
		// read the starting frame
		fread(xx0,sizeof(double),natom,fptrj);
		fread(yy0,sizeof(double),natom,fptrj);
		fread(zz0,sizeof(double),natom,fptrj);

		// calculate the molecular coordinates
		awidstart = 0;
		iMole = 0;
		for (mm=0;mm<nspecie;mm++)
		{
			for (jj=iMole;jj<iMole+nmole_per_specie[mm];jj++)
			{
				mole_xx0[jj] = 0.0;
				mole_yy0[jj] = 0.0;
				mole_zz0[jj] = 0.0;
				for (kk=0;kk<sample_natom_per_mole[mm];kk++)
				{
					atom_idx = jj*sample_natom_per_mole[mm] + kk;
					mole_xx0[jj] += xx0[atom_idx]*sample_aw[awidstart+kk];
					mole_yy0[jj] += yy0[atom_idx]*sample_aw[awidstart+kk];
					mole_zz0[jj] += zz0[atom_idx]*sample_aw[awidstart+kk];
				}
				mole_xx0[jj] /= sample_mw[mm];
				mole_yy0[jj] /= sample_mw[mm];
				mole_zz0[jj] /= sample_mw[mm];
			}
			awidstart += sample_natom_per_mole[mm];
			iMole += nmole_per_specie[mm];
		}

		// read other frames inside the span
		for (ff=1;ff<nframe_span;ff++) // starting from 1 since 0 is already read
		{
			fread(nmole_per_specie, sizeof(int), nspecie, fptrj); 
			natom = 0; 
			nmole = 0;
			for (kk=0;kk<nspecie;kk++)
			{
				natom += nmole_per_specie[kk]*sample_natom_per_mole[kk]; 
				nmole += nmole_per_specie[kk];
			} 
			assert(natom<NATOM_MAX);
			assert(nmole<NMOLE_MAX);
			fread(xx,sizeof(double),natom,fptrj);
			fread(yy,sizeof(double),natom,fptrj);
			fread(zz,sizeof(double),natom,fptrj);
			// calculate the molecular coordinates
			awidstart = 0;
			iMole = 0;
			for (mm=0;mm<nspecie;mm++)
			{
				for (jj=iMole;jj<iMole+nmole_per_specie[mm];jj++)
				{
					mole_xx[jj] = 0.0;
					mole_yy[jj] = 0.0;
					mole_zz[jj] = 0.0;
					for (kk=0;kk<sample_natom_per_mole[mm];kk++)
					{
						atom_idx = jj*sample_natom_per_mole[mm] + kk;
						mole_xx[jj] += xx[atom_idx]*sample_mw[awidstart+kk];
						mole_yy[jj] += yy[atom_idx]*sample_mw[awidstart+kk];
						mole_zz[jj] += zz[atom_idx]*sample_mw[awidstart+kk];
					}
					mole_xx[jj] /= sample_mw[mm];
					mole_yy[jj] /= sample_mw[mm];
					mole_zz[jj] /= sample_mw[mm];
				}
				awidstart += sample_natom_per_mole[mm];
				iMole += nmole_per_specie[mm];
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
				// printf("D0x=%lf\n",D0x[ff]);
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





