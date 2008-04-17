/**
 * Project: mspms2
 * File: startend.c
 * 
 * Copyright (C) 2008    Yang Wang <ywangd@gmail.com>
 * Created @ 15/04/2008
 * 
 * Description:
 * 
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <ctype.h>
#include <time.h>
#include "mspms2.h"

int opening()
{
	// open output file at the very beginning to keep log of the run
	fpouts = fopen(OUTPUT,"w");

	fprintf(fpouts, "WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW\n");
	fprintf(fpouts, "WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW\n");
	fprintf(fpouts, "WW......WWW.....WW.....WWW......WWW.....WWW\n");
	fprintf(fpouts, "WW..W.W..W..WWWWWW..WW..WW..W.W..W..WWWWWWW\n");
	fprintf(fpouts, "WW..W.W..W..WWWWWW..WW..WW..W.W..W..WWWWWWW\n");
	fprintf(fpouts, "WW..W.W..WW....WWW..WW..WW..W.W..WW....WWWW\n");
	fprintf(fpouts, "WW..W.W..WWWWW..WW..WW..WW..W.W..WWWWW..WWW\n");
	fprintf(fpouts, "WW..W.W..WWWWW..WW..WW..WW..W.W..WWWWW..WWW\n");
	fprintf(fpouts, "WW..WWW..W.....WWW.....WWW..WWW..W.....WWW2\n");
	fprintf(fpouts, "WWWWWWWWWWWWWWWWWW..WWWWWWWWWWWWWWWWWWWWWWW\n");
	fprintf(fpouts, "WWWWWWWWWWWWWWWWWW..WWWWWWWWWWWWWWWWWWWWWWW\n");
	fprintf(fpouts, "WWWWWWWWWWWWWWWWWW..WWWWWWWWWWWWWWWWWWWWWWW\n");
	fprintf(fpouts, "WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWYWANG\n");

	return 0;
}



int ending()
{
	int ii;
	double fExecutionTime;

	fprintf(stderr,"%d frames in the trajectory file.\n",nframe);
	fprintf(fpouts, "%d frames in the trajectory file.\n", nframe);

	fprintf(fpouts,
			"=========================================================\n");
	if (bEquilibrium)
	{
		fprintf(fpouts," EQUILIBRIUM RUN\n");
	}
	else
	{
		fprintf(fpouts, " DATA TAKING RUN\n");
	}
	fprintf(fpouts,
			"=========================================================\n");
	fprintf(fpouts, "Total energy                %15.6le\n", accum[0][5]);
	fprintf(fpouts, "   std. dev.                %15.4lf\n", accum[0][6]);
	fprintf(fpouts, "   fluctuation              %15.4lf\n", accum[0][7]);

	fprintf(fpouts, "Potentail energy            %15.6le\n", accum[1][5]);
	fprintf(fpouts, "   std. dev.                %15.4lf\n", accum[1][6]);
	fprintf(fpouts, "   fluctuation              %15.4lf\n", accum[1][7]);

	fprintf(fpouts, "Kinetic energy              %15.6le\n", accum[2][5]);
	fprintf(fpouts, "   std. dev.                %15.4lf\n", accum[2][6]);
	fprintf(fpouts, "   fluctuation              %15.4lf\n", accum[2][7]);

	fprintf(fpouts, "Inter potential energy      %15.6le\n", accum[3][5]);
	fprintf(fpouts, "   std. dev.                %15.4lf\n", accum[3][6]);
	fprintf(fpouts, "   fluctuation              %15.4lf\n", accum[3][7]);

	fprintf(fpouts, "Intra potential energy      %15.6le\n", accum[4][5]);
	fprintf(fpouts, "   std. dev.                %15.4lf\n", accum[4][6]);
	fprintf(fpouts, "   fluctuation              %15.4lf\n", accum[4][7]);

	fprintf(fpouts, "LJ energy                   %15.6le\n", accum[5][5]);
	fprintf(fpouts, "   std. dev.                %15.4lf\n", accum[5][6]);
	fprintf(fpouts, "   fluctuation              %15.4lf\n", accum[5][7]);

	fprintf(fpouts, "Bond energy                 %15.6le\n", accum[6][5]);
	fprintf(fpouts, "   std. dev.                %15.4lf\n", accum[6][6]);
	fprintf(fpouts, "   fluctuation              %15.4lf\n", accum[6][7]);

	fprintf(fpouts, "Angle energy                %15.6le\n", accum[7][5]);
	fprintf(fpouts, "   std. dev.                %15.4lf\n", accum[7][6]);
	fprintf(fpouts, "   fluctuation              %15.4lf\n", accum[7][7]);

	fprintf(fpouts, "Dihedral energy             %15.6le\n", accum[8][5]);
	fprintf(fpouts, "   std. dev.                %15.4lf\n", accum[8][6]);
	fprintf(fpouts, "   fluctuation              %15.4lf\n", accum[8][7]);

	fprintf(fpouts, "Improper energy             %15.6le\n", accum[9][5]);
	fprintf(fpouts, "   std. dev.                %15.4lf\n", accum[9][6]);
	fprintf(fpouts, "   fluctuation              %15.4lf\n", accum[9][7]);

	fprintf(fpouts, "Ewald energy                %15.6le\n", accum[10][5]);
	fprintf(fpouts, "   std. dev.                %15.4lf\n", accum[10][6]);
	fprintf(fpouts, "   fluctuation              %15.4lf\n", accum[10][7]);

	fprintf(fpouts, "Real part energy            %15.6le\n", accum[11][5]);
	fprintf(fpouts, "   std. dev.                %15.4lf\n", accum[11][6]);
	fprintf(fpouts, "   fluctuation              %15.4lf\n", accum[11][7]);

	fprintf(fpouts, "Fourier part energy         %15.6le\n", accum[12][5]);
	fprintf(fpouts, "   std. dev.                %15.4lf\n", accum[12][6]);
	fprintf(fpouts, "   fluctuation              %15.4lf\n", accum[12][7]);

	fprintf(fpouts, "Self part energy            %15.6le\n", accum[13][5]);
	fprintf(fpouts, "   std. dev.                %15.4lf\n", accum[13][6]);
	fprintf(fpouts, "   fluctuation              %15.4lf\n", accum[13][7]);

	fprintf(fpouts, "Vaccum energy               %15.6le\n", accum[16][5]);
	fprintf(fpouts, "   std. dev.                %15.4lf\n", accum[16][6]);
	fprintf(fpouts, "   fluctuation              %15.4lf\n", accum[16][7]);

	fprintf(fpouts, "Solid fluid energy          %15.6le\n", accum[14][5]);
	fprintf(fpouts, "   std. dev.                %15.4lf\n", accum[14][6]);
	fprintf(fpouts, "   fluctuation              %15.4lf\n", accum[14][7]);

	fprintf(fpouts, "Temperature                 %15.6le\n", accum[15][5]);
	fprintf(fpouts, "   std. dev.                %15.4lf\n", accum[15][6]);
	fprintf(fpouts, "   fluctuation              %15.4lf\n", accum[15][7]);

	fprintf(fpouts, "Wolf energy                 %15.6le\n", accum[17][5]);
	fprintf(fpouts, "   std. dev.                %15.4lf\n", accum[17][6]);
	fprintf(fpouts, "   fluctuation              %15.4lf\n", accum[17][7]);

	fprintf(fpouts, "Pressure                    %15.6le\n", accum[18][5]);
	fprintf(fpouts, "   std. dev.                %15.4lf\n", accum[18][6]);
	fprintf(fpouts, "   fluctuation              %15.4lf\n", accum[18][7]);

	fprintf(fpouts, "Box volume                  %15.6le\n", accum[19][5]);
	fprintf(fpouts, "   std. dev.                %15.4lf\n", accum[19][6]);
	fprintf(fpouts, "   fluctuation              %15.4lf\n", accum[19][7]);

	fprintf(fpouts, "Ideal pressure              %15.6le\n", accum[20][5]);
	fprintf(fpouts, "   std. dev.                %15.4lf\n", accum[20][6]);
	fprintf(fpouts, "   fluctuation              %15.4lf\n", accum[20][7]);

	if (what_simulation == HYBRID_MONTE_CARLO)
	{
		fprintf(fpouts, "Total canonical moves       %15d\n", counts[20]);
		fprintf(fpouts, "accepted canonical moves    %15d\n", counts[21]);
		fprintf(fpouts, "   ratio                    %15.4lf\n", counts[21]
				*1.0/counts[20]);
		fprintf(fpouts, "   delt                     %15.4lf\n", delt);

		if (pvolm > 0.0)
		{
			fprintf(fpouts, "Total volume changes        %15d\n", counts[23]);
			fprintf(fpouts, "accepted volume changes     %15d\n", counts[24]);
			fprintf(fpouts, "   ratio                    %15.4lf\n",
					counts[24] *1.0/counts[23]);
			fprintf(fpouts, "   delv                     %15.4lf\n", delv);
		}
	}

	fprintf(fpouts,
			"=========================================================\n");

    // If this is still in equilibrium, do not free memories
	if (bEquilibrium)
	{
		return 0;
	}
	
	// release the dynamically allocated memory for saving old positions for HMC simulation
	if (what_simulation == HYBRID_MONTE_CARLO || what_simulation==SIMULATED_ANNEALING)
	{
		free(xx_old);
		free(yy_old);
		free(zz_old);
	}

	if (iSF_type==SF_NANOTUBE_HYPERGEO)
	{
		free(hgntc_xx);
		free(hgntc_yy);
		free(hgnt_radius);
	}
	else if (iSF_type==SF_NANOTUBE_ATOM_EXPLICIT)
	{
		free(solid_sigma);
		free(solid_epsilon);
		free(solid_charge);
		free(solid_xx);
		free(solid_yy);
		free(solid_zz);
	}
	else if (iSF_type==SF_NANOTUBE_MY_INTERP)
	{
		// free memories
		free(interp_vector);
		for (ii=0; ii<NUNIQUE_ATOM_MAX; ii++)
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

	// Calculate the execution duration of the program
	fExecutionTime = clock()*1.0/CLOCKS_PER_SEC;
	fprintf(stderr, "Total execution time: %20.2lf seconds.\n", fExecutionTime);
	fprintf(fpouts, "Total execution time: %20.2lf seconds.\n", fExecutionTime);

	// close files
	fclose(fplog);
	fclose(fptrj);
	fclose(fpouts);

	return 0;
}
