/**
 * Project: mspms2
 * File: misc.c
 * 
 * Copyright (C) 2008    Yang Wang <ywangd@gmail.com>
 * Created @ 15/04/2008
 * 
 * Description:
 * 		Anything that does not have a perfect place to go goes here.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <ctype.h>
#include "mspms2.h"

int calculate_ljlrc()
{
	int mm, nn;

	// calculate the lrc for the real system
	uljlrc = 0.0;
	pljlrc = 0.0;
	for (mm=0; mm<nspecie; mm++)
	{
		for (nn=0; nn<nspecie; nn++)
		{
			uljlrc += (uljlrc_term[mm][nn]*nmole_per_specie[mm]
					*nmole_per_specie[nn]/boxv);
			pljlrc += (pljlrc_term[mm][nn]*nmole_per_specie[mm]
					*nmole_per_specie[nn]/boxv/boxv);
		}
	}
	pljlrc = pljlrc*J_PER_MOL_A3_TO_PA; // convert J/mol/A^3 to Pascal

	return 0;
}


