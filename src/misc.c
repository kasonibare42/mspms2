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

void get_specie_and_relative_atom_id(int abs_atom_id, int *specie_id, int *sample_atom_id)
{
	int ii, iAtom, iSampleAtom;
	
	iAtom = 0;
	for (ii=0; ii<nspecie; ii++)
	{
		if (abs_atom_id < iAtom+natom_per_specie[ii])
		{
			// We know this atom is for this specie
			*specie_id = ii;
			*sample_atom_id = (abs_atom_id-iAtom)%sample_mole[ii].natom;
			break;
		}
		iAtom += natom_per_specie[ii];
	}
}

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
			uljlrc += (uljlrc_term[mm][nn]*nmole_per_specie[mm]*nmole_per_specie[nn]/boxv);
			pljlrc += (pljlrc_term[mm][nn]*nmole_per_specie[mm]*nmole_per_specie[nn]/boxv/boxv);
		}
	}

	return 0;
}

