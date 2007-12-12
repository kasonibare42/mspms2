/***************************************************************************
 *            specie.cpp
 *
 *  Sun Jul 17 2005
 *  Created by: Yang Wang
 ****************************************************************************/

#include <stdlib.h>
#include <string.h>
#include <iostream>
#include "clsdefs.h"
#include "atom.h"
#include "mole.h"
#include "specie.h"

_simspecie::_simspecie()
{
    strcpy(name, "UNDEF");

    alpha_list			= NULL;
    mole_list 			= NULL;
    // set the address of the pointers for the sample molecule
    // other molecule will copy the sample molecule later
    sample_mole.pBond_list 	= bond_list;
    sample_mole.pAngle_list	= angle_list;
    sample_mole.pDih_list	= dih_list;
    sample_mole.pImp_list	= imp_list;
    sample_mole.pNonpair_list	= nonpair_list;
    sample_mole.pSample_mole	= &sample_mole;

    int ii;
    for (ii=0;ii<_MAX_COUNTER_SPECIE;ii++)
	counter[ii] = 0;
    for (ii=0;ii<_MAX_ACCUMU_SPECIE;ii++)
	accumu[ii].acc1 
	    = accumu[ii].acc2
	    = accumu[ii].acc3
	    = accumu[ii].acc4
	    = accumu[ii].acc5
	    = 0.0;
    /*
    for (ii=0;ii<_MAX_ATOM;ii++)
    {
	sample_mole.atom_list[ii].pNbr_list = nbr_list[ii];
	sample_mole.atom_list[ii].id = ii; // set the atom id properly
    }
    */
}

_simspecie::~_simspecie()
{
    strcpy(name, "UNDEF");
    if (mole_list)
    {
	delete []mole_list;
	mole_list = NULL;
    }
    if (alpha_list)
    {
	delete []alpha_list;
	alpha_list = NULL;
    }

}


void _simspecie::init()
{
    sample_mole.pBond_list 	= bond_list;
    sample_mole.pAngle_list	= angle_list;
    sample_mole.pDih_list	= dih_list;
    sample_mole.pNonpair_list	= nonpair_list;
    sample_mole.pSample_mole	= &sample_mole;
}

