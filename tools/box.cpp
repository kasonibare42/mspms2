/***************************************************************************
 *            box.cpp
 *
 *  Sun Jul 17 2005
 *  Created by: Yang Wang
 ****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include "clsdefs.h"
#include "atom.h"
#include "mole.h"
#include "specie.h"
#include "nanotube.h"
#include "box.h"

_simbox::_simbox()
{
    specie_list = NULL;
    nanotube_list = NULL;
    groove_list = NULL;
    uljlrc_term = NULL;
    pljlrc_term = NULL;
    int ii;
    for (ii=0;ii<_MAX_COUNTER_BOX;ii++)
	counter[ii] = 0;
    for (ii=0;ii<_MAX_ACCUMU_BOX;ii++)
	accumu[ii].acc1 
	    = accumu[ii].acc2
	    = accumu[ii].acc3
	    = accumu[ii].acc4
	    = accumu[ii].acc5
	    = 0.0;
    utot = 0.0;
    utotcs = 0.0;
    ukin = 0.0;
    upot = 0.0;
    upotcs = 0.0;
    ulj = 0.0;
    uljcs = 0.0;
    uljlrc = 0.0;
    ucl = 0.0;
    ucllrc = 0.0;
    usflj = 0.0;
    usfljcs = 0.0;
    ubond = 0.0;
    uangle = 0.0;
    udih = 0.0;
    unplj = 0.0;
    unpljcs = 0.0;
    unpcl = 0.0;

    pressure = 0.0;
    pideal = 0.0;
    plj = 0.0;
    pcl = 0.0;
    psflj = 0.0;
    pbond = 0.0;
    pangle = 0.0;
    pdih = 0.0;
    pnplj = 0.0;
    pnpcl = 0.0;
    pljlrc = 0.0;

}

_simbox::~_simbox()
{
    if (specie_list)
    {
	delete []specie_list;
	specie_list = NULL;
    }
    if (nanotube_list)
    {
	delete []nanotube_list;
	nanotube_list = NULL;
    }
    if (groove_list)
    {
	delete []groove_list;
	groove_list = NULL;
    }
    if (uljlrc_term)
    {
	delete []uljlrc_term;
	uljlrc_term = NULL;
    }
    if (pljlrc_term)
    {
	delete []pljlrc_term;
	pljlrc_term = NULL;
    }
}

PSIMMOLE _simbox::get_mole(int index)
{
    int		ii;
    PSIMSPECIE	pSpecie;
    for (ii=0;ii<nspecie;ii++)
    {
	pSpecie = specie_list + ii;
	if (index < pSpecie->nmole)
	{
	    return (pSpecie->mole_list+index);
	}
	index = index - pSpecie->nmole;
    }
}
