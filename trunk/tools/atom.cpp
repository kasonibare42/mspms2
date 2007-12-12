/***************************************************************************
 *            atom.cpp
 *
 *  Sun Jul 17 2005
 *  Created by: Yang Wang
 ****************************************************************************/
#include <stdlib.h>
#include <string.h>
#include "clsdefs.h"
#include "atom.h"

_simatom::_simatom()
{
    strcpy(name, "UNDEF");
    strcpy(type, "UNDEF");

    pNbr_list 		= NULL;
    pAlpha_list 	= NULL;

    xdisp = 0.0;
    ydisp = 0.0;
    zdisp = 0.0;

    fx = fy = fz = 0.0;
    fxra = fyra = fzra = 0.0;
    fxer = fyer = fzer = 0.0;
    fx_lj = fy_lj = fz_lj = 0.0;
    fx_nplj = fy_nplj = fz_nplj = 0.0;
    fx_cl = fy_cl = fz_cl = 0.0;
    fx_npcl = fy_npcl = fz_npcl = 0.0;
    fx_bond = fy_bond = fz_bond = 0.0;
    fx_angle = fy_angle = fz_angle = 0.0;
    fx_dih = fy_dih = fz_dih = 0.0;
    fx_sflj = fy_sflj = fz_sflj = 0.0;

    fPosition = 0; // UNDEFINED
}

_simatom::~_simatom()
{}

void _simatom::check_atom_weight()
{}

void _simatom::check_atom_props()
{}

void _simatom::check_atom_alpha()
{}

