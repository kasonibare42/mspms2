/***************************************************************************
 *            mole.cpp
 *
 *  Sun Jul 17 2005
 *  Created by: Yang Wang
 ****************************************************************************/

#include <stdlib.h>
#include <iostream>
#include <string.h>
#include <math.h>
#include "clsdefs.h"
#include "atom.h"
#include "mole.h"

using namespace std;

_simmole::_simmole()
{
    pXx			= NULL;
    pYy			= NULL;
    pZz			= NULL;
    pBond_list		= NULL;
    pAngle_list		= NULL;
    pDih_list		= NULL;
    pNonpair_list	= NULL;
    pSample_mole	= NULL;

    atom_list		= NULL;

    fPosition		= 0; // UNDEFINED
}

_simmole::~_simmole()
{
    pXx			= NULL;
    pYy			= NULL;
    pZz			= NULL;
    pBond_list		= NULL;
    pAngle_list		= NULL;
    pDih_list		= NULL;
    pNonpair_list	= NULL;
    pSample_mole	= NULL;

    if (atom_list)
    {
	delete []atom_list;
	atom_list = NULL;
    }
}

void _simmole::cal_com()
{
    int		ii;
    PSIMATOM	pAtom;
    double	wgt;

    xx = yy = zz = 0.0;
    for (ii=0;ii<natom;ii++)
    {
	pAtom = atom_list + ii;
	wgt = pAtom->weight;
	xx += pAtom->xx*wgt;
	yy += pAtom->yy*wgt;
	zz += pAtom->zz*wgt;
    }
    xx = xx/weight;
    yy = yy/weight;
    zz = zz/weight;
}

void _simmole::cal_all_atom_xyz2efg()
{
    // always use (1,0,0), (0,1,0), (0,0,1) as the orientation vectors.
    // so the calculations are reduced as follows:
    int 	ii;
    PSIMATOM	pAtom;
    for (ii=0;ii<natom;ii++)
    {
	pAtom = atom_list + ii;
	if (ii == iCenter_atom)
	    pAtom->ee = pAtom->ff = pAtom->gg = 0.0;
	else
	{
	    pAtom->ee = pAtom->xx - *pXx;
	    pAtom->ff = pAtom->yy - *pYy;
	    pAtom->gg = pAtom->zz - *pZz;
	}
    }
}

void _simmole::cal_all_atom_efg2xyz()
{
    // always use (1,0,0), (0,1,0), (0,0,1) as the orientation vectors.
    // so the calculations are reduced as follows
    int 	ii;
    PSIMATOM	pAtom;
    for (ii=0;ii<natom;ii++)
    {
	pAtom = atom_list + ii;
	if (ii != iCenter_atom)
	{
	    pAtom->xx = pAtom->ee + *pXx;
	    pAtom->yy = pAtom->ff + *pYy;
	    pAtom->zz = pAtom->gg + *pZz;
	}
    }
}
void _simmole::cal_all_atom_efg2xyz_ortn(double ex, double ey, double ez, double fx, double fy, double fz, double gx, double gy, double gz)
{
    int		ii;
    PSIMATOM	pAtom;
    for (ii=0;ii<natom;ii++)
    {
	pAtom = atom_list + ii;
	if (ii != iCenter_atom)
	{
	    // e2x
	    pAtom->xx  = ex*pAtom->ee + *pXx;
	    pAtom->xx += fx*pAtom->ff;
	    pAtom->xx += gx*pAtom->gg;
	    // f2y
	    pAtom->yy  = ey*pAtom->ee + *pYy;
	    pAtom->yy += fy*pAtom->ff;
	    pAtom->yy += gy*pAtom->gg;
	    // g2z
	    pAtom->zz  = ez*pAtom->ee + *pZz;
	    pAtom->zz += fz*pAtom->ff;
	    pAtom->zz += gz*pAtom->gg;
	}
    }
}

void _simmole::analyze_mole_structure()
{
    if (natom > _MAX_ATOM)
    {
	std::cout << "ERROR: Too many atoms in one molecule. Try increase _MAX_ATOM.\n";
	exit(1);
    }

    // check the neighbours for every atom
    int 	ii, jj;
    int		nNbr;
    PSIMATOM 	pAtom_1, pAtom_2;
    double	xxi, yyi, zzi;
    double	xij, yij, zij;
    double 	rij;
    for (ii=0;ii<natom;ii++)
    {
	nNbr = 0;
	pAtom_1 = atom_list + ii; // get the atom
	xxi = pAtom_1->xx;
	yyi = pAtom_1->yy;
	zzi = pAtom_1->zz;
	for (jj=0;jj<natom;jj++)
	{
	    if (ii==jj) continue;
	    pAtom_2 = atom_list + jj;
	    // calculate the distance
	    xij = xxi - pAtom_2->xx;
	    yij = yyi - pAtom_2->yy;
	    zij = zzi - pAtom_2->zz;
	    rij = sqrt(xij*xij + yij*yij + zij*zij);
	    if (rij <= std_bonding_dist)
	    {
		pAtom_1->pNbr_list[nNbr] = jj;
		nNbr++;
	    }
	    if (nNbr > _MAX_NBR)
	    {
		std::cout << "ERROR: Too many neighbours for one atom. Try increase _MAX_NBR.\n";
		exit(1);
	    }
	}
	pAtom_1->nnbr = nNbr;
    }
    // end of checking neighbours

    // check the bond, angle, dihedral
    PSIMATOM	pAtom_3, pAtom_4;
    int 	index1, index2, index3, index4;
    int		kk, ll;
    nbond = nangle = ndih = nimp = 0;
    for (ii=0;ii<natom;ii++)
    {
	pAtom_1 = atom_list + ii;
	index1 = pAtom_1->id;
	for (jj=0;jj<pAtom_1->nnbr;jj++)
	{
	    pAtom_2 = atom_list + pAtom_1->pNbr_list[jj]; // get a neighbour
	    index2 = pAtom_2->id;

	    if (index2 > index1) // avoid duplicate bonds
	    {
		pBond_list[nbond].comp_id[0] = index1;
		pBond_list[nbond].comp_id[1] = index2;
		nbond++;
		if (nbond > _MAX_BOND)
		{
		    std::cout << "ERROR: Too many bonds for one molecule. Try increase _MAX_BOND.\n";
		    exit(1);
		}

	    }

	    if (pAtom_2->nnbr > 1) // make sure pAtom_2 is not an end point
	    {
		for (kk=0;kk<pAtom_2->nnbr;kk++)
		{
		    pAtom_3 = atom_list + pAtom_2->pNbr_list[kk];
		    index3 = pAtom_3->id;

		    if (index3 > index1)
		    {
			pAngle_list[nangle].comp_id[0] = index1;
			pAngle_list[nangle].comp_id[1] = index2;
			pAngle_list[nangle].comp_id[2] = index3;
			nangle++;
			if (nangle > _MAX_ANGLE)
			{
			    std::cout << "ERROR: Too many angles for one molecule. Try increase _MAX_ANGLE.\n";
			    exit(1);
			}
		    }

		    if (index3 != index1 && pAtom_3->nnbr > 1) // make sure pAtom_3 is not an end point
		    {
			for (ll=0;ll<pAtom_3->nnbr;ll++)
			{
			    pAtom_4 = atom_list + pAtom_3->pNbr_list[ll];
			    index4 = pAtom_4->id;
			    if (index4 != index2 && index4 > index1)
			    {
				pDih_list[ndih].comp_id[0] = index1;
				pDih_list[ndih].comp_id[1] = index2;
				pDih_list[ndih].comp_id[2] = index3;
				pDih_list[ndih].comp_id[3] = index4;
				ndih++;
				if (ndih > _MAX_DIHEDRAL)
				{
				    std::cout << "ERROR: Too many dihedrals for one molecule. Try increase _MAX_DIHEDRAL.\n";
				    exit(1);
				}
			    }
			} // for (ll=0;ll<pAtom_3->nnbr;ll++)
		    } // if (pAtom_3->nnbr > 1)

		} // for (kk=0;kk<pAtom_2->nnbr;kk++)
	    } // if (pAtom_2->nnbr > 1)

	} // for (jj=0;jj<pAtom_1->nnbr;jj++)

    } // for (ii=0;ii<natom;ii++)
    // end of checking bond, angle, dihedral

    // check impropers
    // improper defined as 1-2-3-4
    // where 1 is the central atom which 2,3,4 bonded to
    // the angle is defined as the angle between two planes of
    // 1,2,3 and 2,3,4
    nimp = 0;
    for (ii=0;ii<natom;ii++)
    {
	pAtom_1 = atom_list + ii;
	index1 = pAtom_1->id;
	// if an atom has 3 or more neighbour
	// it is a center of an improper dihedral
	// atom 1
	if (pAtom_1->nnbr >= 3)
	{
	    // loop through the neighbour for atom 2
	    for (jj=0;jj<pAtom_1->nnbr-2;jj++)
	    {
		pAtom_2 = atom_list + pAtom_1->pNbr_list[jj]; // get atom 2
		index2 = pAtom_2->id;
		// loop through the neighbour for atom 3
		for (kk=jj+1;kk<pAtom_1->nnbr-1;kk++)
		{
		    pAtom_3 = atom_list + pAtom_1->pNbr_list[kk]; // get atom 3
		    index3 = pAtom_3->id;
		    // loop through the neighbour for atom 4
		    for (ll=kk+1;ll<pAtom_1->nnbr;ll++)
		    {
			pAtom_4 = atom_list + pAtom_1->pNbr_list[ll]; // get atom 4
		       	index4 = pAtom_4->id;

			pImp_list[nimp].comp_id[0] = index1;
			pImp_list[nimp].comp_id[1] = index2;
			pImp_list[nimp].comp_id[2] = index3;
			pImp_list[nimp].comp_id[3] = index4;
			nimp++;
			if (nimp > _MAX_IMPROPER)
			{
			    std::cout << "ERROR: Too many impropers for one molecule. Try increase _MAX_IMPROPER.\n";
			    exit(1);
			}
		    }
		}
	    }
	} // if an atom has three or more neighbours

    } // loop through all atoms for checking impropers


    // checking nonpair list
    nnonpair = 0;
    bool	isnonpair;
    PBOND	pBond;
    PANGLE	pAngle;
    PDIHEDRAL	pDih;
    for (ii=0;ii<natom-1;ii++)
    {
	for (jj=ii+1;jj<natom;jj++)
	{
	    isnonpair = true;

	    if (fNonpair_mode!=_1_2_NONPAIR) // 1,2  nonpair
	    {
		for (kk=0;kk<nbond;kk++) // bond checking
		{
		    pBond = pBond_list + kk;
		    if (ii==pBond->comp_id[0] && jj==pBond->comp_id[1])
			isnonpair = false;
		    if (ii==pBond->comp_id[1] && jj==pBond->comp_id[0])
			isnonpair = false;
		}
		if (fNonpair_mode!=_1_3_NONPAIR) // 1,3  nonpair
		{
		    for (kk=0;kk<nangle;kk++) // angle checking
		    {
			pAngle = pAngle_list + kk;
			if (ii==pAngle->comp_id[0] && jj==pAngle->comp_id[2])
			    isnonpair = false;
			if (ii==pAngle->comp_id[2] && jj==pAngle->comp_id[0])
			    isnonpair = false;
		    }
		    if (fNonpair_mode!=_1_4_NONPAIR) // 1,4 nonpair
		    {
			for (kk=0;kk<ndih;kk++) // dihedral checking
			{ 
			    pDih = pDih_list + kk;
			    if (ii==pDih->comp_id[0] && jj==pDih->comp_id[3])
				isnonpair = false;
			    if (ii==pDih->comp_id[3] && jj==pDih->comp_id[0])
				isnonpair = false;
			}
		    } // if (fNonpair_mode!=_1_4_NONPAIR) 
		} // if (fNonpair_mode!=_1_3_NONPAIR)
	    } // if (fNonpair_mode!=_1_2_NONPAIR)
	    if (isnonpair)
	    {
		pNonpair_list[nnonpair].comp_id[0] = ii;
		pNonpair_list[nnonpair].comp_id[1] = jj;
		nnonpair++;
		if (nnonpair > _MAX_NONPAIR)
		{
		    std::cout << "ERROR: Too many non-bonded pairs. Try increase _MAX_NONPAIR.\n";
		    exit(1);
		}
	    }

	} // for (jj=ii+1;jj<natom;jj++)
    } // for (ii=0;ii<natom-1;ii++)
}


void _simmole::print_structure()
{
    int ii, jj;
    PSIMATOM pAtom;
    PBOND pBond;
    PANGLE pAngle;
    PDIHEDRAL pDih;
    PIMPROPER pImp;
    PNONPAIR pNonpair;

    cout << "Number of atoms = " << natom << endl
	<< "Number of bonds = " << nbond << endl
	<< "Number of angles = " << nangle << endl
	<< "Number of dihedrals = " << ndih << endl
	<< "Number of nonpairs = " << nnonpair << endl;

    for (ii=0;ii<natom;ii++)
    {
	pAtom = atom_list + ii;
	cout << "Neighbours of atom " << ii << " = ";
	for (jj=0;jj<pAtom->nnbr;jj++)
	{
	    cout << pAtom->pNbr_list[jj];
	    if (jj<pAtom->nnbr-1)
		cout << ',';
	}
	cout << '\n';
    }

    for (ii=0;ii<nbond;ii++)
    {
	pBond = pBond_list + ii;
	cout << "Bond " << ii << " = "
	    << pBond->comp_id[0] << "-" 
	    << pBond->comp_id[1] << "    Pars: "
	    << pBond->kr << ", "
	    << pBond->req << ", "
	    << pBond->alpha << endl;
    }

    for (ii=0;ii<nangle;ii++)
    {
	pAngle = pAngle_list + ii;
	cout << "Angle " << ii << " = "
	    << pAngle->comp_id[0] << "-" 
	    << pAngle->comp_id[1] << "-" 
	    << pAngle->comp_id[2] << "    Pars: "
	    << pAngle->k_theta << ", "
	    << pAngle->theta_eq << ", "
	    << pAngle->agl_para_3 << ", "
	    << pAngle->agl_para_4 << ", "
	    << pAngle->agl_para_5 << endl;
    }

    for (ii=0;ii<ndih;ii++)
    {
	pDih = pDih_list + ii;
	cout << "Dihedral " << ii << " = "
	    << pDih->comp_id[0] << "-"
	    << pDih->comp_id[1] << "-"
	    << pDih->comp_id[2] << "-"
	    << pDih->comp_id[3] << "    Pars: "
	    << pDih->v1 << ", "
	    << pDih->v2 << ", "
	    << pDih->v3 << ", "
	    << pDih->v4 << endl;
    }

    for (ii=0;ii<nimp;ii++)
    {
	pImp = pImp_list + ii;
	cout << "Improper " << ii << " = "
	    << pImp->comp_id[0] << "-"
	    << pImp->comp_id[1] << "-"
	    << pImp->comp_id[2] << "-"
	    << pImp->comp_id[3] << "    Pars: "
	    << pImp->v1 << ", "
	    << pImp->v2 << ", "
	    << pImp->v3 << ", "
	    << pImp->v4 << endl;
    }



    for (ii=0;ii<nnonpair;ii++)
    {
	pNonpair = pNonpair_list + ii;
	cout << "Non-bonded pair " << ii << " = "
	    << pNonpair->comp_id[0] << ","
	    << pNonpair->comp_id[1] << endl;
    }


}





