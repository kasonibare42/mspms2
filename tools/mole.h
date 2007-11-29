/***************************************************************************
 *            mole.h
 *
 *  Fri July 15 2005
 *  Created by: Yang Wang
 ****************************************************************************/

#ifndef _MOLE_H
#define _MOLE_H

typedef struct _bond
{
    int		comp_id[2]; // the id of the two atoms of the bond
    double	kr, req; // bond parameters for harmonic bonding
    double	alpha; // Morse potential, De is kr
} BOND, *PBOND;

typedef struct _angle
{
    int		comp_id[3]; // the id of the three atoms of the angle
    double	k_theta, theta_eq; // angle parameters
    // following parameters are used for Toukan and Rahman's water angle potential
    // Toukan, K.;Rahman, A. Phys. Rev. B 1985, 31, 2643
    // Matej Praprotnik et. al. JPCA 2004, 108, 11056
    double	agl_para_3, agl_para_4, agl_para_5;
    double	k_r_theta, k_r_rprime; 
    double	req_HH, req_OH;
} ANGLE, *PANGLE;

typedef struct _dihedral
{
    int		comp_id[4];
    double	v1, v2, v3, v4;
} DIHEDRAL, *PDIHEDRAL;

typedef struct _improper
{
    int		comp_id[4];
    double	v1, v2, v3, v4;
} IMPROPER, *PIMPROPER;

typedef struct _nonpair
{
    int		comp_id[2];
    double	fij;
} NONPAIR, *PNONPAIR;

typedef class _simmole
{
    public:
	int		natom;
	int		id;
	double		weight; // kg/mol
	double		radius; // angstrom
	double		*pXx, *pYy, *pZz; // coordinates of a certain atom, used as the center of the molecule
	double		xx, yy, zz; // can still be used as center of mass if the center is not a specific atom
	int		iCenter_atom;
	// Variables for nanotube systems only                //
	int		fPosition; // where is the molecule, outside, inside, etc., only for nanotube system

	// atom list
	PSIMATOM	atom_list;

	// structural variables
	int		nbond;
	int		nangle;
	int		ndih;
	int		nimp;
	int		nnonpair;
	PBOND		pBond_list;
	PANGLE		pAngle_list;
	PDIHEDRAL	pDih_list;
	PIMPROPER	pImp_list;
	PNONPAIR	pNonpair_list;

	// sample molecule
	_simmole*	pSample_mole;

	// variables used for analyze molecular structure
	double		std_bonding_dist;
	int		fNonpair_mode; // 1,4 nonpair or 1,5 nonpair

	// Functions
	_simmole();
	~_simmole();

	// calculate the xyz coordinates to efg
	void		cal_all_atom_xyz2efg();
	// calculate the efg coordinates to xyz
	void		cal_all_atom_efg2xyz();
	// calculate the efg coordinates to xyz according to orientations
	void		cal_all_atom_efg2xyz_ortn(double, double, double, double, double, double, double, double, double);
	// analyze the structure of the molecule
	void		analyze_mole_structure();
	// calculate the center of mass
	void 		cal_com();
	// print structures
	void 		print_structure();
} SIMMOLE, *PSIMMOLE;

#endif
