/***************************************************************************
 *            atom.h
 *
 *  Mon May 23 2005
 *  Created by: Yang Wang
 ****************************************************************************/

#ifndef _ATOM_H
#define _ATOM_H

//---------------------------------------------------------
// struct of smooth nanotube related variables
// only useful for smooth nanotubes
// polynomial parameters
// A struct variable contain the interaction parameters
// for a specific atom type with various nanotubes(1,1)-(20,7).
//---------------------------------------------------------
// atom_type    the type name of a specific atom
// inalpha      the polynomial parameter for interactions
//              inside the nanotube. The first subscript is
//              the type of the nanotube, ie. (1,1)...
//              The second subscript is the order of the
//              parameters, 1-8.
// outalpha     same as inalpha, but is for the interactions
//              outside the nanotubes
//---------------------------------------------------------
// Check and set the variables in the variables initialization
// procedure, after read in or analysis of the atom type.
// Then read in following parameters from the potential database
typedef struct _alpha
{
    double          alpha_inner[_MAX_NANOTUBE_TYPE][_POLYNOMIAL_ORDER];
    double          alpha_outter[_MAX_NANOTUBE_TYPE][_POLYNOMIAL_ORDER];
} ALPHA, *PALPHA;

typedef class _simatom
{
    public:
	int             id;
	int             parent_mole_id; // the id of it's parent molecule
	char            name[20]; // H, C, O, etc.
	char            type[20]; // HC, CT, OS, etc
	double          weight; // kg/mol
	double          sigma; // angstrom
	double          epsilon; // J/mol
	int		fGhost; // flag for whether count the atom for LJ interactions or not
	int 		tasostype;
	double          charge;
	double		xx, yy, zz; // coordinates
	double		ee, ff, gg; // inner coordinates
	double		vx, vy, vz; // velocities
	double		px, py, pz; // momenta
	double		xdisp, ydisp, zdisp; // displacements
	double		fx, fy, fz; // total forces
	double		fxra, fyra, fzra; // total intra forces
	double		fxer, fyer, fzer; // total inter forces
	double		fx_lj, fy_lj, fz_lj; // inter-molecular LJ forces
	double		fx_nplj, fy_nplj, fz_nplj; // nonpair LJ forces
	double		fx_cl, fy_cl, fz_cl; // inter-molecular coulomb forces
	double		fx_npcl, fy_npcl, fz_npcl; // nonpair coulomb forces??
	double		fx_bond, fy_bond, fz_bond; // bond forces
	double		fx_angle, fy_angle, fz_angle; // angle forces
	double		fx_dih, fy_dih, fz_dih; // dihedral forces
	double		fx_sflj, fy_sflj, fz_sflj; // sf lj forces
	// Neighbours of an atom                       //
	int	        nnbr; // number of neighbours
	int*            pNbr_list; // neighbour list, pointer to a list of atoms
	// Variables for nanotubes only                //
	int		fPosition; // where is the atom, outside, inside, etc.
	bool            bOutside; // If the atom is outside nanotubes.
	PALPHA          pAlpha_list; // alpha values for smooth nanotubes

	// Functions
	_simatom();
	~_simatom();
	void 			check_atom_weight(); // check atom's MW
	void			check_atom_props(); // check atom's charge, sigma, epsilon, etc.
	void			check_atom_alpha(); // check smooth nanotube parameters for this atom

} SIMATOM, *PSIMATOM;

#endif /* _ATOM_H */

