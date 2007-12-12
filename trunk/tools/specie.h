/***************************************************************************
 *            specie.h
 *
 *  Sun Jul 17 2005
 *  Created by: Yang Wang
 ****************************************************************************/
#ifndef _SPECIE_H
#define _SPECIE_H   

typedef struct _accumu
{
    long double acc1, acc2, acc3, acc4, acc5;
}ACCUMU, *PACCUMU;

// specie is a collection of molecules of the same type
typedef class _simspecie
{
    public:
	int 		id;
	int 		nmole;
	int		nmole_max;
	int 		natom;
	int		natom_per_mole;
	char 		name[_MAX_NAME_LENGTH];

	double		zact;

	// structure variables
	int		nbond; // for molecules
	int		nangle; // for molecules
	int		ndih; // for molecules
	int		nnonpair; // for molecules
	int		nnbr[_MAX_ATOM]; // number of neighbours for atoms

	BOND		bond_list[_MAX_BOND];
	ANGLE		angle_list[_MAX_ANGLE];
	DIHEDRAL	dih_list[_MAX_DIHEDRAL];
	IMPROPER	imp_list[_MAX_IMPROPER];
	NONPAIR		nonpair_list[_MAX_NONPAIR];
	int		nbr_list[_MAX_ATOM][_MAX_NBR];

	// alpha parameters for smooth nanotubes
	PALPHA		alpha_list;

	// molecule list
	PSIMMOLE	mole_list;

	// sample molecule
	SIMMOLE		sample_mole;

	// counters and accumulators
	// counter[0] insertions
	// counter[1] accepted insertions
	// counter[2] deletions
	// counter[3] accepted deletions
	long		counter[_MAX_COUNTER_SPECIE];
	// 0 nmole
	ACCUMU		accumu[_MAX_ACCUMU_SPECIE];

	// output data
	// natom_ave, natom_err, natom_fluc, etc.;

	_simspecie();
	~_simspecie();
	void init();

} SIMSPECIE, *PSIMSPECIE;

#endif /* SPECIE_H */
