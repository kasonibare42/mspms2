#ifndef MOLEVARS_H_
#define MOLEVARS_H_

// definitions
/// ATOM_LIST_LENGTH = NATOM_MAX + SAMPLE_NATOM_MAX 
#define ATOM_LIST_LENGTH			4200
/// BOND_LIST_LENGTH = NBOND_MAX + SAMPLE_NBOND_MAX
#define BOND_LIST_LENGTH			2100
#define ANGLE_LIST_LENGTH			2100
#define DIHEDRAL_LIST_LENGTH		2100
#define IMPROPER_LIST_LENGTH		2100
#define NBP_LIST_LENGTH				4100

// atom, molecule, specie, bond, angle, dihedral, improper, non-bonded pair		
#define NATOM_MAX	4000
#define NBOND_MAX	2000
#define NANGLE_MAX	2000
#define NDIH_MAX	2000
#define NIMP_MAX	2000
#define NNBP_MAX	4000

// Sample species/molecules
#define SAMPLE_NATOM_MAX	200
#define SAMPLE_NBOND_MAX	100
#define SAMPLE_NANGLE_MAX	100
#define SAMPLE_NDIH_MAX		100
#define SAMPLE_NIMP_MAX		100
#define SAMPLE_NNBP_MAX		100

// Molecules
#define NMOLE_MAX	1000
// Species
#define NSPECIE_MAX	3

#define MOLE_STATUS_UNINIT	-1
#define MOLE_STATUS_VACANCY	0
#define MOLE_STATUS_NORMAL	1

#define EXCLUDE_LIST_MAX	1000
#define SAMPLE_NATOM_PER_MOLE_MAX	100

/// The index of the first atom which is NOT initialized (used by any molecules)
#define idAtomUninit	natom_hist_max 
#define idBondUninit	nbond_hist_max
#define idAngleUninit	nangle_hist_max
#define idDihUninit		ndih_hist_max
#define idImpUninit		nimp_hist_max
#define idNbpUninit		nnbp_hist_max

#define Kb			bnd_para_0
#define Req			bnd_para_1
#define alpha		bnd_para_2

#define Ktheta			agl_para_1
#define Thetaeq			agl_para_2

#define c1			dih_para_1
#define c2			dih_para_2
#define c3			dih_para_3
#define c4			dih_para_4

#define komega		imp_para_0
#define omega0		imp_para_1

#define sample_atomname	(&atomname[NATOM_MAX])
#define sample_xx		(xx+NATOM_MAX)
#define sample_yy		(yy+NATOM_MAX)
#define sample_zz		(zz+NATOM_MAX)
#define sample_ee		(ex+NATOM_MAX)
#define sample_ff		(fy+NATOM_MAX)
#define sample_gg		(gz+NATOM_MAX)
#define sample_aw		(aw+NATOM_MAX)
#define sample_epsilon	(epsilon+NATOM_MAX)
#define sample_sigma	(sigma+NATOM_MAX)
#define sample_charge	(charge+NATOM_MAX)
#define sample_ghost_type	(ghost_type+NATOM_MAX)
#define sample_tasos_type 	(tasos_type+NATOM_MAX)

#define sample_bond_idx		(&bond_idx[NBOND_MAX])
#define sample_bond_type	(bond_type+NBOND_MAX)
#define sample_Kb			(Kb+NBOND_MAX)
#define sample_Req			(Req+NBOND_MAX)
#define sample_alpha		(alpha+NBOND_MAX)

#define sample_angle_idx 	(&angle_idx[NANGLE_MAX])
#define sample_angle_type 	(angle_type+NANGLE_MAX)
#define sample_Ktheta 		(Ktheta+NANGLE_MAX)
#define sample_Thetaeq 		(Thetaeq+NANGLE_MAX)
#define sample_agl_para_3 	(agl_para_3+NANGLE_MAX)
#define sample_agl_para_4 	(agl_para_4+NANGLE_MAX)
#define sample_agl_para_5 	(agl_para_5+NANGLE_MAX)
#define sample_isAngle_unique 	(isAngle_unique+NANGLE_MAX)

#define sample_dih_idx 		(&dih_idx[NDIH_MAX])
#define sample_dih_type 	(dih_type+NDIH_MAX)
#define sample_c1 		(c1+NDIH_MAX)
#define sample_c2 		(c2+NDIH_MAX)
#define sample_c3  		(c3+NDIH_MAX)
#define sample_c4 		(c4+NDIH_MAX)
#define sample_isDih_unique 	(isDih_unique+NDIH_MAX)

#define sample_imp_idx 		(&imp_idx[NIMP_MAX])
#define sample_imp_type 	(imp_type+NIMP_MAX)
#define sample_komega 		(komega+NIMP_MAX)
#define sample_omega0 		(omega0+NIMP_MAX)

#define sample_nbp_idx 		(&nbp_idx[NNBP_MAX])

// atoms
int natom; ///< total number of atoms in the system
/** max total number of atoms ever achieved in the system
 * it is used to loop through the atom list, since vacant molecules could be in the middle
 * of the molecule list. Hence their atoms could be in the middle of the atom list too.
 * To loop through the entire atom list, we can not just use natom.
 * OOOOOO..OOOO..OOOO     (O is valid atom, . is vacancy)
 * so natom = 14, natom_hist_max = 18
 * To loop through the entire atom list, we have to use 18 and skip the vacant atoms.
 * Skipping the vacant atoms, this is done by using mole_status variable. Find out which
 * molecule the atom belongs to and then check the molecule's status. Same for bond, angle,
 * dihedral, improper, non-bonded pair.
 */
int natom_hist_max, nbond_hist_max, nangle_hist_max, ndih_hist_max,
		nimp_hist_max, nnbp_hist_max, nmole_hist_max;

char atomname[ATOM_LIST_LENGTH][5];
int ghost_type[ATOM_LIST_LENGTH]; ///< flag of ghost atom, including ghost LJ or ghost electrostatic
int tasos_type[ATOM_LIST_LENGTH]; ///< atom type for tasos interpolation 
double aw[ATOM_LIST_LENGTH], epsilon[ATOM_LIST_LENGTH],
		sigma[ATOM_LIST_LENGTH], charge[ATOM_LIST_LENGTH]; ///< atom weight etc. 
double xx[ATOM_LIST_LENGTH], yy[ATOM_LIST_LENGTH], zz[ATOM_LIST_LENGTH]; ///< atom coordinates
double ex[ATOM_LIST_LENGTH], fy[ATOM_LIST_LENGTH], gz[ATOM_LIST_LENGTH];///< inner coordinates relative to the center of mass, also used for PBC reconstruction of the molecule
double vx[NATOM_MAX], vy[NATOM_MAX], vz[NATOM_MAX]; ///< atom velocity
double fxl[NATOM_MAX], fyl[NATOM_MAX], fzl[NATOM_MAX]; ///< inter force
double fxs[NATOM_MAX], fys[NATOM_MAX], fzs[NATOM_MAX]; ///< intra force

/// which molecule (the physical ID) this atom belongs to. Physical atom id to physical mole id.
int atom2mole[NATOM_MAX], bnd2mole[NBOND_MAX], agl2mole[NANGLE_MAX],
		dih2mole[NDIH_MAX], imp2mole[NIMP_MAX], nbp2mole[NNBP_MAX];

// molecule
int mole2specie[NMOLE_MAX]; ///< which specie this molecule belongs to
/// Get the physical ID of a molecule using its meta ID within its specie.
int iMIDS_to_PMID[NSPECIE_MAX][NMOLE_MAX];

int nmole; ///< total number of molecules in the system
double mw[NMOLE_MAX]; ///< molecule weight
int mole_status[NMOLE_MAX]; ///< the status of the molecule, e.g. vacancy etc.
double mole_xx[NMOLE_MAX], mole_yy[NMOLE_MAX], mole_zz[NMOLE_MAX];
int mole_first_atom_idx[NMOLE_MAX]; ///< index of the first atom in a molecule
int mole_first_bond_idx[NMOLE_MAX]; ///< index of the first bond in a molecule
int mole_first_angle_idx[NMOLE_MAX]; ///< index of the first angle in a molecule
int mole_first_dih_idx[NMOLE_MAX]; ///< index of the first dihedral in a molecule
int mole_first_imp_idx[NMOLE_MAX]; ///< index of the first improper in a molecule
int mole_first_nbp_idx[NMOLE_MAX]; ///< index of the first nonbonded pair in a molecule
int mole_last_atom_idx[NMOLE_MAX]; ///< index of the last atom in a molecule
int mole_last_bond_idx[NMOLE_MAX]; ///< index of the last bond in a molecule
int mole_last_angle_idx[NMOLE_MAX]; ///< index of the last angle in a molecule
int mole_last_dih_idx[NMOLE_MAX]; ///< index of the last dihedral in a molecule
int mole_last_imp_idx[NMOLE_MAX]; ///< index of the last improper in a molecule
int mole_last_nbp_idx[NMOLE_MAX]; ///< index of the last nonbonded pair in a molecule
double fOrientE_x, fOrientE_y, fOrientE_z; ///< the E axis of the orientation vectors
double fOrientF_x, fOrientF_y, fOrientF_z;
double fOrientG_x, fOrientG_y, fOrientG_z;

// specie
// All the index number are based on the entire molecule list (not restricted to a single specie)
int nspecie; ///< total number of species in the system
char szSpecieName[NSPECIE_MAX][200]; ///< the name of each specie
int nmole_per_specie[NSPECIE_MAX]; ///< number of molecules in a certain specie
int natom_per_specie[NSPECIE_MAX];
/// The physical index of first vancant molecule for a specie
int specie_first_vacancy_idx[NSPECIE_MAX];
/// The first, last atom and molecule info for a specie. May NOT be useful since the molecule list for a specie can be discontinuous.
int specie_first_atom_idx[NSPECIE_MAX];
int specie_last_atom_idx[NSPECIE_MAX];
int specie_first_mole_idx[NSPECIE_MAX];
int specie_last_mole_idx[NSPECIE_MAX];

// System
char sysname[200]; ///< name of the object system
double system_mass; ///< mass of the system, unit is kg/mol

// bond, angle, dihedral, improper, non-bonded pair				
int nbond, nangle, ndih, nimp, nnbp; ///< bond
int bond_idx[BOND_LIST_LENGTH][2];
int bond_type[BOND_LIST_LENGTH];
double bnd_para_0[BOND_LIST_LENGTH], bnd_para_1[BOND_LIST_LENGTH],
		bnd_para_2[BOND_LIST_LENGTH];

int angle_idx[ANGLE_LIST_LENGTH][3]; ///< angle
int angle_type[ANGLE_LIST_LENGTH];
double agl_para_1[ANGLE_LIST_LENGTH], agl_para_2[ANGLE_LIST_LENGTH],
		agl_para_3[ANGLE_LIST_LENGTH], agl_para_4[ANGLE_LIST_LENGTH],
		agl_para_5[ANGLE_LIST_LENGTH];
bool isAngle_unique[ANGLE_LIST_LENGTH]; ///< make sure 1,3 atoms do not make mutiple angles, this is for possible ring structures

int dih_idx[DIHEDRAL_LIST_LENGTH][4]; ///< dihedral
int dih_type[DIHEDRAL_LIST_LENGTH];
double dih_para_1[DIHEDRAL_LIST_LENGTH], dih_para_2[DIHEDRAL_LIST_LENGTH],
		dih_para_3[DIHEDRAL_LIST_LENGTH], dih_para_4[DIHEDRAL_LIST_LENGTH];
bool isDih_unique[DIHEDRAL_LIST_LENGTH]; ///< make sure 1,4 atoms do not make mutiple dihedrals, this is for possible ring structures

int imp_idx[IMPROPER_LIST_LENGTH][4]; ///< improper
int imp_type[IMPROPER_LIST_LENGTH];
double imp_para_0[IMPROPER_LIST_LENGTH], imp_para_1[IMPROPER_LIST_LENGTH];

int nbp_idx[NBP_LIST_LENGTH][2]; ///< non-bonded pairs

// Sample species/molecules				
double sample_mw[NSPECIE_MAX]; ///< molcule weight for Samples
int sample_natom_per_mole[NSPECIE_MAX]; ///< number of atoms in a Sample
int sample_nbond_per_mole[NSPECIE_MAX];
int sample_nangle_per_mole[NSPECIE_MAX];
int sample_ndih_per_mole[NSPECIE_MAX];
int sample_nimp_per_mole[NSPECIE_MAX];
int sample_nnbp_per_mole[NSPECIE_MAX];

int sample_mole_first_atom_idx[NSPECIE_MAX];
int sample_mole_first_bond_idx[NSPECIE_MAX];
int sample_mole_first_angle_idx[NSPECIE_MAX];
int sample_mole_first_dih_idx[NSPECIE_MAX];
int sample_mole_first_imp_idx[NSPECIE_MAX];
int sample_mole_first_nbp_idx[NSPECIE_MAX];
int sample_mole_last_atom_idx[NSPECIE_MAX];
int sample_mole_last_bond_idx[NSPECIE_MAX];
int sample_mole_last_angle_idx[NSPECIE_MAX];
int sample_mole_last_dih_idx[NSPECIE_MAX];
int sample_mole_last_imp_idx[NSPECIE_MAX];
int sample_mole_last_nbp_idx[NSPECIE_MAX];

int excllist[EXCLUDE_LIST_MAX]; ///< The list of excluding pair of atoms
int pointexcl_atom[NSPECIE_MAX][SAMPLE_NATOM_PER_MOLE_MAX+1]; ///< the index of exclude list starts and ends for an atom


#endif /*MOLEVARS_H_*/
