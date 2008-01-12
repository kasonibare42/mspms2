#ifndef MOLEVARS_H_
#define MOLEVARS_H_

// definitions
// atom, molecule, specie, bond, angle, dihedral, improper, non-bonded pair		
#define NATOM_MAX	4000
#define NMOLE_MAX	1000
#define NSPECIE_MAX	3
#define NBOND_MAX	2000
#define NANGLE_MAX	2000
#define NDIH_MAX	2000
#define NIMP_MAX	2000
#define NNBP_MAX	4000

#define MOLE_STATUS_UNINIT	-1
#define MOLE_STATUS_VACANCY	0
#define MOLE_STATUS_NORMAL	1

// Sample species/molecules
#define SAMPLE_NATOM_MAX	200
#define SAMPLE_NBOND_MAX	100
#define SAMPLE_NANGLE_MAX	100
#define SAMPLE_NDIH_MAX		100
#define SAMPLE_NIMP_MAX		100
#define SAMPLE_NNBP_MAX		100

#define EXCLUDE_LIST_MAX	1000
#define SAMPLE_NATOM_PER_MOLE_MAX	100

// variables
// atoms
int natom; ///< total number of atoms in the system
char atomname[NATOM_MAX][5];
int atom2mole[NATOM_MAX]; ///< which molecule this atom belongs to
int isghost[NATOM_MAX]; ///< flag of ghost atom, including ghost LJ or ghost electrostatic
int tasostype[NATOM_MAX]; ///< atom type for tasos interpolation 
double aw[NATOM_MAX], epsilon[NATOM_MAX], sigma[NATOM_MAX], charge[NATOM_MAX]; ///< atom weight etc. 
double xx[NATOM_MAX], yy[NATOM_MAX], zz[NATOM_MAX]; ///< atom coordinates
double ex[NATOM_MAX], fy[NATOM_MAX], gz[NATOM_MAX];///< inner coordinates relative to the center of mass, also used for PBC reconstruction of the molecule
double vx[NATOM_MAX], vy[NATOM_MAX], vz[NATOM_MAX]; ///< atom velocity
double fxl[NATOM_MAX], fyl[NATOM_MAX], fzl[NATOM_MAX]; ///< inter force
double fxs[NATOM_MAX], fys[NATOM_MAX], fzs[NATOM_MAX]; ///< intra force

// molecule
int nmole; ///< total number of molecules in the system
int iPhysicalMoleIDFromMetaID[NMOLE_MAX]; ///< get the physical ID of a molecule using its meta ID
int iPhysicalMoleIDFromMetaIDinSpecie[NSPECIE_MAX][NMOLE_MAX]; ///< get the physical ID of a molecule using its meta ID within its specie
int mole2specie[NMOLE_MAX]; ///< which specie this molecule belongs to
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
int nspecie; ///< total number of species in the system
char szSpecieName[NSPECIE_MAX][200]; ///< the name of each specie
int nmole_per_specie[NSPECIE_MAX]; ///< number of molecules in a certain specie
int specie_first_atom_idx[NSPECIE_MAX];
int specie_last_atom_idx[NSPECIE_MAX];
int specie_first_mole_idx[NSPECIE_MAX];
int specie_last_mole_idx[NSPECIE_MAX];
int specie_first_vacancy_idx[NSPECIE_MAX]; ///< the index of first vancant molecule for a specie

// System
char sysname[200]; ///< name of the object system

// bond, angle, dihedral, improper, non-bonded pair				
int nbond, nangle, ndih, nimp, nnbp; ///< bond
int bond_idx[NBOND_MAX][2];
int bond_type[NBOND_MAX];
double Kb[NBOND_MAX], Req[NBOND_MAX], alpha[NBOND_MAX];
int angle_idx[NANGLE_MAX][3]; ///< angle
int angle_type[NANGLE_MAX];
double Ktheta[NANGLE_MAX], Thetaeq[NANGLE_MAX];
double agl_para_3[NANGLE_MAX], agl_para_4[NANGLE_MAX], agl_para_5[NANGLE_MAX];
int isAngle_unique[NANGLE_MAX]; ///< make sure 1,3 atoms do not make mutiple angles, this is for possible ring structures
int dih_idx[NDIH_MAX][4]; ///< dihedral
int dih_type[NDIH_MAX];
double c1[NDIH_MAX], c2[NDIH_MAX], c3[NDIH_MAX], c4[NDIH_MAX];
int isDih_unique[NDIH_MAX]; ///< make sure 1,4 atoms do not make mutiple dihedrals, this is for possible ring structures
int imp_idx[NIMP_MAX][4]; ///< improper
int imp_type[NIMP_MAX];
double komega[NIMP_MAX], omega0[NIMP_MAX];
int nbp_idx[NNBP_MAX][2]; ///< non-bonded pairs

// Sample species/molecules				
double sample_mw[NSPECIE_MAX]; ///< molcule weight for Samples
int sample_natom_per_mole[NSPECIE_MAX]; ///< number of atoms in a Sample
int sample_nbond_per_mole[NSPECIE_MAX];
int sample_nangle_per_mole[NSPECIE_MAX];
int sample_ndih_per_mole[NSPECIE_MAX];
int sample_nimp_per_mole[NSPECIE_MAX];
int sample_nnbp_per_mole[NSPECIE_MAX];
char sample_atomname[SAMPLE_NATOM_MAX][5];
double sample_xx[SAMPLE_NATOM_MAX];
double sample_yy[SAMPLE_NATOM_MAX];
double sample_zz[SAMPLE_NATOM_MAX];
double sample_ee[SAMPLE_NATOM_MAX];
double sample_ff[SAMPLE_NATOM_MAX];
double sample_gg[SAMPLE_NATOM_MAX];
double sample_aw[SAMPLE_NATOM_MAX]; 
double sample_epsilon[SAMPLE_NATOM_MAX]; 
double sample_sigma[SAMPLE_NATOM_MAX]; 
double sample_charge[SAMPLE_NATOM_MAX]; 
double sample_isghost[SAMPLE_NATOM_MAX]; 
double sample_tasostype[SAMPLE_NATOM_MAX];
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
int sample_bond_idx[SAMPLE_NBOND_MAX][2];
int sample_angle_idx[SAMPLE_NANGLE_MAX][3];
int sample_dih_idx[SAMPLE_NDIH_MAX][4];
int sample_imp_idx[SAMPLE_NIMP_MAX][4];
int sample_nbp_idx[SAMPLE_NNBP_MAX][2];
int sample_bond_type[SAMPLE_NBOND_MAX];
double sample_Kb[SAMPLE_NBOND_MAX];
double sample_Req[SAMPLE_NBOND_MAX];
double sample_alpha[SAMPLE_NBOND_MAX];
int sample_angle_type[SAMPLE_NANGLE_MAX];
double sample_Ktheta[SAMPLE_NANGLE_MAX];
double sample_Thetaeq[SAMPLE_NANGLE_MAX];
double sample_agl_para_3[SAMPLE_NANGLE_MAX];
double sample_agl_para_4[SAMPLE_NANGLE_MAX];
double sample_agl_para_5[SAMPLE_NANGLE_MAX];
int sample_isAngle_unique[SAMPLE_NANGLE_MAX];
int sample_dih_type[SAMPLE_NDIH_MAX];
double sample_c1[SAMPLE_NDIH_MAX];
double sample_c2[SAMPLE_NDIH_MAX];
double sample_c3[SAMPLE_NDIH_MAX];
double sample_c4[SAMPLE_NDIH_MAX];
int sample_isDih_unique[SAMPLE_NDIH_MAX];
int sample_imp_type[SAMPLE_NIMP_MAX];
double sample_komega[SAMPLE_NIMP_MAX];
double sample_omega0[SAMPLE_NIMP_MAX];
int excllist[EXCLUDE_LIST_MAX]; ///< The list of excluding pair of atoms
int pointexcl_atom[NSPECIE_MAX][SAMPLE_NATOM_PER_MOLE_MAX+1]; ///< the index of exclude list starts and ends for an atom


#endif /*MOLEVARS_H_*/
