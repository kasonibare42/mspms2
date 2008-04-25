#ifndef MOLEVARS_H_
#define MOLEVARS_H_

// Definitions
// atom, molecule, specie, bond, angle, dihedral, improper, non-bonded pair		
#define NATOM_MAX	4000
#define NBOND_MAX	2000
#define NANGLE_MAX	2000
#define NDIH_MAX	2000
#define NIMP_MAX	2000
#define NNBP_MAX	4000
// Molecules
#define NMOLE_MAX	1000
// Species
#define NSPECIE_MAX	3

// Sample species/molecules
#define SAMPLE_NATOM_MAX	200
#define SAMPLE_NBOND_MAX	100
#define SAMPLE_NANGLE_MAX	100
#define SAMPLE_NDIH_MAX		100
#define SAMPLE_NIMP_MAX		100
#define SAMPLE_NNBP_MAX		100

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

typedef struct _SAMPLE_MOLECULE
{
    // Atom related properties
	char atom_name[SAMPLE_NATOM_MAX][SHORT_STRING_LENGTH];
	int type[SAMPLE_NATOM_MAX], ghost_type[SAMPLE_NATOM_MAX], 
		interp_type[SAMPLE_NATOM_MAX];
	double aw[SAMPLE_NATOM_MAX], epsilon[SAMPLE_NATOM_MAX], 
	       sigma[SAMPLE_NATOM_MAX], charge[SAMPLE_NATOM_MAX];
	double xx[SAMPLE_NATOM_MAX], yy[SAMPLE_NATOM_MAX], zz[SAMPLE_NATOM_MAX];
	double ee[SAMPLE_NATOM_MAX], ff[SAMPLE_NATOM_MAX], gg[SAMPLE_NATOM_MAX];
	// Molecule related properties
	char mole_name[LONG_STRING_LENGTH];
	double mw; // molecular weight
	int natom, nbond, nangle, ndih, nimp, nnbp;
	int bnd_idx[SAMPLE_NBOND_MAX][2];
	int bnd_type[SAMPLE_NBOND_MAX];
	double bnd_para_0[SAMPLE_NBOND_MAX], bnd_para_1[SAMPLE_NBOND_MAX], 
	       bnd_para_2[SAMPLE_NBOND_MAX];
	int agl_idx[SAMPLE_NANGLE_MAX][3];
	int agl_type[SAMPLE_NANGLE_MAX];
	double agl_para_1[SAMPLE_NANGLE_MAX], agl_para_2[SAMPLE_NANGLE_MAX], 
	       agl_para_3[SAMPLE_NANGLE_MAX], agl_para_4[SAMPLE_NANGLE_MAX], 
	       agl_para_5[SAMPLE_NANGLE_MAX];
	bool isAngle_unique[SAMPLE_NANGLE_MAX];
	int dih_idx[SAMPLE_NDIH_MAX][4];
	int dih_type[SAMPLE_NDIH_MAX];
	double dih_para_1[SAMPLE_NDIH_MAX], dih_para_2[SAMPLE_NDIH_MAX], 
	       dih_para_3[SAMPLE_NDIH_MAX], dih_para_4[SAMPLE_NDIH_MAX];
	bool isDih_unique[SAMPLE_NDIH_MAX];
	int imp_idx[SAMPLE_NIMP_MAX][4];
	int imp_type[SAMPLE_NIMP_MAX];
	double imp_para_0[SAMPLE_NIMP_MAX], imp_para_1[SAMPLE_NIMP_MAX];
	double nbp_idx[SAMPLE_NNBP_MAX][2];
} SAMPLE_MOLECULE, *PSAMPLE_MOLECULE;

SAMPLE_MOLECULE sample_mole[NSPECIE_MAX]; // sample molecules

// atoms
int natom; ///< total number of atoms in the system
double xx[NATOM_MAX], yy[NATOM_MAX], zz[NATOM_MAX]; ///< atom coordinates
double ex[NATOM_MAX], fy[NATOM_MAX], gz[NATOM_MAX];///< inner coordinates relative to the center of mass, also used for PBC reconstruction of the molecule
double vx[NATOM_MAX], vy[NATOM_MAX], vz[NATOM_MAX]; ///< atom velocity
double fxl[NATOM_MAX], fyl[NATOM_MAX], fzl[NATOM_MAX]; ///< inter force
double fxs[NATOM_MAX], fys[NATOM_MAX], fzs[NATOM_MAX]; ///< intra force

// molecules
int nmole; ///< total number of molecules in the system
double mole_xx[NMOLE_MAX], mole_yy[NMOLE_MAX], mole_zz[NMOLE_MAX];
double fOrientE_x, fOrientE_y, fOrientE_z; ///< the E axis of the orientation vectors
double fOrientF_x, fOrientF_y, fOrientF_z;
double fOrientG_x, fOrientG_y, fOrientG_z;

// specie
// All the index number are based on the entire molecule list (not restricted to a single specie)
int nspecie; ///< total number of species in the system
int nmole_per_specie[NSPECIE_MAX]; ///< number of molecules in a certain specie
int natom_per_specie[NSPECIE_MAX];

// System
char sysname[200]; ///< name of the object system
double system_mass; ///< mass of the system, unit is kg/mol

// bond, angle, dihedral, improper, non-bonded pair				
int nbond, nangle, ndih, nimp, nnbp; ///< bond

#endif /*MOLEVARS_H_*/
