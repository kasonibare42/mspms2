#ifndef DEFS_H_
#define DEFS_H_

#define INPUT		"in.mspms"
#define OUTPUT		"out.mspms"
#define CONFIG		"cfg.mspms"
#define COORDSIN	"xyz.mspms"
#define LOG			"log.mspms"
#define SNAPSHOT	"ss.mspms"
#define TRAJECTORY	"trj.mspms"
#define SAVE		"sav.mspms"
#define LOAD		"old.mspms"

#define MOLECULAR_DYNAMICS		0
#define HYBRID_MONTE_CARLO		1
#define SIMULATED_ANNEALING		2

#define NVE		0
#define NVT		1
#define NPT		2

#define DYNAMIC_ID			-1

#define GHOST_NONE	0
#define GHOST_LJ	1
#define GHOST_FULL	2

#define NEW			1
#define CONTINUE	2
#define CONFIG_ONLY	3

/// Inter-molecular potentials
#define INTER_MOLE_LJ	0

#define BOND_NONE			0
#define BOND_HARMONIC		1
#define BOND_MORSE			2
#define BOND_FENE			3
#define BOND_SPRING_PATH_INTEGRAL	4

#define ANGLE_NONE			0
#define ANGLE_HARMONIC		1
#define ANGLE_TR_WATER		2  ///< Toukan and Rahman water potentials

#define DIH_NONE			0
#define DIH_OPLS_COSIN		1
#define DIH_CHARMM			2 ///< charmm type dihedral potential

#define IMP_NONE		0
#define IMP_CHARMM		1 ///< charmm type improper potential

#define SF_NONE						0
#define SF_NANOTUBE_HYPERGEO  		1
#define SF_NANOTUBE_ATOM_EXPLICIT	2
#define SF_NANOTUBE_TASOS			3
#define SF_NANOTUBE_MY_INTERP		4 ///< my interpolation grid

#define ELECTROSTATIC_NONE				0
#define ELECTROSTATIC_EWALD				1
#define ELECTROSTATIC_WOLF				2
#define ELECTROSTATIC_SIMPLE_COULOMB	3

#define EWALD_BC_TINFOIL	1
#define EWALD_BC_VACUUM		2

#define EWALD_1D		1
#define EWALD_2D		2 ///< not yet in use
#define EWALD_3D		3

#define SOLID_HETERO		0  ///< solid type, e.g. mof?
#define SOLID_UNIFORM		1 ///< e.g. nanotoubes

#define FF_NONE					0
#define FF_DFT_METAL_CLUSTER	1

#define NCOUNTS_MAX	40
#define NUNIQUE_ATOM_MAX	5 ///< number of unqiue atoms for grid interpolation

// For numerical derivatives
#define X_AXIS	0
#define Y_AXIS	1
#define Z_AXIS	2

#define SHORT_STRING_LENGTH		5
#define LONG_STRING_LENGTH		200

#define STEP_SIZE	1.0e-8
#define TOLERANCE	1.0e-8
#define AVOGADRO	6.0221415e23    ///< mol^-1
#define BOLTZMAN_CONSTANT   1.3806505e-23  ///< J/K
#define PLANCK_CONSTANT		6.6260693e-34		///< J*s or (kg*m^2/s)
#define HBAR		1.05457168e-34		///< J*s or (kg*m^2/s)
#define COULOMB_PREFACTOR	167100.80840  ///< unit is J.A, (1.602177e-19)^2/(4*pi*epsilon0)/kB*1e10 
#define RGAS		8.314472 ///< J/mol/K 
#define R_RGAS		0.120272219 ///< reciprocal of RGAS
#define EV_TO_K		11604.5112770 /// ev to K. 1 ev = 1.6021773e-19 joule
#define pi		3.141592653589793
/// 3*pi^2*theta, assume the same charge density as graphite, theta=0.382 A^-2, used for hypergeo nanotube
#define C3_PI_SQ_THETA 	11.3105666436484 

#endif /*DEFS_H_*/
