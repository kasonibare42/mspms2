#ifndef VARS_H
#define VARS_H   

#define _safealloc(pt,num,size)         pt=calloc(num,size); assert(pt!=NULL)
#define _safefree(pt)	if ((pt)!=NULL) {free(pt); (pt)=NULL;}

#define INPUT		"in.mspms"
#define OUTPUT		"out.mspms"
#define CONFIG		"cfg.mspms"
#define COORDSIN	"xyz.mspms"
#define LOGFILE		"log.mspms"
#define SNAPSHOT	"ss.mspms"
#define MOVIE		"trj.mspms"
#define SAVEFILE	"sav.mspms"
#define LOADFILE	"old.mspms"

#define true		1
#define false		0

#define md_run		0
#define hmc_run		1

#define nve_run		0
#define nvt_run		1
#define npt_run		2

#define not_ghost	0
#define lj_ghost	1
#define all_ghost	2

#define new_run		1
#define continue_run	2
#define new_from_old	3

#define bond_none		0
#define bond_harmonic		1
#define bond_morse		2
#define bond_fene		3

#define angle_none		0
#define angle_harmonic		1
#define angle_TRwater		2  ///< Toukan and Rahman water potentials
#define dih_none		0
#define dih_opls_cosin		1
#define dih_charmm		2 ///< charmm type dihedral potential
#define imp_none		0
#define imp_charmm		1 ///< charmm type improper potential
#define _NO_SF_POTENTIAL	0
#define nanotube_hypergeo  	1
#define nanotube_atom_explicit	2
#define nanotube_tasos		3
#define nanotube_my_interp	4 ///< my interpolation grid
#define _NO_ELECTROSTATIC_INTERACTION	0
#define elec_ewald		1
#define elec_wolf		2
#define elec_simple_coulomb	3

#define ewald_bc_tinfoil	1
#define ewald_bc_vacuum		2

#define ewald_1D		1
#define ewald_2D		2 ///< not yet in use
#define ewald_3D		3

#define solid_hetero		0  ///< solid type, e.g. mof?
#define solid_uniform		1 ///< e.g. nanotoubes
#define nspecie_max	3
#define nmole_max	1000
#define natom_max	4000
#define nbond_max	2000
#define nangle_max	2000
#define ndih_max	2000
#define nimp_max	2000
#define nnbp_max	4000
#define exclude_max	4000
#define natom_per_mole_max	100
#define num_counter_max	40
#define nunique_atom_max	5 ///< number of unqiue atoms for grid interpolation
#define Rgas		8.314472 ///< J/mol/K 
#define rRgas		0.120272219 ///< reciprocal of Rgas
#define Avogadro	6.0221415e23    ///< mol^-1
#define const_columb	1389355.1051  ///< unit is J/mol. Na*(1.602177e-19)^2/4/PI/epsilon0/(1.0e-10)
#define kb_1e30		1.3806505e7   ///< kb/1e-30  
#define virial_to_pressure	553512.954376   ///< 1.0/(3.0*6.0221415e-7)
#define pi		3.141592653589793
#define Euler_const	0.577215665
#define tolerant_err	1.0e-8
#define parallel_to_z_err		1.0e-6
/// 3*pi^2*theta, assume the same charge density as graphite, theta=0.382 A^-2, used for hypergeo nanotube
#define c3_pisq_theta 	11.3105666436484 
#define PascalA3_to_J_mol	6.0221415e-7 ///< Na*1e-30 turn preq*boxv to J/mol
#define J_mol_A3_to_Pascal	1660538.86313 ///< 1/6.0221415e-7 = 1.0/(6.0221415e23*1.0e-30)

int nspecie, nmole, natom; ///< total number of molecules, atoms, species 
int nmole_per_specie[nspecie_max]; ///< number of molecules in a certain specie
int natom_per_mole[nspecie_max]; ///< number of atoms in a molecule, which is belong to a certain specie
int nbond_per_mole[nspecie_max];
int nangle_per_mole[nspecie_max];
int ndih_per_mole[nspecie_max];
int nimp_per_mole[nspecie_max];
int nnbp_per_mole[nspecie_max];
int specie_first_atom_idx[nspecie_max+1];
int specie_first_mole_idx[nspecie_max+1];
int mole2specie[nmole_max]; ///< which specie this molecule belongs to
int mole_first_atom_idx[nmole_max+1]; ///< index of the first atom in a molecule
double mole_xx[nmole_max], mole_yy[nmole_max], mole_zz[nmole_max];
int vacancy_idx[nspecie_max]; ///< the index of first vancant molecule for a specie
int mole_status[nmole_max]; ///< the status of the molecule, e.g. vacancy etc.
int mole_first_bond_idx[nmole_max+1]; ///< index of the first bond in a molecule
int mole_first_angle_idx[nmole_max+1]; ///< index of the first angle in a molecule
int mole_first_dih_idx[nmole_max+1]; ///< index of the first dihedral in a molecule
int mole_first_imp_idx[nmole_max+1]; ///< index of the first improper in a molecule
int mole_first_nbp_idx[nmole_max+1]; ///< index of the first nonbonded pair in a molecule
double mw[nmole_max]; ///< molecule weight
int atom2mole[natom_max]; ///< which molecule this atom belongs to
double xx[natom_max], yy[natom_max], zz[natom_max]; ///< atom coordinates
/// inner coordinates relative to the center of mass, also used for PBC reconstruction of the molecule
double ex[natom_max], fy[natom_max], gz[natom_max];
double vx[natom_max], vy[natom_max], vz[natom_max]; ///< atom velocity
double fxl[natom_max], fyl[natom_max], fzl[natom_max]; ///< inter force
double fxs[natom_max], fys[natom_max], fzs[natom_max]; ///< intra force
double aw[natom_max], epsilon[natom_max], sigma[natom_max], charge[natom_max]; ///< atom weight etc. 
int isghost[natom_max]; ///< flag of ghost atom, including ghost LJ or ghost electrostatic
int tasostype[natom_max]; ///< atom type for tasos interpolation 
char atomname[natom_max][5];

/// for atom explicit solid (nanotubes)
int solid_natom;
double *solid_xx, *solid_yy, *solid_zz;
int fSolid_type; ///< 0-heterogeneous; 1-uniform (e.g. nanotubes); 
double *solid_sigma, *solid_epsilon, *solid_charge;

///< for hypergeometric nanotubes
int ntube;
double *hgntc_xx, *hgntc_yy, *hgnt_radius; ///< (h)yper(g)eometric (n)ano(t)ube (c)enter

int nbond, nangle, ndih, nimp, nnbp;
int bond_idx[nbond_max][2];
int bond_type[nbond_max];
double Kb[nbond_max], Req[nbond_max], alpha[nbond_max];
int angle_idx[nangle_max][3];
int angle_type[nangle_max];
double Ktheta[nangle_max], Thetaeq[nangle_max];
double agl_para_3[nangle_max], agl_para_4[nangle_max], agl_para_5[nangle_max];
int isAngle_unique[nangle_max]; ///< make sure 1,3 atoms do not make mutiple angles, this is for possible ring structures
int dih_idx[ndih_max][4];
int dih_type[ndih_max];
double c1[ndih_max], c2[ndih_max], c3[ndih_max], c4[ndih_max];
int isDih_unique[ndih_max]; ///< make sure 1,4 atoms do not make mutiple dihedrals, this is for possible ring structures
int imp_idx[nimp_max][4];
int imp_type[nimp_max];
double komega[nimp_max], omega0[nimp_max];
int nbp_idx[nnbp_max][2];
int excllist[exclude_max];
int pointexcl_specie[nspecie_max+1]; ///< the index of exclude list starts and ends for a specie
int pointexcl_atom[nspecie_max][natom_per_mole_max+1]; ///< the index of exclude list starts and ends for an atom

/// name of the object system
char sysname[200];

FILE *fpins, *fpouts, *fpcfg, *fplog;
FILE *fpss, *fptrj, *fpsave, *fpload;

int nstep, nstep_eq, nstep_start; ///< nstep_start is the starting step for continue runs 
int fStart_option; ///< 1-new run, 2-continue run, 3-new from old
/// steps for averages, print out, save, snapshots, trajectory
int nstep_ave, nstep_print, nstep_save, nstep_ss, nstep_trj;
int nstep_inner; ///< inner step for multi time step
double delt; ///< time step
double deltby2, delts, deltsby2;
int ij, jk;
double treq, preq; ///< input required temperature, pressure
double boxv, boxlx, boxly, boxlz; ///< box size
int nconstraint; ///< constraint, 3 for periodic, 6 for aperiodic
double rcutoff, rcutoffsq, rcuton, rcutonsq;
double rcutoffelec, rcutoffelecsq;
double f0; ///< 1,4 LJ potential modifier for OPLS, set to 1.0 for no modification or 0.5 for OPLS or 0.0 for TraPPE.
int isLJswitchOn; ///< use switch potential for LJ or not
/** 
 * \brief Electrostatic interaction type
 * 
 * 0 - No electrostatic interaction
 * 1 - Ewald Summation
 * 2 - Wolf's method
 * 3 - Simple direct coulomb interaction
 */
int iChargeType; 
int isEwaldOn; ///< Ewald summation electrostatic interactions
int fEwald_BC; ///< flag of boundary condition for ewald summation. 1-tinfoil; 2-vacuum
int fEwald_Dim; ///< flag of Ewald method dimension. 1, 2 or 3 dimension.
int isWolfOn; ///< wolf method for electrostatic interactions
int isSimpleCoulomb; ///< simple coulomb interactions on/off
double kappa, kappasq; ///< sqrt(alpha) in ewald summation. 
int KMAXX, KMAXY, KMAXZ; ///< ewald parameters
int KSQMAX; ///< ewald parameter
int what_simulation; ///< simulation type, MD, HMC, etc.
int what_ensemble; ///< what type of ensemble, NVT, NPT etc.
int whichNH; ///< which nose hoover subroutine to use? usually 3 for molecule, 2 for atoms, see more details in nvtnh.c
/** 
 * \brief solid-fluid type. for different nanotube potentials and future possible other materials
 * 
 * 0 - No solid-fluid interaction\n
 * 1 - Hypergeo nanotube potential\n
 * 2 - Atom explicit nanotube/other materials potentials
 * 3 - Tasos interpolation
 * 4 - Yang's interpolation
 */
int sf_type; 
char atomname[natom_max][5];

int istep; ///< counter of step, current step
double utot; ///< calculated in printit()
double upot; ///< calculated in printit()
double ukin;
double uinter, uintra; ///< inter and intra molecular energy
double uvdw; ///< van der wall energy, LJ energy, include unbp
double unbp_vdw; ///< nonbonded pair energy
double ubond, uangle, udih, uimp;
double uewald; ///< total ewald energy, refer to Frenkel and Smit, eq. 12.1.25
double ureal; ///< real part of ewald, term 3 in 12.1.25
double ufourier; ///< fourier part of ewald, term 1 in 12.1.25
double uself; ///< self interaction correction part of ewald, term 2 in 12.1.25
double uexcl; ///< excluding energy for ewald summation
double uvacuum; ///< vacuum boundary for ewald
double uGz0; ///< 1D ewald Gz=0 term
double LJswitch; ///< switch factor for LJ
double usflj; ///< solid-fluid LJ energy
double ucoulomb; ///< direct coulomb energy
double virial; ///< calculated in printit()
double virial_inter;
double virial_intra;
double pideal;
double pinst; ///< instantaneous pressure

/** 
 * The Session variables are used to store the results from one calculation session.
 * E.g., one loop_ij() session. They will be used to calculate the final total energies
 * or pair energies.
 */
double gUtotSession, gUpotSession, gUkinSession, gUinterSession, gUintraSession;
double gUvdwSession, gUvdwNbpSession, gUbondSession, gUangleSession,
		gUdihSession;
double gUimpSession, gUewaldSession, gUrealSession, gUfourierSession;
double gUselfSession, gUexclSession, gUvacuumSession, gUGz0Session,
		gUsfljSession;
double gUcoulombSession, gVirialSession, gVirialInterSession,
		gVirialIntraSession;
double gPidealSession, gPinstSession;
double gUwolfSession, gUwolfrealSession, gUwolfconSession;

int nfree; ///< freedom
double tinst; ///< instantaneous temperature

int nframe; ///< number of frames in the trajectory file

double TWOPI_LX, TWOPI_LY, TWOPI_LZ; ///< ewald
double Bfactor_ewald, Vfactor_ewald;
double twopi_over_3v; ///< constant for vacuum boundary

/// variables for wolf method
double uwolf, uwolf_real, uwolf_con, wolfvcon1, wolfvcon2, wolffcon1, wolffcon2;

double roff2_minus_ron2_cube; ///< used for switch potential

/// variables for nose hoover method
double Gts, Qts, vts, rts, dt_outer2, dt_outer4;
/// NPT 
double Gbs, Qbs, vbs, rbs, dt_outer8;
double utsbs; ///< extra energy for the barostat NPT

/// frenkel and smit's nose hoover method
double qq, ps, gg, ss;
double delt_sqby2, delts_sqby2;
double vxo[natom_max], vyo[natom_max], vzo[natom_max];
double vxn[natom_max], vyn[natom_max], vzn[natom_max];
double bx[natom_max], by[natom_max], bz[natom_max];
double unhts;
double qqs, pss, ggs, sss;
double unhtss;

// for HMC simulation
double prob_cm, ///< canonical move probability
		prob_vc, ///< volume change probability
		prob_id; ///< insertion/deletion probability
double ratio_cm_req, ///< required canonical move accept ratio
		ratio_vc_req; ///< required volume change accept ratio
int nstep_md_per_hmc; ///< steps of md moves per hmc cycle
int nstep_delt_adj_cycle; ///< steps between two delt adjustment
int nstep_delv_adj_cycle; ///< steps between two delv adjustment
double delv;

// my interpolations
// energy and forces interpolation parameters
double *ene0[nunique_atom_max], *ene1[nunique_atom_max],
		*ene2[nunique_atom_max], *ene3[nunique_atom_max],
		*ene4[nunique_atom_max], *ene5[nunique_atom_max],
		*ene6[nunique_atom_max], *ene7[nunique_atom_max],
		*ene8[nunique_atom_max], *ene9[nunique_atom_max],
		*ene10[nunique_atom_max], *ene11[nunique_atom_max],
		*ene12[nunique_atom_max], *ene13[nunique_atom_max],
		*ene14[nunique_atom_max], *ene15[nunique_atom_max],
		*ene16[nunique_atom_max], *ene17[nunique_atom_max],
		*ene18[nunique_atom_max], *ene19[nunique_atom_max],
		*ene20[nunique_atom_max], *ene21[nunique_atom_max],
		*ene22[nunique_atom_max], *ene23[nunique_atom_max],
		*ene24[nunique_atom_max], *ene25[nunique_atom_max],
		*ene26[nunique_atom_max], *ene27[nunique_atom_max],
		*ene28[nunique_atom_max], *ene29[nunique_atom_max],
		*ene30[nunique_atom_max], *ene31[nunique_atom_max];
double *fxa0[nunique_atom_max], *fxa1[nunique_atom_max],
		*fxa2[nunique_atom_max], *fxa3[nunique_atom_max],
		*fxa4[nunique_atom_max], *fxa5[nunique_atom_max],
		*fxa6[nunique_atom_max], *fxa7[nunique_atom_max],
		*fxa8[nunique_atom_max], *fxa9[nunique_atom_max],
		*fxa10[nunique_atom_max], *fxa11[nunique_atom_max],
		*fxa12[nunique_atom_max], *fxa13[nunique_atom_max],
		*fxa14[nunique_atom_max], *fxa15[nunique_atom_max],
		*fxa16[nunique_atom_max], *fxa17[nunique_atom_max],
		*fxa18[nunique_atom_max], *fxa19[nunique_atom_max],
		*fxa20[nunique_atom_max], *fxa21[nunique_atom_max],
		*fxa22[nunique_atom_max], *fxa23[nunique_atom_max],
		*fxa24[nunique_atom_max], *fxa25[nunique_atom_max],
		*fxa26[nunique_atom_max], *fxa27[nunique_atom_max],
		*fxa28[nunique_atom_max], *fxa29[nunique_atom_max],
		*fxa30[nunique_atom_max], *fxa31[nunique_atom_max];
double *fya0[nunique_atom_max], *fya1[nunique_atom_max],
		*fya2[nunique_atom_max], *fya3[nunique_atom_max],
		*fya4[nunique_atom_max], *fya5[nunique_atom_max],
		*fya6[nunique_atom_max], *fya7[nunique_atom_max],
		*fya8[nunique_atom_max], *fya9[nunique_atom_max],
		*fya10[nunique_atom_max], *fya11[nunique_atom_max],
		*fya12[nunique_atom_max], *fya13[nunique_atom_max],
		*fya14[nunique_atom_max], *fya15[nunique_atom_max],
		*fya16[nunique_atom_max], *fya17[nunique_atom_max],
		*fya18[nunique_atom_max], *fya19[nunique_atom_max],
		*fya20[nunique_atom_max], *fya21[nunique_atom_max],
		*fya22[nunique_atom_max], *fya23[nunique_atom_max],
		*fya24[nunique_atom_max], *fya25[nunique_atom_max],
		*fya26[nunique_atom_max], *fya27[nunique_atom_max],
		*fya28[nunique_atom_max], *fya29[nunique_atom_max],
		*fya30[nunique_atom_max], *fya31[nunique_atom_max];
double *fza0[nunique_atom_max], *fza1[nunique_atom_max],
		*fza2[nunique_atom_max], *fza3[nunique_atom_max],
		*fza4[nunique_atom_max], *fza5[nunique_atom_max],
		*fza6[nunique_atom_max], *fza7[nunique_atom_max],
		*fza8[nunique_atom_max], *fza9[nunique_atom_max],
		*fza10[nunique_atom_max], *fza11[nunique_atom_max],
		*fza12[nunique_atom_max], *fza13[nunique_atom_max],
		*fza14[nunique_atom_max], *fza15[nunique_atom_max],
		*fza16[nunique_atom_max], *fza17[nunique_atom_max],
		*fza18[nunique_atom_max], *fza19[nunique_atom_max],
		*fza20[nunique_atom_max], *fza21[nunique_atom_max],
		*fza22[nunique_atom_max], *fza23[nunique_atom_max],
		*fza24[nunique_atom_max], *fza25[nunique_atom_max],
		*fza26[nunique_atom_max], *fza27[nunique_atom_max],
		*fza28[nunique_atom_max], *fza29[nunique_atom_max],
		*fza30[nunique_atom_max], *fza31[nunique_atom_max];
double *interp_vector;
double uclx, ucly, uclz;
double grid_itvl_x, grid_itvl_y, grid_itvl_z;

int ngrid_x, ngrid_y, ngrid_z, ngrid_total;
int ncube_x, ncube_y, ncube_z, ncube_total;
double xcenter, ycenter, zcenter;
double xmin, xmax, ymin, ymax, zmin, zmax;

/// long range correction related variables
/// every unique molecule should have two types of long range corrections
/// one is intra-specie and the other is cross-specie(inter-specie)
double uljlrc, pljlrc;
double uljlrc_term[nspecie_max][nspecie_max];
double pljlrc_term[nspecie_max][nspecie_max];

// counters and accumulators
int icounter[num_counter_max];
double accumulator[num_counter_max][8];
// 0-4 for accumulator, 5-ave, 6-err, 7-fluc

/*! \var accumulator[num_counter_max][8]
 \brief Accumulators for various system states.
 
 accumulator[][0]	Raw accumulator for a variable within a block, e.g. utot	\n
 accumulator[][1]	Raw accumulator for the square of a variable within a block, e.g. utot*utot	\n
 accumulator[][2]	Block accumulator for variables average	\n
 accumulator[][3]	Block accumulator for averages for squared varaibles	\n
 accumulator[][4]	Block accumulator of (in block average)*(in block average)	\n
 accumulator[][5]	Final average	\n
 accumulator[][6]	Standard deviation	\n
 accumulator[][7]	Fluctuation	\n
 
 
 accumulator[0][]   utot	\n
 accumulator[1][]   upot	\n
 accumulator[2][]   ukin	\n
 accumulator[3][]   uinter	\n
 accumulator[4][]   uintra	\n
 accumulator[5][]   uvdw	\n
 accumulator[6][]   ubond	\n
 accumulator[7][]   uangle	\n
 accumulator[8][]   udih	\n
 accumulator[9][]   uimp	\n
 accumulator[10][]   uewald	\n
 accumulator[11][]   ureal	\n
 accumulator[12][]   ufourier	\n
 accumulator[13][]   uself	\n
 accumulator[14][]   usflj	\n
 accumulator[15][]   tinst	\n
 accumulator[16][]   uvacuum	\n
 accumulator[17][]   uwolf	\n
 accumulator[18][]   pressure	\n
 accumulator[19][]   boxv	\n
 accumulator[20][]   pideal	\n
 */

/*! \var icounter[num_counter_max]
 * \brief Counters.
 * 
 icounter[10]   number of average cycles	\n
 icounter[11]	decrease, for equilibrium	\n
 icounter[20]	number of canonical moves	\n
 icounter[21]   number of accepted canonical moves	\n
 icounter[22]   number of volume change moves	\n
 icounter[23]	number of accepted volume change moves	\n
 icounter[24]   number of insertion	\n
 icounter[25]   number of accepted insertion	\n
 icounter[26]   number of deletion	\n
 icounter[27]   number of accepted deletion	\n
 */

#endif /* VARS_H */
