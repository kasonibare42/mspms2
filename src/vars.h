#ifndef VARS_H
#define VARS_H   

/// File operators
FILE *fpins, *fpouts, *fpcfg, *fplog, *fpss, *fptrj, *fpsave, *fpload, *fpcoords;

/*************************************************************************************
 * Input parameters and their direct related variables
 *************************************************************************************/
char title[200]; ///< Title of the simulation
int ij, jk; ///< Random seeds
double sigma_base, epsilon_base, mass_base; ///< Base units for reduced units
double time_base, pressure_base; ///< Derived reduced units
double treq, preq; ///< input required temperature (K), pressure (bar)
double boxlx, boxly, boxlz; ///< box size
double boxv;
double rcutoff, rcuton, rcutoffelec;
double rcutoffsq, rcutonsq, rcutoffelecsq;
int isLJlrcOn; ///< if L-J long range correction is ON
int isLJswitchOn; ///< use switch potential for LJ or not
int iStart_option; ///< 1-new run, 2-continue run, 3-new from old
/// Number of steps for total, equilibrium, averages, print, save, snapshots, trajectory
int nstep, nstep_eq, nstep_ave, nstep_print, nstep_save, nstep_ss, nstep_trj;
int nstep_start, nstep_end; ///< nstep_start is the starting step for continue runs 
bool bEquilibrium; ///< Are we still in equilibrium run?
double delt; ///< time step
int nstep_inner; ///< inner step for multi time step
double deltby2, delts, deltsby2;
double f0; ///< 1,4 LJ potential modifier (1.0-no modification;0.5-OPLS;0.0-TraPPE)
int what_simulation; ///< simulation type, MD, HMC, etc.
int what_ensemble; ///< what type of ensemble, NVT, NPT etc.
int nconstraint; ///< constraint, 3 for periodic, 6 for aperiodic
int iInterMolePotType; ///< Inter molecular potential model type
int iChargeType; ///< Electrostatic interaction type
int iSF_type; ///< Solid-fluid interaction type
int iExternal_FF_type; ///< flag of other non-standard force fields
/*************************************************************************************
 * End of Input parameters and their direct related variables
 *************************************************************************************/


/*************************************************************************************
 * HMC Input parameters and their direct related variables
 *************************************************************************************/
double pdisp, ///< canonical move probability
		pvolm, ///< volume change probability
		pmake, ///< insertion probability
		pkill; ///< deletion probability
double delv;
double rreq_disp, ///< required canonical move accept ratio
		rreq_volm; ///< required volume change accept ratio
int nstep_md_per_hmc; ///< steps of md moves per hmc cycle
int nstep_delt_adj_cycle; ///< steps between two delt adjustment
int nstep_delv_adj_cycle; ///< steps between two delv adjustment
double pcomp[NSPECIE_MAX]; ///< The probability for a specie to be selected for insertion/deletion
double cpt[NSPECIE_MAX]; ///< pseudo chemcial potential
double zact[NSPECIE_MAX]; ///< activity


double pvolm_upper; ///< upper limit of volume change probability, i.e. pdisp+pvolm
double probability_to_insert[NSPECIE_MAX]; ///< The probability to insert a molecule for a specie
double fugacity_required[NSPECIE_MAX]; ///< required Fugacity of each specie 
/*************************************************************************************
 * End of HMC Input parameters and their direct related variables
 *************************************************************************************/

/** 
 * \brief Electrostatic interaction type
 * 
 * 0 - No electrostatic interaction\n
 * 1 - Ewald Summation\n
 * 2 - Wolf's method\n
 * 3 - Simple direct coulomb interaction\n
 */
int isEwaldOn; ///< Ewald summation electrostatic interactions
int fEwald_BC; ///< flag of boundary condition for ewald summation. 1-tinfoil; 2-vacuum
int fEwald_Dim; ///< flag of Ewald method dimension. 1, 2 or 3 dimension.
int isWolfOn; ///< wolf method for electrostatic interactions
int isSimpleCoulomb; ///< simple coulomb interactions on/off
double kappa, kappasq; ///< sqrt(alpha) in ewald summation. 
int KMAXX, KMAXY, KMAXZ; ///< ewald parameters
int KSQMAX; ///< ewald parameter
int whichNH; ///< which nose hoover subroutine to use? usually 3 for molecule, 2 for atoms, see more details in nvtnh.c
double TWOPI_LX, TWOPI_LY, TWOPI_LZ; ///< ewald
double Bfactor_ewald, Vfactor_ewald;
double twopi_over_3v; ///< constant for vacuum boundary


/// variables for nose hoover method
double Gts, Qts, vts, rts, dt_outer2, dt_outer4;
/// NPT 
double Gbs, Qbs, vbs, rbs, dt_outer8;
double utsbs; ///< extra energy for the barostat NPT

/// frenkel and smit's nose hoover method
double qq, ps, gg, ss;
double delt_sqby2, delts_sqby2;
double vxo[NATOM_MAX], vyo[NATOM_MAX], vzo[NATOM_MAX];
double vxn[NATOM_MAX], vyn[NATOM_MAX], vzn[NATOM_MAX];
double bx[NATOM_MAX], by[NATOM_MAX], bz[NATOM_MAX];
double unhts;
double qqs, pss, ggs, sss;
double unhtss;


/** 
 * \brief solid-fluid type. for different nanotube potentials and future possible other materials
 * 
 * 0 - No solid-fluid interaction\n
 * 1 - Hypergeo nanotube potential\n
 * 2 - Atom explicit nanotube/other materials potentials\n
 * 3 - Tasos interpolation\n
 * 4 - Yang's interpolation\n
 */

int istep; ///< counter of step, current step
double utot; ///< calculated in printit()
double upot; ///< calculated in printit()
double ukin;
double uinter; 
double uintra; ///< inter and intra molecular energy
double uvdw; ///< van der wall energy, LJ energy, include unbp
double usg; ///< Silvera-Goldman potential
double unbp_vdw; ///< nonbonded pair energy
double ubond; 
double uangle; 
double udih; 
double uimp;
double uelec; ///< Electrostatic interaction energy. This can be either uewald, uwolf, ucoulomb.
#define uewald 		uelec ///< total ewald energy, refer to Frenkel and Smit, eq. 12.1.25
#define uwolf 		uelec ///< total wolf energy.
#define ucoulomb 	uelec ///< direct coulomb energy
double ureal; ///< real part of ewald, term 3 in 12.1.25
double ufourier; ///< fourier part of ewald, term 1 in 12.1.25
double uself; ///< self interaction correction part of ewald, term 2 in 12.1.25
double uexcl; ///< excluding energy for ewald summation
double uvacuum; ///< vacuum boundary for ewald
double uGz0; ///< 1D ewald Gz=0 term
/// variables for wolf method
double uwolf_real; 
double uwolf_con; 
double wolfvcon1; 
double wolfvcon2; 
double wolfvcon3; 
double wolffcon1; 
double wolffcon2;
double LJswitch; ///< switch factor for LJ
double usflj; ///< solid-fluid LJ energy
double uotherff; ///< Other force field energy
#define udftmcff		uotherff ///< DFT metal cluster FF energy
double virial; ///< calculated in printit()
double virial_inter;
double virial_intra;
double pideal;
double pinst; ///< instantaneous pressure

double ushift; ///< shift energy for cut and shift

double coulomb_prefactor; ///< This is the reduced value of COULOMB_PREFACTOR

/// the constant terms for shift energies
/// They are related to rcutoff
/// shift1 = (1/rcutoff)^6
/// shift4 = 4*shift;
double shift4, shift1; 
double sgshift; ///< shift energy for SG potential

int ncut;

/**
 * Variables used to save the old states
 */
double *xx_old, *yy_old, *zz_old;
double upot_old, utot_old;
double upot_new, utot_new;
double uljlrc_old, pljlrc_old;
double ukin_old, tinst_old, uinter_old, uintra_old, uvdw_old, ubond_old,
		uangle_old, udih_old, uimp_old, uewald_old, uwolf_old, ucoulomb_old, usflj_old, 
		unhts_old, unhtss_old, virial_inter_old, virial_intra_old,
		utsbs_old, pinst_old, boxlx_old, boxly_old, boxlz_old, boxv_old;
double ushift_old;

/** 
 * The Session variables are used to store the results from one calculation session.
 * E.g., one loop_ij() session. They will be used to calculate the final total energies
 * or pair energies.
 */
double gUtotSession, gUpotSession, gUkinSession, gUinterSession, gUintraSession;
double gUvdwSession, gUvdwNbpSession, gUbondSession, gUangleSession,
		gUdihSession;
double gUsgSession; ///< Silvera-goldman potential
double gUimpSession, gUewaldSession, gUrealSession, gUfourierSession;
double gUselfSession, gUexclSession, gUvacuumSession, gUGz0Session,
		gUsfljSession;
double gUcoulombSession, gVirialSession, gVirialInterSession,
		gVirialIntraSession;
double gPidealSession, gPinstSession;
double gUwolfSession, gUwolfrealSession, gUwolfconSession;
double gUMetalClusterSession;
double gUShiftSession;
int gNcut;

int nfree; ///< freedom
double tinst; ///< instantaneous temperature

int nframe; ///< number of frames in the trajectory file

/// Variables for Silvera-Goldman potential
double sg_alpha, sg_beta, sg_gama, sg_c6,
       sg_c8, sg_c9, sg_c10, sg_rc;
double spring[NSPECIE_MAX]; ///< spring constant for polymer bead ring

double roff2_minus_ron2_cube; ///< used for switch potential


/* set up parameters for this simulated annealing run */
/* how many points do we try before stepping (move to the next temperature) */
double n_tries;
/* how many iterations for each T? (how many MD steps before use the acceptance critirea) */
int iters_fixed_t;
/* initial temperature */
double t_initial;
/* damping factor for temperature */
double mu_t, t_min;

typedef struct _dftmcffparam
{
	int iWhichAxis;
	int index;
} DFTMCFFPARAM;

int isAccept;
double rndnum[3];
double dH;
double fDeltaU;
double fDeltaUljlrc, fDeltaPljlrc;
double fVolumeNew;
double fLengthNew, fWidthNew, fHeightNew;
double fRatioNewV2OldV;
double fRatioNewL2OldL;
double fMinHalf;

/// for atom explicit solid (nanotubes)
int solid_natom;
double *solid_xx, *solid_yy, *solid_zz;
int fSolid_type; ///< 0-heterogeneous; 1-uniform (e.g. nanotubes); 
double *solid_sigma, *solid_epsilon, *solid_charge;

///< for hypergeometric nanotubes
int ntube;
double *hgntc_xx, *hgntc_yy, *hgnt_radius; ///< (h)yper(g)eometric (n)ano(t)ube (c)enter

// my interpolations
// energy and forces interpolation parameters
double *ene0[NUNIQUE_ATOM_MAX], *ene1[NUNIQUE_ATOM_MAX],
		*ene2[NUNIQUE_ATOM_MAX], *ene3[NUNIQUE_ATOM_MAX],
		*ene4[NUNIQUE_ATOM_MAX], *ene5[NUNIQUE_ATOM_MAX],
		*ene6[NUNIQUE_ATOM_MAX], *ene7[NUNIQUE_ATOM_MAX],
		*ene8[NUNIQUE_ATOM_MAX], *ene9[NUNIQUE_ATOM_MAX],
		*ene10[NUNIQUE_ATOM_MAX], *ene11[NUNIQUE_ATOM_MAX],
		*ene12[NUNIQUE_ATOM_MAX], *ene13[NUNIQUE_ATOM_MAX],
		*ene14[NUNIQUE_ATOM_MAX], *ene15[NUNIQUE_ATOM_MAX],
		*ene16[NUNIQUE_ATOM_MAX], *ene17[NUNIQUE_ATOM_MAX],
		*ene18[NUNIQUE_ATOM_MAX], *ene19[NUNIQUE_ATOM_MAX],
		*ene20[NUNIQUE_ATOM_MAX], *ene21[NUNIQUE_ATOM_MAX],
		*ene22[NUNIQUE_ATOM_MAX], *ene23[NUNIQUE_ATOM_MAX],
		*ene24[NUNIQUE_ATOM_MAX], *ene25[NUNIQUE_ATOM_MAX],
		*ene26[NUNIQUE_ATOM_MAX], *ene27[NUNIQUE_ATOM_MAX],
		*ene28[NUNIQUE_ATOM_MAX], *ene29[NUNIQUE_ATOM_MAX],
		*ene30[NUNIQUE_ATOM_MAX], *ene31[NUNIQUE_ATOM_MAX];
double *fxa0[NUNIQUE_ATOM_MAX], *fxa1[NUNIQUE_ATOM_MAX],
		*fxa2[NUNIQUE_ATOM_MAX], *fxa3[NUNIQUE_ATOM_MAX],
		*fxa4[NUNIQUE_ATOM_MAX], *fxa5[NUNIQUE_ATOM_MAX],
		*fxa6[NUNIQUE_ATOM_MAX], *fxa7[NUNIQUE_ATOM_MAX],
		*fxa8[NUNIQUE_ATOM_MAX], *fxa9[NUNIQUE_ATOM_MAX],
		*fxa10[NUNIQUE_ATOM_MAX], *fxa11[NUNIQUE_ATOM_MAX],
		*fxa12[NUNIQUE_ATOM_MAX], *fxa13[NUNIQUE_ATOM_MAX],
		*fxa14[NUNIQUE_ATOM_MAX], *fxa15[NUNIQUE_ATOM_MAX],
		*fxa16[NUNIQUE_ATOM_MAX], *fxa17[NUNIQUE_ATOM_MAX],
		*fxa18[NUNIQUE_ATOM_MAX], *fxa19[NUNIQUE_ATOM_MAX],
		*fxa20[NUNIQUE_ATOM_MAX], *fxa21[NUNIQUE_ATOM_MAX],
		*fxa22[NUNIQUE_ATOM_MAX], *fxa23[NUNIQUE_ATOM_MAX],
		*fxa24[NUNIQUE_ATOM_MAX], *fxa25[NUNIQUE_ATOM_MAX],
		*fxa26[NUNIQUE_ATOM_MAX], *fxa27[NUNIQUE_ATOM_MAX],
		*fxa28[NUNIQUE_ATOM_MAX], *fxa29[NUNIQUE_ATOM_MAX],
		*fxa30[NUNIQUE_ATOM_MAX], *fxa31[NUNIQUE_ATOM_MAX];
double *fya0[NUNIQUE_ATOM_MAX], *fya1[NUNIQUE_ATOM_MAX],
		*fya2[NUNIQUE_ATOM_MAX], *fya3[NUNIQUE_ATOM_MAX],
		*fya4[NUNIQUE_ATOM_MAX], *fya5[NUNIQUE_ATOM_MAX],
		*fya6[NUNIQUE_ATOM_MAX], *fya7[NUNIQUE_ATOM_MAX],
		*fya8[NUNIQUE_ATOM_MAX], *fya9[NUNIQUE_ATOM_MAX],
		*fya10[NUNIQUE_ATOM_MAX], *fya11[NUNIQUE_ATOM_MAX],
		*fya12[NUNIQUE_ATOM_MAX], *fya13[NUNIQUE_ATOM_MAX],
		*fya14[NUNIQUE_ATOM_MAX], *fya15[NUNIQUE_ATOM_MAX],
		*fya16[NUNIQUE_ATOM_MAX], *fya17[NUNIQUE_ATOM_MAX],
		*fya18[NUNIQUE_ATOM_MAX], *fya19[NUNIQUE_ATOM_MAX],
		*fya20[NUNIQUE_ATOM_MAX], *fya21[NUNIQUE_ATOM_MAX],
		*fya22[NUNIQUE_ATOM_MAX], *fya23[NUNIQUE_ATOM_MAX],
		*fya24[NUNIQUE_ATOM_MAX], *fya25[NUNIQUE_ATOM_MAX],
		*fya26[NUNIQUE_ATOM_MAX], *fya27[NUNIQUE_ATOM_MAX],
		*fya28[NUNIQUE_ATOM_MAX], *fya29[NUNIQUE_ATOM_MAX],
		*fya30[NUNIQUE_ATOM_MAX], *fya31[NUNIQUE_ATOM_MAX];
double *fza0[NUNIQUE_ATOM_MAX], *fza1[NUNIQUE_ATOM_MAX],
		*fza2[NUNIQUE_ATOM_MAX], *fza3[NUNIQUE_ATOM_MAX],
		*fza4[NUNIQUE_ATOM_MAX], *fza5[NUNIQUE_ATOM_MAX],
		*fza6[NUNIQUE_ATOM_MAX], *fza7[NUNIQUE_ATOM_MAX],
		*fza8[NUNIQUE_ATOM_MAX], *fza9[NUNIQUE_ATOM_MAX],
		*fza10[NUNIQUE_ATOM_MAX], *fza11[NUNIQUE_ATOM_MAX],
		*fza12[NUNIQUE_ATOM_MAX], *fza13[NUNIQUE_ATOM_MAX],
		*fza14[NUNIQUE_ATOM_MAX], *fza15[NUNIQUE_ATOM_MAX],
		*fza16[NUNIQUE_ATOM_MAX], *fza17[NUNIQUE_ATOM_MAX],
		*fza18[NUNIQUE_ATOM_MAX], *fza19[NUNIQUE_ATOM_MAX],
		*fza20[NUNIQUE_ATOM_MAX], *fza21[NUNIQUE_ATOM_MAX],
		*fza22[NUNIQUE_ATOM_MAX], *fza23[NUNIQUE_ATOM_MAX],
		*fza24[NUNIQUE_ATOM_MAX], *fza25[NUNIQUE_ATOM_MAX],
		*fza26[NUNIQUE_ATOM_MAX], *fza27[NUNIQUE_ATOM_MAX],
		*fza28[NUNIQUE_ATOM_MAX], *fza29[NUNIQUE_ATOM_MAX],
		*fza30[NUNIQUE_ATOM_MAX], *fza31[NUNIQUE_ATOM_MAX];
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
double uljlrc_term[NSPECIE_MAX][NSPECIE_MAX];
double pljlrc_term[NSPECIE_MAX][NSPECIE_MAX];

// counters and accumulators
int counts[NCOUNTS_MAX];
double accum[NCOUNTS_MAX][8];
// 0-4 for accumulator, 5-ave, 6-err, 7-fluc

/*! \var accum[NCOUNTS_MAX][8]
 \brief Accumulators for various system states.
 
 accum[][0]	Raw accumulator for a variable within a block, e.g. utot	\n
 accum[][1]	Raw accumulator for the square of a variable within a block, e.g. utot*utot	\n
 accum[][2]	Block accumulator for variables average	\n
 accum[][3]	Block accumulator for averages for squared varaibles	\n
 accum[][4]	Block accumulator of (in block average)*(in block average)	\n
 accum[][5]	Final average	\n
 accum[][6]	Standard deviation	\n
 accum[][7]	Fluctuation	\n
 
 
 accum[0][]   utot	\n
 accum[1][]   upot	\n
 accum[2][]   ukin	\n
 accum[3][]   uinter	\n
 accum[4][]   uintra	\n
 accum[5][]   uvdw	\n
 accum[6][]   ubond	\n
 accum[7][]   uangle	\n
 accum[8][]   udih	\n
 accum[9][]   uimp	\n
 accum[10][]   uewald	\n
 accum[11][]   ureal	\n
 accum[12][]   ufourier	\n
 accum[13][]   uself	\n
 accum[14][]   usflj	\n
 accum[15][]   tinst	\n
 accum[16][]   uvacuum	\n
 accum[17][]   uwolf	\n
 accum[18][]   pressure	\n
 accum[19][]   boxv	\n
 accum[20][]   pideal	\n
 */

/*! \var counts[NCOUNTS_MAX]
 * \brief Counters.
 * 
 counts[10]   number of average cycles	\n
 counts[11]	decrease, for equilibrium	\n
 
 Following counters will be set to zero periodically during equilibrium run.
 After equilibrium, they will be used to calculate the final acceptance ratio.\n
 counts[20]	number of canonical moves	\n
 counts[21]   number of accepted canonical moves	\n
 
 counts[23]   number of volume change moves	\n
 counts[24]	number of accepted volume change moves	\n
 
 counts[26]   number of insertion	\n
 counts[27]   number of accepted insertion	\n
 counts[28]   number of deletion	\n
 counts[29]   number of accepted deletion	\n
 */

#endif /* VARS_H */
