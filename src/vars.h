#ifndef VARS_H
#define VARS_H   

/// File operators
FILE *fpins; 
FILE *fpouts; 
FILE *fpcfg; 
FILE *fplog; 
FILE *fpss; 
FILE *fptrj; 
FILE *fpsave; 
FILE *fpload;
FILE *fpcoords;

/*************************************************************************************
 * Input parameters and their direct related variables
 *************************************************************************************/
char title[LONG_STRING_LENGTH]; ///< Title of the simulation
int ij; 
int jk; ///< Random seeds
double sigma_base;
double epsilon_base; 
double mass_base; ///< Base units for reduced units
double time_base;
double pressure_base; ///< Derived reduced units
double force_base;
double velocity_base;
double treq;
double preq; ///< input required temperature (K), pressure (bar)
double boxlx; 
double boxly; 
double boxlz; ///< box size
double boxv;
double rcutoff; 
double rcuton;
double rcutoffelec;
double rcutoffsq;
double rcutonsq;
double rcutoffelecsq;
int isLJlrcOn; ///< if L-J long range correction is ON
int isLJswitchOn; ///< use switch potential for LJ or not
int iStart_option; ///< 1-new run, 2-continue run, 3-new from old
/// Number of steps for total, equilibrium, averages, print, save, snapshots, trajectory
int nstep;
int nstep_eq;
int nstep_ave;
int nstep_print;
int nstep_save;
int nstep_ss;
int nstep_trj;
int nstep_start;
int nstep_end; ///< nstep_start is the starting step for continue runs 
bool bEquilibrium; ///< Are we still in equilibrium run?
double delt; ///< time step
int nstep_inner; ///< inner step for multi time step
double deltby2;
double delts;
double deltsby2;
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
double pdisp; ///< canonical move probability
double pvolm; ///< volume change probability
double pmake; ///< insertion probability
double pkill; ///< deletion probability
double delv;
double rreq_disp; ///< required canonical move accept ratio
double rreq_volm; ///< required volume change accept ratio
double rinst_disp; ///< Instantaneous ratio for displacement
double rinst_volm; ///< Instantaneous ratio for volume change
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
int fEwald_BC; ///< flag of boundary condition for ewald summation. 1-tinfoil; 2-vacuum
int fEwald_Dim; ///< flag of Ewald method dimension. 1, 2 or 3 dimension.
double kappa; 
double kappasq; ///< sqrt(alpha) in ewald summation. 
int KMAXX; 
int KMAXY; 
int KMAXZ; ///< ewald parameters
int KSQMAX; ///< ewald parameter
int whichNH; ///< which nose hoover subroutine to use? usually 3 for molecule, 2 for atoms, see more details in nvtnh.c
double TWOPI_LX;
double TWOPI_LY;
double TWOPI_LZ; ///< ewald
double Bfactor_ewald; 
double Vfactor_ewald;
double twopi_over_3v; ///< constant for vacuum boundary

/// variables for nose hoover method
double Gts; 
double Qts;
double vts;
double rts;
double dt_outer2;
double dt_outer4;
/// NPT 
double Gbs; 
double Qbs;
double vbs;
double rbs;
double dt_outer8;
double utsbs; ///< extra energy for the barostat NPT

/// frenkel and smit's nose hoover method
double qq;
double ps;
double gg;
double ss;
double delt_sqby2;
double delts_sqby2;
double vxo[NATOM_MAX];
double vyo[NATOM_MAX];
double vzo[NATOM_MAX];
double vxn[NATOM_MAX];
double vyn[NATOM_MAX];
double vzn[NATOM_MAX];
double bx[NATOM_MAX];
double by[NATOM_MAX];
double bz[NATOM_MAX];
double unhts;
double qqs;
double pss;
double ggs;
double sss;
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

/**
 * Variables used to save the old states
 */
double *xx_old; 
double *yy_old;
double *zz_old;
double upot_old;
double utot_old;
double upot_new;
double utot_new;
double uljlrc_old;
double pljlrc_old;
double ukin_old;
double tinst_old;
double uinter_old;
double uintra_old;
double uvdw_old;
double ubond_old;
double uangle_old;
double udih_old;
double uimp_old;
double uewald_old; 
double uwolf_old;
double ucoulomb_old;
double usflj_old;
double unhts_old;
double unhtss_old;
double virial_inter_old;
double virial_intra_old;
double utsbs_old;
double pinst_old;
double boxlx_old;
double boxly_old;
double boxlz_old;
double boxv_old;
double ushift_old;

/** 
 * The Session variables are used to store the results from one calculation session.
 * E.g., one loop_ij() session. They will be used to calculate the final total energies
 * or pair energies.
 */

int nfree; ///< freedom
double tinst; ///< instantaneous temperature

int nframe; ///< number of frames in the trajectory file

double roff2_minus_ron2_cube; ///< used for switch potential


/* set up parameters for this simulated annealing run */
/* how many points do we try before stepping (move to the next temperature) */
double n_tries;
/* how many iterations for each T? (how many MD steps before use the acceptance critirea) */
int iters_fixed_t;
/* initial temperature */
double t_initial;
/* damping factor for temperature */
double mu_t;
double t_min;

typedef struct _dftmcffparam
{
	int iWhichAxis;
	int index;
} DFTMCFFPARAM;

bool isAccept;
double rndnum[3];
double dH;
double fDeltaU;
double fDeltaUljlrc;
double fDeltaPljlrc;
double fVolumeNew;
double fLengthNew;
double fWidthNew;
double fHeightNew;
double fRatioNewV2OldV;
double fRatioNewL2OldL;
double fMinHalf;

/// for atom explicit solid (nanotubes)
int solid_natom;
double *solid_xx;
double *solid_yy;
double *solid_zz;
int fSolid_type; ///< 0-heterogeneous; 1-uniform (e.g. nanotubes); 
double *solid_sigma;
double *solid_epsilon;
double *solid_charge;

///< for hypergeometric nanotubes
int ntube;
double *hgntc_xx;
double *hgntc_yy;
double *hgnt_radius; ///< (h)yper(g)eometric (n)ano(t)ube (c)enter
double const_3pisq_theta;

// my interpolations
// energy and forces interpolation parameters
double *ene0[NUNIQUE_ATOM_MAX];
double *ene1[NUNIQUE_ATOM_MAX];
double *ene2[NUNIQUE_ATOM_MAX];
double *ene3[NUNIQUE_ATOM_MAX];
double *ene4[NUNIQUE_ATOM_MAX];
double *ene5[NUNIQUE_ATOM_MAX];
double *ene6[NUNIQUE_ATOM_MAX];
double *ene7[NUNIQUE_ATOM_MAX];
double *ene8[NUNIQUE_ATOM_MAX];
double *ene9[NUNIQUE_ATOM_MAX];
double *ene10[NUNIQUE_ATOM_MAX];
double *ene11[NUNIQUE_ATOM_MAX];
double *ene12[NUNIQUE_ATOM_MAX];
double *ene13[NUNIQUE_ATOM_MAX];
double *ene14[NUNIQUE_ATOM_MAX];
double *ene15[NUNIQUE_ATOM_MAX];
double *ene16[NUNIQUE_ATOM_MAX];
double *ene17[NUNIQUE_ATOM_MAX];
double *ene18[NUNIQUE_ATOM_MAX];
double *ene19[NUNIQUE_ATOM_MAX];
double *ene20[NUNIQUE_ATOM_MAX];
double *ene21[NUNIQUE_ATOM_MAX];
double *ene22[NUNIQUE_ATOM_MAX];
double *ene23[NUNIQUE_ATOM_MAX];
double *ene24[NUNIQUE_ATOM_MAX];
double *ene25[NUNIQUE_ATOM_MAX];
double *ene26[NUNIQUE_ATOM_MAX];
double *ene27[NUNIQUE_ATOM_MAX];
double *ene28[NUNIQUE_ATOM_MAX];
double *ene29[NUNIQUE_ATOM_MAX]; 
double *ene30[NUNIQUE_ATOM_MAX];
double *ene31[NUNIQUE_ATOM_MAX];
double *fxa0[NUNIQUE_ATOM_MAX];
double *fxa1[NUNIQUE_ATOM_MAX];
double *fxa2[NUNIQUE_ATOM_MAX];
double *fxa3[NUNIQUE_ATOM_MAX];
double *fxa4[NUNIQUE_ATOM_MAX];
double *fxa5[NUNIQUE_ATOM_MAX];
double *fxa6[NUNIQUE_ATOM_MAX];
double *fxa7[NUNIQUE_ATOM_MAX];
double *fxa8[NUNIQUE_ATOM_MAX];
double *fxa9[NUNIQUE_ATOM_MAX];
double *fxa10[NUNIQUE_ATOM_MAX];
double *fxa11[NUNIQUE_ATOM_MAX];
double *fxa12[NUNIQUE_ATOM_MAX];
double *fxa13[NUNIQUE_ATOM_MAX];
double *fxa14[NUNIQUE_ATOM_MAX];
double *fxa15[NUNIQUE_ATOM_MAX];
double *fxa16[NUNIQUE_ATOM_MAX];
double *fxa17[NUNIQUE_ATOM_MAX];
double *fxa18[NUNIQUE_ATOM_MAX];
double *fxa19[NUNIQUE_ATOM_MAX];
double *fxa20[NUNIQUE_ATOM_MAX];
double *fxa21[NUNIQUE_ATOM_MAX];
double *fxa22[NUNIQUE_ATOM_MAX];
double *fxa23[NUNIQUE_ATOM_MAX];
double *fxa24[NUNIQUE_ATOM_MAX];
double *fxa25[NUNIQUE_ATOM_MAX];
double *fxa26[NUNIQUE_ATOM_MAX];
double *fxa27[NUNIQUE_ATOM_MAX];
double *fxa28[NUNIQUE_ATOM_MAX];
double *fxa29[NUNIQUE_ATOM_MAX];
double *fxa30[NUNIQUE_ATOM_MAX];
double *fxa31[NUNIQUE_ATOM_MAX];
double *fya0[NUNIQUE_ATOM_MAX];
double *fya1[NUNIQUE_ATOM_MAX];
double *fya2[NUNIQUE_ATOM_MAX];
double *fya3[NUNIQUE_ATOM_MAX];
double *fya4[NUNIQUE_ATOM_MAX];
double *fya5[NUNIQUE_ATOM_MAX];
double *fya6[NUNIQUE_ATOM_MAX];
double *fya7[NUNIQUE_ATOM_MAX];
double *fya8[NUNIQUE_ATOM_MAX];
double *fya9[NUNIQUE_ATOM_MAX];
double *fya10[NUNIQUE_ATOM_MAX];
double *fya11[NUNIQUE_ATOM_MAX];
double *fya12[NUNIQUE_ATOM_MAX];
double *fya13[NUNIQUE_ATOM_MAX];
double *fya14[NUNIQUE_ATOM_MAX];
double *fya15[NUNIQUE_ATOM_MAX];
double *fya16[NUNIQUE_ATOM_MAX];
double *fya17[NUNIQUE_ATOM_MAX];
double *fya18[NUNIQUE_ATOM_MAX];
double *fya19[NUNIQUE_ATOM_MAX];
double *fya20[NUNIQUE_ATOM_MAX];
double *fya21[NUNIQUE_ATOM_MAX];
double *fya22[NUNIQUE_ATOM_MAX];
double *fya23[NUNIQUE_ATOM_MAX];
double *fya24[NUNIQUE_ATOM_MAX];
double *fya25[NUNIQUE_ATOM_MAX];
double *fya26[NUNIQUE_ATOM_MAX];
double *fya27[NUNIQUE_ATOM_MAX];
double *fya28[NUNIQUE_ATOM_MAX];
double *fya29[NUNIQUE_ATOM_MAX];
double *fya30[NUNIQUE_ATOM_MAX];
double *fya31[NUNIQUE_ATOM_MAX];
double *fza0[NUNIQUE_ATOM_MAX];
double *fza1[NUNIQUE_ATOM_MAX];
double *fza2[NUNIQUE_ATOM_MAX];
double *fza3[NUNIQUE_ATOM_MAX];
double *fza4[NUNIQUE_ATOM_MAX];
double *fza5[NUNIQUE_ATOM_MAX];
double *fza6[NUNIQUE_ATOM_MAX];
double *fza7[NUNIQUE_ATOM_MAX];
double *fza8[NUNIQUE_ATOM_MAX];
double *fza9[NUNIQUE_ATOM_MAX];
double *fza10[NUNIQUE_ATOM_MAX];
double *fza11[NUNIQUE_ATOM_MAX];
double *fza12[NUNIQUE_ATOM_MAX];
double *fza13[NUNIQUE_ATOM_MAX];
double *fza14[NUNIQUE_ATOM_MAX];
double *fza15[NUNIQUE_ATOM_MAX];
double *fza16[NUNIQUE_ATOM_MAX];
double *fza17[NUNIQUE_ATOM_MAX];
double *fza18[NUNIQUE_ATOM_MAX];
double *fza19[NUNIQUE_ATOM_MAX];
double *fza20[NUNIQUE_ATOM_MAX];
double *fza21[NUNIQUE_ATOM_MAX];
double *fza22[NUNIQUE_ATOM_MAX];
double *fza23[NUNIQUE_ATOM_MAX];
double *fza24[NUNIQUE_ATOM_MAX];
double *fza25[NUNIQUE_ATOM_MAX];
double *fza26[NUNIQUE_ATOM_MAX];
double *fza27[NUNIQUE_ATOM_MAX];
double *fza28[NUNIQUE_ATOM_MAX];
double *fza29[NUNIQUE_ATOM_MAX];
double *fza30[NUNIQUE_ATOM_MAX];
double *fza31[NUNIQUE_ATOM_MAX];

double *interp_vector;
double uclx;
double ucly;
double uclz;
double grid_itvl_x;
double grid_itvl_y;
double grid_itvl_z;

int ngrid_x;
int ngrid_y;
int ngrid_z;
int ngrid_total;
int ncube_x;
int ncube_y;
int ncube_z;
int ncube_total;
double xcenter;
double ycenter;
double zcenter;
double xmin;
double xmax;
double ymin;
double ymax;
double zmin;
double zmax;

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
