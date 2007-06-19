#ifndef VARS_H
#define VARS_H   

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

#define not_ghost	0
#define lj_ghost	1
#define all_ghost	2

#define num_counter_max	30

#define new_run		1
#define continue_run	2
#define new_from_old	3

#define bond_none		0
#define bond_harmonic		1
#define bond_morse		2
#define bond_fene		3

#define angle_none		0
#define angle_harmonic		1
#define angle_TRwater		2  // Toukan and Rahman water potentials

#define dih_none		0
#define dih_opls_cosin		1
#define dih_charmm		2 // charmm type dihedral potential

#define imp_none		0
#define imp_charmm		1 // charmm type improper potential

#define nanotube_hypergeo  	1
#define nanotube_atom_explicit	2
#define nanotube_tasos		3
#define nanotube_my_interp	4 // my interpolation grid

#define ewald_bc_tinfoil	1
#define ewald_bc_vaccum		2

#define ewald_1D		1
#define ewald_2D		2
#define ewald_3D		3

#define solid_hetero		0  // solid type, e.g. mof?
#define solid_uniform		1 // e.g. nanotoubes

#define nspecie_max	3
#define nmole_max	1000
#define natom_max	4000
#define nbond_max	2000
#define nangle_max	2000
#define ndih_max	2000
#define nimp_max	2000
#define nnbp_max	2000
#define exclude_max	4000

#define nunique_atom_max	5 // number of unqiue atoms for grid interpolation

#define Rgas		8.314472 /* J/mol/K */
#define rRgas		0.120272219 /* reciprocal of Rgas */
#define const_columb	1389355.1051  // unit is J/mol. Na*(1.602177e-19)^2/4/PI/epsilon0/(1.0e-10)
#define pi		3.141592653589793
#define Euler_const	0.577215665
#define parallel_to_z_err		1.0e-6
#define c3_pisq_theta 	11.3105666436484 // 3*pi^2*theta, assume the same charge density as graphite, theta=0.382 A^-2


int nspecie, nmole, natom; /* total number of molecules, atoms, species */
int mole_first_atom_idx[nmole_max+1]; // index of the first atom in a molecule
double mole_xx[nmole_max], mole_yy[nmole_max], mole_zz[nmole_max];
int mole_status[nmole_max]; // the status of the molecule, e.g. vacancy etc.
int mole_first_bond_idx[nmole_max+1]; // index of the first bond in a molecule
int mole_first_angle_idx[nmole_max+1]; // index of the first angle in a molecule
int mole_first_dih_idx[nmole_max+1]; // index of the first dihedral in a molecule
int mole_first_imp_idx[nmole_max+1]; // index of the first improper in a molecule
double mw[nmole_max]; // molecule weight
double xx[natom_max], yy[natom_max], zz[natom_max]; /* position */
// inner coordinates relative to the center of mass, also used for PBC reconstruction of the molecule
double ex[natom_max], fy[natom_max], gz[natom_max]; 
double vx[natom_max], vy[natom_max], vz[natom_max]; /* velocity */
double fxl[natom_max], fyl[natom_max], fzl[natom_max]; /* inter force */
double fxs[natom_max], fys[natom_max], fzs[natom_max]; /* intra force */
double aw[natom_max], epsilon[natom_max], sigma[natom_max], charge[natom_max]; /* atom weight etc. */
int isghost[natom_max]; /* flag of ghost atom, including ghost LJ or ghost electrostatic */
int tasostype[natom_max]; /* atom type for tasos interpolation */
char atomname[natom_max][5];
int atom_to_mole_idx[natom_max]; // the index of the molecule which the atom is belong to

// for atom explicit solid (nanotubes)
int solid_natom;
double *solid_xx, *solid_yy, *solid_zz;
int fSolid_type; // 0-heterogeneous; 1-uniform (e.g. nanotubes); 
double *solid_sigma, *solid_epsilon, *solid_charge;

// for hypergeometric nanotubes
int ntube;
double *hgntc_xx, *hgntc_yy, *hgnt_radius; // (h)yper(g)eometric (n)ano(t)ube (c)enter

int nbond, nangle, ndih, nimp, nnbp;
int bond_idx[nbond_max][2];
int bond_type[nbond_max]; 
double Kb[nbond_max], Req[nbond_max], alpha[nbond_max];
int angle_idx[nangle_max][3];
int angle_type[nangle_max];
double Ktheta[nangle_max], Thetaeq[nangle_max];
double agl_para_3[nangle_max], agl_para_4[nangle_max], agl_para_5[nangle_max];
int isAngle_unique[nangle_max]; // make 1,3 atoms do not make mutiple angles, this is for possible ring structures
int dih_idx[ndih_max][4];
int dih_type[ndih_max];
double c1[ndih_max], c2[ndih_max], c3[ndih_max], c4[ndih_max];
int isDih_unique[ndih_max]; // make sure 1,4 atoms do not make mutiple dihedrals, this is for possible ring structures
int imp_idx[nimp_max][4];
int imp_type[nimp_max];
double komega[nimp_max], omega0[nimp_max];
int nbp_idx[nnbp_max][2];
int excllist[exclude_max];
int pointexcl[natom_max+1]; // the index of exclude list starts and ends

char sysname[200];

FILE *fpins, *fpouts, *fpcfg, *fplog;
FILE *fpss, *fptrj, *fpsave, *fpload;

int nstep, nstep_eq, nstep_start; /* nstep_start is the starting step for continue runs */
int fStart_option; // 1-new run, 2-continue run, 3-new from old
// steps for averages, print out, save, snapshots, trajectory
int nstep_ave, nstep_print, nstep_save, nstep_ss, nstep_trj;
int nstep_inner; /* inner step for multi time step */
double delt; /* time step */
double deltby2, delts, deltsby2;
int ij, jk;
double treq, preq; /* input required temperature, pressure */
double boxlx, boxly, boxlz; /* box size */
int nconstraint; // constraint, 3 for periodic, 6 for aperiodic
double rcutoff, rcutoffsq, rcuton, rcutonsq;
double rcutoffelec, rcutoffelecsq;
double f0; // 1,4 LJ potential modifier for OPLS, set to 1.0 for no modification or 0.5 for OPLS or 0.0 for TraPPE.
int isLJswitchOn; // use switch potential for LJ or not
int isEwaldOn; // Ewald summation electrostatic interactions
int fEwald_BC; // flag of boundary condition for ewald summation. 1-tinfoil; 2-vaccum
int fEwald_Dim; // flag of Ewald method dimension. 1, 2 or 3 dimension.
int isWolfOn; // wolf method for electrostatic interactions
double kappa, kappasq; // sqrt(alpha) in ewald summation. 
int KMAXX, KMAXY, KMAXZ; // ewald parameters
int KSQMAX; // ewald parameter
char coords_file[100]; // name of the coordinates file
int isNVTnh; // use nose hoover for NVT ensemble or not
int whichNH; // which nose hoover subroutine to use? usually 3 for molecule, 2 for atoms, see more details in nvtnh.c
int isSFon; // is solid-fluid interaction on
int sf_type; // solid-fluid type. for different nanotube potentials and future possible other materials
char atomname[natom_max][5];


int istep; // counter of step, current step
double utot;
double upot, ukin;
double uinter, uintra; // inter and intra molecular energy
double uvdw; // van der wall energy, LJ energy
double ubond, uangle, udih, uimp;
double uewald; // total ewald energy, refer to Frenkel and Smit, eq. 12.1.25
double ureal; // real part of ewald, term 3 in 12.1.25
double ufourier; // fourier part of ewald, term 1 in 12.1.25
double uself; // self interaction correction part of ewald, term 2 in 12.1.25
double uexcl; // excluding energy for ewald summation
double uvaccum; // vaccum boundary for ewald
double uGz0; // 1D ewald Gz=0 term
double LJswitch; // switch factor for LJ
double usflj; // solid-fluid LJ energy

int nfree; // freedom
double tinst; // instantaneous temperature

int nframe; // number of frames in the trajectory file

double TWOPI_LX, TWOPI_LY, TWOPI_LZ; // ewald
double Bfactor_ewald, Vfactor_ewald;
double twopi_over_3v; // constant for vaccum boundary

// variables for wolf method
double uwolf, uwolf_real, uwolf_con;
double wolfvcon1, wolfvcon2, wolffcon1, wolffcon2;

double roff2_minus_ron2_cube; // used for switch potential

// following variables are for nose hoover method
double Gts, Qts, vts, rts, AA;
double dt_outer2, dt_outer4, NRT;
double ukin_nhts, upot_nhts;
// frenkel and smit's nose hoover method
double qq, ps, gg, ss;
double delt_sqby2, delts_sqby2;
double vxo[natom_max], vyo[natom_max], vzo[natom_max]; 
double vxn[natom_max], vyn[natom_max], vzn[natom_max]; 
double bx[natom_max], by[natom_max], bz[natom_max];
double unhts;
double qqs, pss, ggs, sss;
double unhtss;

// my interpolations
double *Ene[nunique_atom_max]; 
double *dEx[nunique_atom_max], *dEy[nunique_atom_max], *dEz[nunique_atom_max];
double *dFxx[nunique_atom_max], *dFxy[nunique_atom_max], *dFxz[nunique_atom_max];
double *dFyx[nunique_atom_max], *dFyy[nunique_atom_max], *dFyz[nunique_atom_max];
double *dFzx[nunique_atom_max], *dFzy[nunique_atom_max], *dFzz[nunique_atom_max];
double uclx, ucly, uclz;
double grid_itvl_x, grid_itvl_y, grid_itvl_z;;
int ngrid_x, ngrid_y, ngrid_z, ngrid_total;
double xcenter, ycenter, zcenter;
double xmin, xmax, ymin, ymax, zmin, zmax;

// counters and accumulators
int icounter[num_counter_max];
double accumulator[num_counter_max][8];
// 0-4 for accumulator, 5-ave, 6-err, 7-fluc

/*
accumulator[0]   utot
accumulator[1]   upot
accumulator[2]   ukin
accumulator[3]   uinter
accumulator[4]   uintra
accumulator[5]   uvdw
accumulator[6]   ubond
accumulator[7]   uangle
accumulator[8]   udih
accumulator[9]   uimp
accumulator[10]   uewald
accumulator[11]   ureal
accumulator[12]   ufourier
accumulator[13]   uself
accumulator[14]   usflj
accumulator[15]   tinst
accumulator[16]   uvaccum

icounter[10]   number of average cycles

*/


#endif /* VARS_H */
