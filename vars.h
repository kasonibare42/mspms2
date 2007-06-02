#ifndef VARS_H
#define VARS_H   

#define INPUT		"in.mspms"
#define OUTPUT		"out.mspms"
#define CONFIG		"cfg.mspms"
#define LOGFILE		"log.mspms"
#define SNAPSHOT	"ss.mspms"
#define MOVIE		"trj.mspms"
#define SAVEFILE	"sav.mspms"

#define true		1
#define false		0
#define not_ghost	0
#define lj_ghost	1
#define all_ghost	2

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

#define nspecie_max	3
#define nmole_max	1000
#define natom_max	4000
#define nbond_max	2000
#define nangle_max	2000
#define ndih_max	2000
#define nimp_max	2000
#define nnbp_max	2000
#define exclude_max	4000
#define Rgas		8.314472 /* J/mol/K */
#define rRgas		0.120272219 /* reciprocal of Rgas */
#define const_columb	1389355.1051  // unit is J/mol. Na*(1.602177e-19)^2/4/PI/epsilon0/(1.0e-10)
#define pi		3.141592653589793


int nspecie, nmole, natom; /* total number of molecules, atoms, species */
double xx[natom_max], yy[natom_max], zz[natom_max]; /* position */
double vx[natom_max], vy[natom_max], vz[natom_max]; /* velocity */
double fxl[natom_max], fyl[natom_max], fzl[natom_max]; /* inter force */
double fxs[natom_max], fys[natom_max], fzs[natom_max]; /* intra force */
double aw[natom_max], epsilon[natom_max], sigma[natom_max], charge[natom_max]; /* atom weight etc. */
int isghost[natom_max]; /* flag of ghost atom, including ghost LJ or ghost electrostatic */
int tasostype[natom_max]; /* atom type for tasos interpolation */

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
int nbp_idx[nnbp_max][2];
int excllist[exclude_max];
int pointexcl[natom_max+1]; // the index of exclude list starts and ends

char sysname[200];

FILE *fpins, *fpouts, *fpcfg, *fplog;
FILE *fpss, *fptrj, *fpsave;

int nstep, nstep_eq, nstep_start; /* nstep_start is the starting step for continue runs */
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
int isEwaldOn; // switch for Ewald summation electrostatic interactions
int isWolfOn; // switch of wolf method for electrostatic interactions
double kappa; // sqrt(alpha) in ewald summation. 
int KMAXX, KMAXY, KMAXZ; // ewald parameters
int KSQMAX; // ewald parameter
char coords_file[100]; // name of the coordinates file


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

int nfree; // freedom
double tinst; // instantaneous temperature

double TWOPI_LX, TWOPI_LY, TWOPI_LZ; // ewald
double Bfactor_ewald, Vfactor_ewald;



#endif /* VARS_H */
