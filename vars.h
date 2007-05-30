#ifndef VARS_H
#define VARS_H   

#define INPUT		"in.mspms"
#define OUTPUT		"out.mspms"
#define CONFIG		"cfg.mspms"

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


int nspecie, nmole, natom; /* total number of molecules, atoms, species */
float xx[natom_max], yy[natom_max], zz[natom_max]; /* position */
float vx[natom_max], vy[natom_max], vz[natom_max]; /* velocity */
float fxl[natom_max], fyl[natom_max], fzl[natom_max]; /* inter force */
float fxs[natom_max], fys[natom_max], fzs[natom_max]; /* intra force */
float aw[natom_max], epsilon[natom_max], sigma[natom_max], charge[natom_max]; /* atom weight etc. */
int isghost[natom_max]; /* flag of ghost atom, including ghost LJ or ghost electrostatic */
int tasostype[natom_max]; /* atom type for tasos interpolation */

int nbond, nangle, ndih, nimp, nnbp;
int bond_idx[nbond_max][2];
int bond_type[nbond_max]; 
float Kb[nbond_max], Req[nbond_max], alpha[nbond_max];
int angle_idx[nangle_max][3];
int angle_type[nangle_max];
float Ktheta[nbond_max], Thetaeq[nbond_max];
int dih_idx[ndih_max][4];
int dih_type[nbond_max];
float c1[nbond_max], c2[nbond_max], c3[nbond_max], c4[nbond_max];
int imp_idx[nimp_max][4];
int nbp_idx[nnbp_max][2];
int excllist[exclude_max];
int pointexcl[natom_max+1]; // the index of exclude list starts and ends

char sysname[200];

FILE *fpins, *fpouts, *fpcfg;

int nstep, nstep_eq, nstep_start; /* nstep_start is the starting step for continue runs */
int nstep_inner; /* inner step for multi time step */
float dt; /* time step */
int ij, jk;
float treq, preq; /* input required temperature, pressure */
float boxlx, boxly, boxlz; /* box size */
int nconstraint; // constraint, 3 for periodic, 6 for aperiodic
float rcutoff, rcutoffsq;


float upot, ukin;
float uvdw; // van der wall energy, LJ energy
float ubond, uangle, udih, uimp;

int nfree; // freedom
float tinst; // instantaneous temperature



#endif /* VARS_H */
