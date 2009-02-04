#ifndef FUNCS_H_
#define FUNCS_H_

int init_vars();
int velinit();
int rezero();
void get_specie_and_relative_atom_id(int abs_atom_id, int *specie_id, int *sample_atom_id);
int init_sf_hypergeo();
int init_sf_atom_explicit();
int init_tasos_grid();
int init_my_interp();
int init_npt_respa();
int init_nvt();
int init_hmc();
int calculate_ljlrc();
int fnValidateInput();

int frclong();
int ljfrc(double rijsq, double sigmaij, double epsilonij, double *uij, double *fij, double *uijshift);
int ewald_real_frc(double rijsq, double chargeij, double *uij, double *fij);
int wolf_real_frc(double rijsq, double chargeij, double *uij, double *fij);
int xmole_coulomb_frc(int iSpecie, int iMole, int iabs, int iAtom, int jSpecie, int jMole, int jabs, int jAtom);
int xatom_coulomb_frc(double rijsq, double chargeij, double *uij, double *fij);
int cal_com_and_efg_one(int iSpecie, int iMole, int iabs, int iAtom);
int reconstruct_from_com_one(int iMole, int iabs, int iAtom);
int frcshort();
int bndfrc(int iSpecie, int iBond, int iabs, double *uij, double *virial_ij);
int aglfrc(int iSpecie, int iAngle, int iabs, double *uij, double *virial_ij);
int dihfrc(int iSpecie, int iDih, int iabs, double *uij, double *virial_ij);
int impfrc(int iSpecie, int iImp, int iabs, double *uij, double *virial_ij);
int cal_sf_hypergeo(int ii, int iSpecie, int iAtom, double *uij, double *fij);
int cal_sf_atom_explicit(int ii, int iSpecie, int iAtom, double *uij, double *fij);
int get_values_from_grid(double fxx, double fyy, double fzz, int type, double *usf, double *fsf);
int call_tasos_forces(int itype, double fxx, double fyy, double fzz, double *usflj_tasos, double *tasos_force);

bool md_move(int (*pfn_md_scheme)(), int nStepMD);

int vver();
int averages();
int printit();
int calres();
int loadit();
int saveit();
int printit();
int snapshot();
int trajectory();
int echo();
int vver_nh_3();
int npt_respa();
int hmc();
int fnValidateInit();
int volume_change();
int collect_aves();
int fnMetalClusterFF_Cu();
int fnMetalClusterFF_Ag();
int fnInsDelMole();
int opening();
int ending();
int readins();
int init_siman();
int siman();

#endif /*FUNCS_H_*/
