#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string.h>
#include <math.h>
#include "inirw.h"
#include "clsdefs.h"
#include "atom.h"
#include "mole.h"
#include "specie.h"
#include "nanotube.h"
#include "box.h"

#define Rgas 8.314472 // J/mol/K

#define bond_none               0
#define bond_harmonic           1
#define bond_morse              2
#define bond_fene               3

#define angle_none              0
#define angle_harmonic          1
#define angle_TRwater           2  // Toukan and Rahman water potentials

#define dih_none                0
#define dih_opls_cosin          1
#define dih_charmm              2 // charmm type dihedral potential

#define imp_none                0
#define imp_charmm              1 // charmm type improper potential

#define nanotube_polynomial     1
#define nanotube_atom_explicit  2
#define nanotube_tasos          3

#define pi              3.141592653589793

using namespace std;

char sysname[200];
int nspecie;
int nmole_per_specie[_MAX_SPECIE];
int nmole;
char potential_form[_MAX_SPECIE][200];
int bond_type[_MAX_SPECIE];
int angle_type[_MAX_SPECIE];
int dih_type[_MAX_SPECIE];
int imp_type[_MAX_SPECIE];
int natom_per_mole[_MAX_SPECIE];
int natom;
int natom_per_specie[_MAX_SPECIE];
SIMSPECIE SampleSpecie[10];
PSIMMOLE pMole;
PSIMATOM pAtom;
ifstream ins;
ofstream outs;
int nbond, nangle, ndih, nimp, nnbp;

int main(int argc, char *argv[])
{
	int ii, jj, kk;
	INISET objPotentials("potentials.ini");
	double data[10];
	PBOND pBond;
	PANGLE pAngle;
	PDIHEDRAL pDih;
	PNONPAIR pNonpair;
	bool isExist;
	char szPotential_model[200];
	char szBond[200];
	char szAngle[200];
	char szDihedral[200];
	int iatom_counter, imole_counter;

	ins.open(argv[1],ios::in);
	if (ins==NULL)
	{
		cout << "Error: need input file." << endl;
		exit(1);
	}

	natom = 0;
	nmole = 0;
	nbond = nangle = ndih = nimp = nnbp = 0;

	ins >> sysname;
	ins >> nspecie;
	ins.ignore(200,'\n');
	for (ii=0;ii<nspecie;ii++)
	{
		// (SampleSpecie+ii)->init();
		pMole = &((SampleSpecie+ii)->sample_mole);
		ins >> nmole_per_specie[ii];
		ins >> natom_per_mole[ii];
		ins >> pMole->std_bonding_dist;
		ins >> potential_form[ii];
		ins >> bond_type[ii];
		ins >> angle_type[ii];
		ins >> dih_type[ii];
		ins >> imp_type[ii];
		nmole += nmole_per_specie[ii];
		natom_per_specie[ii] = nmole_per_specie[ii]*natom_per_mole[ii];
		natom += natom_per_specie[ii];
		pMole->natom = natom_per_mole[ii];
		pMole->atom_list = new SIMATOM[natom_per_mole[ii]];
		pMole->pSample_mole = pMole;
		pMole->fNonpair_mode = _1_5_NONPAIR;
		ins.ignore(200,'\n'); // ignore the rest of the first line
	}

	for (ii=0;ii<nspecie;ii++)
	{
		pMole = &((SampleSpecie+ii)->sample_mole);
		pMole->weight = 0.0;
		for (jj=0;jj<pMole->natom;jj++)
		{
			pAtom = pMole->atom_list + jj;
			pAtom->id = jj;
			pAtom->pNbr_list = SampleSpecie[ii].nbr_list[jj];
			ins >> pAtom->name >> pAtom->xx >> pAtom->yy >> pAtom->zz
				>> pAtom->type >> pAtom->fGhost >> pAtom->tasostype;
			ins.ignore(200,'\n');
			objPotentials.ReadDoubleMany("Atom Basic Properties",pAtom->type,4,data);
			pAtom->weight = 0.001*data[0];
			pAtom->sigma = data[1];
			pAtom->epsilon = data[2]*Rgas; // K to J/mol
			pAtom->charge = data[3];

			pMole->weight += pAtom->weight;
		}
		pMole->iCenter_atom = 0; // always use first atom as the center of mass atom
		pMole->pXx = &((pMole->atom_list + pMole->iCenter_atom)->xx);
		pMole->pYy = &((pMole->atom_list + pMole->iCenter_atom)->yy);
		pMole->pZz = &((pMole->atom_list + pMole->iCenter_atom)->zz);
		// set ee,ff,gg
		pMole->cal_all_atom_xyz2efg();
		// analyze the molecular structure
		pMole->analyze_mole_structure();

		nbond += pMole->nbond*nmole_per_specie[ii];
		nangle += pMole->nangle*nmole_per_specie[ii];
		ndih += pMole->ndih*nmole_per_specie[ii];
		nimp += pMole->nimp*nmole_per_specie[ii];
		nnbp += pMole->nnonpair*nmole_per_specie[ii];

		sprintf(szPotential_model, "Potential model: %s", potential_form[ii]);
		for (jj=0;jj<pMole->nbond;jj++)
		{
			pBond = pMole->pBond_list + jj;
			isExist = false;
			sprintf(szBond, "%s-%s",(pMole->atom_list+pBond->comp_id[0])->type,
					(pMole->atom_list+pBond->comp_id[1])->type);
			if (objPotentials.SearchData(szPotential_model, szBond))
				isExist = true;
			else
			{
				sprintf(szBond, "%s-%s", (pMole->atom_list+pBond->comp_id[1])->type,
						(pMole->atom_list+pBond->comp_id[0])->type);
				if (objPotentials.SearchData(szPotential_model, szBond))
					isExist = true;
			}
			if (isExist)
			{
				switch(bond_type[ii])
				{
					case bond_none:
						pBond->kr = 0.0;
						pBond->req = 0.0;
						pBond->alpha = 0.0;
						break;
					case bond_harmonic:
						objPotentials.ReadDoubleMany(szPotential_model,szBond,2,data);
						pBond->kr = data[0]*Rgas; // convert from K to J/mol
						pBond->req = data[1];
						pBond->alpha = 0.0; // dummy
						break;
					case bond_morse:
						objPotentials.ReadDoubleMany(szPotential_model,szBond,3,data);
						pBond->kr = data[0]*Rgas; // convert from K to J/mol
						pBond->req = data[1];
						pBond->alpha = data[2]; 
						break;
					case bond_fene:
						objPotentials.ReadDoubleMany(szPotential_model,szBond,2,data);
						pBond->kr = data[0]*Rgas; // convert from K to J/mol
						pBond->req = data[1];
						pBond->alpha = 0.0; // dummy
						break;
				}
			}
			else
			{
				cout << "Error: cannot find parameters for bond " << szBond << endl;
				exit(1);
			}
		}

		for (jj=0;jj<pMole->nangle;jj++)
		{
			pAngle = pMole->pAngle_list + jj;
			isExist = false;
			sprintf(szAngle, "%s-%s-%s", (pMole->atom_list+pAngle->comp_id[0])->type,
					(pMole->atom_list+pAngle->comp_id[1])->type,
					(pMole->atom_list+pAngle->comp_id[2])->type);
			if (objPotentials.SearchData(szPotential_model, szAngle)) // found
				isExist = true;
			else
			{
				sprintf(szAngle, "%s-%s-%s", (pMole->atom_list+pAngle->comp_id[2])->type,
						(pMole->atom_list+pAngle->comp_id[1])->type,
						(pMole->atom_list+pAngle->comp_id[0])->type);
				if (objPotentials.SearchData(szPotential_model, szAngle)) // found
					isExist = true;
			}
			if (isExist)
			{
				switch(angle_type[ii])
				{
					case angle_none:
						pAngle->k_theta = 0.0;
						pAngle->theta_eq = 0.0;
						pAngle->agl_para_3 = 0.0;
						pAngle->agl_para_4 = 0.0;
						pAngle->agl_para_5 = 0.0;
						break;
					case angle_harmonic:
						objPotentials.ReadDoubleMany(szPotential_model, szAngle, 2, data);
						pAngle->k_theta = data[0]*Rgas; // K to J/mol
						pAngle->theta_eq = data[1]*pi/180.0; // convert from degree to radian
						pAngle->agl_para_3 = 0.0; // dummy
						pAngle->agl_para_4 = 0.0;
						pAngle->agl_para_5 = 0.0;
						break;
					case angle_TRwater:
						objPotentials.ReadDoubleMany(szPotential_model, szAngle, 5, data);
						pAngle->k_theta = data[0]*Rgas; // K to J/mol
						pAngle->theta_eq = data[1]*Rgas; // it is k_r_theta
						pAngle->agl_para_3 = data[2]*Rgas; // k_r_rprime
						pAngle->agl_para_4 = data[3];
						pAngle->agl_para_5 = data[4];
						break;
				}
			}
			else
			{
				cout << "Error: cannot find parameters for angle " << szAngle << endl;
				exit(1);
			}
		}


		for (jj=0;jj<pMole->ndih;jj++)
		{
			pDih = pMole->pDih_list + jj;
			isExist = false;
			sprintf(szDihedral, "%s-%s-%s-%s", (pMole->atom_list+pDih->comp_id[0])->type,
					(pMole->atom_list+pDih->comp_id[1])->type,
					(pMole->atom_list+pDih->comp_id[2])->type,
					(pMole->atom_list+pDih->comp_id[3])->type);
			if (objPotentials.SearchData(szPotential_model, szDihedral)) // found
				isExist = true;
			else
			{
				sprintf(szDihedral, "%s-%s-%s-%s", (pMole->atom_list+pDih->comp_id[3])->type,
						(pMole->atom_list+pDih->comp_id[2])->type,
						(pMole->atom_list+pDih->comp_id[1])->type,
						(pMole->atom_list+pDih->comp_id[0])->type);
				if (objPotentials.SearchData(szPotential_model, szDihedral)) // found
					isExist = true;
			}
			if (isExist)
			{
				switch(dih_type[ii])
				{
					case dih_none:
						pDih->v1 = 0.0; // dummy
						pDih->v2 = 0.0; // dummy
						pDih->v3 = 0.0; // dummy
						pDih->v4 = 0.0; // dummy
						break;
					case dih_opls_cosin:
						objPotentials.ReadDoubleMany(szPotential_model, szDihedral, 4, data);
						pDih->v1 = data[0]*Rgas;
						pDih->v2 = data[1]*Rgas;
						pDih->v3 = data[2]*Rgas;
						pDih->v4 = data[3]*Rgas;
						break;
					case dih_charmm:
						break;
				}
			}
			else
			{
				cout << "Error: cannot find parameters for dihedral " << szDihedral << endl;
				exit(1);
			}
		}


		// print out molecule structure
		pMole->print_structure();

	}

	ins.close();


	outs.open("cfg.mspms",ios::out);
	outs.setf(ios::fixed);
	outs.precision(4);
	// title
	outs << sysname << endl;
	outs << nspecie << " species" << endl;

	// output loop through species
	for (ii=0;ii<nspecie;ii++)
	{
		outs << "Specie" << ii << endl;
		outs << nmole_per_specie[ii] << " molecules" << endl; 
		pMole = &((SampleSpecie+ii)->sample_mole);
		outs << pMole->natom << " atoms per molecule" << endl;
		// output atoms
		for (kk=0;kk<pMole->natom;kk++)
		{
			pAtom = pMole->atom_list + kk;
			outs << kk << ' '
				<< pAtom->name << ' '
				<< pAtom->xx << ' '
				<< pAtom->yy << ' '
				<< pAtom->zz << ' '
				<< pAtom->ee << ' '
				<< pAtom->ff << ' '
				<< pAtom->zz << ' '
				<< pAtom->weight << ' '
				<< pAtom->epsilon << ' '
				<< pAtom->sigma << ' '
				<< pAtom->charge << ' '
				<< pAtom->fGhost << ' '
				<< pAtom->tasostype << endl;
		}

		// bonds
		outs << pMole->nbond << "  bonds per molecule" << endl;
		for (kk=0;kk<pMole->nbond;kk++)
		{
			pBond = pMole->pBond_list + kk;
			outs << pBond->comp_id[0] << ' '
				<< pBond->comp_id[1] << ' '
				<< bond_type[ii] << ' '
				<< pBond->kr << ' '
				<< pBond->req << ' '
				<< pBond->alpha << endl;
		}


		// angles
		outs << pMole->nangle << "  angles per molecule" << endl;
		for (kk=0;kk<pMole->nangle;kk++)
		{
			pAngle = pMole->pAngle_list + kk;
			outs << pAngle->comp_id[0] << ' '
				<< pAngle->comp_id[1] << ' '
				<< pAngle->comp_id[2] << ' '
				<< angle_type[ii] << ' '
				<< pAngle->k_theta << ' '
				<< pAngle->theta_eq << ' '
				<< pAngle->agl_para_3 << ' '
				<< pAngle->agl_para_4 << ' '
				<< pAngle->agl_para_5 << endl;
		}

		// dihedrals
		outs << pMole->ndih << "  dihedrals per molecule" << endl;
		for (kk=0;kk<pMole->ndih;kk++)
		{
			pDih = pMole->pDih_list + kk;
			outs << pDih->comp_id[0] << ' '
				<< pDih->comp_id[1] << ' '
				<< pDih->comp_id[2] << ' '
				<< pDih->comp_id[3] << ' '
				<< dih_type[ii] << ' '
				<< pDih->v1 << ' '
				<< pDih->v2 << ' '
				<< pDih->v3 << ' '
				<< pDih->v4 << endl;
		}

		// impropers
		outs << pMole->nimp << "  impropers per molecule" << endl;
		// nonbonded
		outs << pMole->nnonpair << "   nonbonded per molecule" << endl;;
		for (kk=0;kk<pMole->nnonpair;kk++)
		{
			pNonpair = pMole->pNonpair_list + kk;
			outs << pNonpair->comp_id[0] << ' '
				<< pNonpair->comp_id[1] << endl;
		}



	}

	cout << "Structures and cfg.mspms file need to be checked before use." << endl;

}



