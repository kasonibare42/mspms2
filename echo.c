#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include "vars.h"

/// Echo the input data and initialized variables
int echo()
{
	int ii, jj, kk, ll;
	char szSimulation[100], szEnsemble[100];
	char szOutput[100];

	fprintf(fpouts, "==================================================\n");
	fprintf(fpouts, "Energies are for all molecules.\n");
	fprintf(fpouts, "Only total energy contains long range corrections.\n");
	fprintf(fpouts, "Pressure is with long range corrections.\n");
	fprintf(fpouts, "==================================================\n");
	fprintf(fpouts, "Units:\n");
	fprintf(fpouts,
			"  Energy(J/mol)\n  Temperature(K)\n  Velocity(m/s)\n  Distance(Angstrom)\n");
	fprintf(fpouts, "  Force(J/Angstrom/mol)\n");
	fprintf(fpouts, "  virial(J/mol)\n");
	fprintf(fpouts, "  Pressure(Pascal)\n");
	fprintf(fpouts, "==================================================\n");
	fprintf(fpouts, "System Name: %s\n", sysname);

	switch (what_simulation)
	{
	case md_run:
		strcpy(szSimulation, "MD");
		break;
	case hmc_run:
		strcpy(szSimulation, "HMC");
		break;
	default:
		strcpy(szSimulation, "Unknown");
		break;
	}
	switch (what_ensemble)
	{
	case nve_run:
		strcpy(szEnsemble, "NVE");
		break;
	case nvt_run:
		strcpy(szEnsemble, "NVT");
		break;
	case npt_run:
		strcpy(szEnsemble, "NPT");
		break;
	default:
		strcpy(szEnsemble, "Unknown");
		break;
	}
	fprintf(fpouts, "Running %s simulation (%d) with %s ensemble (%d).\n",
			szSimulation, what_simulation, szEnsemble, what_ensemble);
	fprintf(
			fpouts,
			"The simulation box is %lf long, %lf wide, %lf high, which is %lf in volume.\n",
			boxlx, boxly, boxlz, boxv);
	fprintf(
			fpouts,
			"The system contains %d specie(s), or %d molecule(s), or %d atom(s).\n",
			nspecie, nmole, natom);
	fprintf(fpouts, "There are totally %d bonds, %d angles, %d dihedrals, ",
			nbond, nangle, ndih);
	fprintf(fpouts, "%d improper dihedrals, and %d non-bonded pairs.\n", nimp,
			nnbp);
	fprintf(fpouts, "Total number of Constraints is %d.\n", nconstraint);
	for (ii=0; ii<nspecie; ii++)
	{
		fprintf(
				fpouts,
				"Specie %d has %d molecules.\n    Each molecule has %d atoms, ",
				ii, nmole_per_specie[ii], natom_per_mole[ii]);
		fprintf(fpouts, "%d bonds, %d angles, %d dihedrals, ",
				nbond_per_mole[ii], nangle_per_mole[ii], ndih_per_mole[ii]);
		fprintf(fpouts, "%d improper dihedrals, and %d non-bonded pairs.\n",
				nimp_per_mole[ii], nnbp_per_mole[ii]);
	}
	switch (fStart_option)
	{
	case new_run:
		strcpy(szOutput, "New Run");
		break;
	case continue_run:
		strcpy(szOutput, "Continue Run");
		break;
	case new_from_old:
		strcpy(szOutput, "New Run From Old Configuration");
		break;
	default:
		strcpy(szOutput, "Unknown");
		break;
	}
	fprintf(fpouts,
			"This is a %s, with %d total steps and %d equilibirum steps.\n",
			szOutput, nstep, nstep_eq);
	fprintf(fpouts, "    The averages are taken every %d steps.\n", nstep_ave);
	fprintf(fpouts, "    The logs are printed every %d steps.\n", nstep_print);
	fprintf(fpouts, "    The save happens every %d steps.\n", nstep_save);
	fprintf(fpouts, "    The snapshots are taken every %d steps.\n", nstep_ss);
	fprintf(fpouts, "    The trajectory is taken every %d steps.\n", nstep_trj);
	fprintf(fpouts, "Random seeds are set to be %d (ij) and %d (jk).\n", ij, jk);
	fprintf(fpouts, "Required Temperature is %lf.\n", treq);
	fprintf(fpouts, "Required Pressure is %lf.\n", preq);
	fprintf(
			fpouts,
			"Time step is %lf fs, with %d inner step(s), which is %lf fs for one inner step.\n",
			delt, nstep_inner, delts);
	fprintf(
			fpouts,
			"    Other time step related variables are: deltby2=%lf, deltsby2=%lf, dt_outer2=%lf, dt_outer4=%lf, dt_outer8=%lf\n",
			deltby2, deltsby2, dt_outer2, dt_outer4, dt_outer8);
	fprintf(
			fpouts,
			"The Cutoff for LJ interaction is %lf, for Coulomb interactions is %lf.\n    The squares are %lf and %lf, respectively.\n",
			rcutoff, rcutoffelec, rcutoffsq, rcutoffelecsq);
	fprintf(fpouts, "The 1,4 LJ interaction modifier is %lf.\n", f0);
	if (f0==0.0)
	{
		fprintf(fpouts,
				"    No 1,4 LJ interaction will be calculated, e.g. TraPPE.\n");
	}
	else if (f0==0.5)
	{
		fprintf(fpouts,
				"    The 1,4 LJ interactions are reduced by half, e.g. OPLS potential.\n");
	}
	else
	{
		fprintf(fpouts, "    The 1,4 LJ interactions are reduced.\n");
	}
	if (isLJswitchOn)
	{
		fprintf(
				fpouts,
				"Switch potential is enabled (%d) with Cuton = %lf (sqare is %lf).\n",
				isLJswitchOn, rcuton, rcutonsq);
	}
	else
	{
		fprintf(fpouts, "Switch potential is disabled (%d).\n", isLJswitchOn);
	}
	if (iChargeType != _NO_ELECTROSTATIC_INTERACTION)
	{
		fprintf(fpouts,
				"Electrostatic interaction calculations are enabled (%d).\n",
				iChargeType);
		if (isEwaldOn)
		{
			fprintf(
					fpouts,
					"    %dD Ewald summation is enabled for electrostatic interactions, ",
					fEwald_Dim);
			switch (fEwald_BC)
			{
			case ewald_bc_tinfoil:
				strcpy(szOutput, "Tinfoil");
				break;
			case ewald_bc_vacuum:
				strcpy(szOutput, "Vacuum");
				break;
			default:
				break;
			}
			fprintf(fpouts, "with %s boundary condition (%d).\n", szOutput,
					fEwald_BC);
			fprintf(fpouts, "    kappa=%lf, ", kappa);
			fprintf(fpouts, "kmaxx=%d, kmaxy=%d, kmaxz=%d, kmaxsq=%d\n", KMAXX,
					KMAXY, KMAXZ, KSQMAX);
		}
		else if (isWolfOn)
		{
			fprintf(fpouts, "    Wolf method is enabled.\n");
			fprintf(fpouts, "    kappa=%lf\n", kappa);
		}
		else if (isSimpleCoulomb)
		{
			fprintf(fpouts, "    Simple coulomb interaction is enabled.\n");
		}
		else
		{
			fprintf(fpouts, "    Unknown electrostatic interaction type.\n");
		}
	}
	else
	{
		fprintf(fpouts,
				"Electrostatic interaction calculations are disabled (%d).\n",
				iChargeType);
	}
	if (what_simulation == md_run)
	{
		if (what_ensemble == npt_run)
		{
			fprintf(fpouts, "MDNPT data section is required:\n");
			fprintf(
					fpouts,
					"    The mass of thermostat is %le. The Mass of barostat is %le. (Larger number means weaker coupling)\n",
					Qts, Qbs);
			fprintf(
					fpouts,
					"    Initialized variables: vts=%lf, vbs=%lf, rts=%lf, utsbs=%lf\n",
					vts, vbs, rts, utsbs);
		}
		else if (what_ensemble == nvt_run)
		{
			fprintf(fpouts, "MDNVT data section is required:\n");
			fprintf(fpouts, "    Required Temperature is %lf.\n", treq);
			fprintf(
					fpouts,
					"    The mass of outer Nose-Hoover thermostat is %le. The mass of inner Nose-Hoover thermostat is %le.\n",
					qq, qqs);
			fprintf(
					fpouts,
					"    Initialized variables: delt_sqby2=%lf, delts_sqby2=%lf, gg=%d, ss=%lf, ps=%lf, unhts=%lf, ",
					delt_sqby2, delts_sqby2, gg, ss, ps, unhts);
			fprintf(fpouts, "ggs=%d, sss=%lf, pss=%lf, unhtss=%lf\n", ggs, sss,
					pss, unhtss);
		}
	}
	else if (what_simulation == hmc_run)
	{
		fprintf(fpouts, "HMC data section is required:\n");
		fprintf(fpouts, "    %d MD step(s) are needed for one HMC step.\n");
		fprintf(
				fpouts,
				"    The probability of MD move is %lf, with required acceptance ratio %lf. The delt (%lf) will be adjusted every %d trials.\n",
				prob_cm, ratio_cm_req, delt, nstep_delt_adj_cycle);
		fprintf(
				fpouts,
				"    The probability of Volume Change move is %lf, with required acceptance ratio %lf. The delv (%lf) will be adjusted every %d trials.\n",
				prob_vc, ratio_vc_req, delv, nstep_delv_adj_cycle);
		fprintf(fpouts,
				"    The probability of insertion/deletion move is %lf.\n",
				prob_id);
	}

	if (sf_type != _NO_SF_POTENTIAL)
	{
		fprintf(fpouts, "Solid-fluid interaction calculation is enabled.\n");
		if (sf_type==nanotube_hypergeo)
		{
			fprintf(
					fpouts,
					"Structureless Nanotube is the absorbent. Hypergeo nanotube potential is enabled.\n");
			fprintf(
					fpouts,
					"There are %d Nanotubes, with sigma = %lf, epsilon = %lf\n",
					ntube, solid_sigma, solid_epsilon);
			for (ii=0; ii<ntube; ii++)
			{
				fprintf(fpouts, "   tube %d: xx=%lf  yy=%lf  radius=%lf\n", ii,
						hgntc_xx[ii], hgntc_yy[ii], hgnt_radius[ii]);
			}
		}
		else if (sf_type==nanotube_atom_explicit)
		{
			fprintf(fpouts,
					"Atom explicit Nanotube (or other absorbent) is enabled.\n");
			fprintf(fpouts, "There are totally %d atoms. ", solid_natom);
			if (fSolid_type==solid_uniform)
			{
				fprintf(
						fpouts,
						"The atoms are uniform, with sigma = %lf, epsilon = %lf, charge = %lf\n",
						*solid_sigma, *solid_epsilon, *solid_charge);
				for (ii=0; ii<solid_natom; ii++)
				{
					fprintf(fpouts, "%d %lf %lf %lf\n", ii, &solid_xx[ii],
							&solid_yy[ii], &solid_zz[ii]);
				}
			}
			else if (fSolid_type==solid_hetero)
			{
				fprintf(fpouts, "The atoms are heterogeneous.\n");
				for (ii=0; ii<solid_natom; ii++)
				{
					fprintf(fpouts, "%d %lf %lf %lf\n", ii, &solid_xx[ii],
							&solid_yy[ii], &solid_zz[ii]);
				}
			}
			else
			{
				fprintf(fpouts, "Unknown solid atom type.\n");
			}
		}
		else if (sf_type==nanotube_tasos)
		{
			fprintf(fpouts,
					"Tasos's interpolation of Nanotube potential is enabled.\n ");
		}
		else if (sf_type==nanotube_my_interp)
		{
			fprintf(fpouts,
					"YWang's interpolation of Nanotube potential is enabled.\n");
			fprintf(
					fpouts,
					"The size of the unit cell is %lf long, %lf wide, %lf high, ",
					uclx, ucly, uclz);
			fprintf(fpouts, "with center at (%lf, %lf, %lf) (x, y, z).\n",
					xcenter, ycenter, zcenter);
			fprintf(fpouts, "xmin=%lf  xmax=%lf\n", xmin, xmax);
			fprintf(fpouts, "ymin=%lf  ymax=%lf\n", ymin, ymax);
			fprintf(fpouts, "zmin=%lf  zmax=%lf\n", zmin, zmax);
		}
	}

	int iAtomPosInMole;
	fprintf(fpouts, "System configuration:\n");
	for (ii=0; ii<nspecie; ii++)
	{
		fprintf(fpouts, "Speice %d\n", ii);
		for (jj=specie_first_mole_idx[ii]; jj<specie_first_mole_idx[ii+1]; jj++)
		{
			fprintf(fpouts, " Mole %d, index %d in specie %d.\n", jj, jj
					-specie_first_mole_idx[ii], ii);
			fprintf(
					fpouts,
					"      ATOM PARAM: name, idx_global, idx_local, coordinates, soft-forces, rigid-forces\n");
			for (kk=mole_first_atom_idx[jj]; kk<mole_first_atom_idx[jj+1]; kk++)
			{
				iAtomPosInMole = kk - mole_first_atom_idx[jj];
				fprintf(fpouts,
						"  %s %d %d %lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n",
						atomname, kk, iAtomPosInMole, xx[kk], yy[kk], zz[kk],
						fxl[kk], fyl[kk], fzl[kk], fxs[kk], fys[kk], fzs[kk]);
				fprintf(fpouts, "      Excluding pairs: ");
				for (ll=pointexcl_atom[ii][iAtomPosInMole]; ll
						<pointexcl_atom[ii][iAtomPosInMole+1]; ll++)
				{
					fprintf(fpouts, "%d ", excllist[ll]+kk);
				}
				fprintf(fpouts, "\n");
			}
			fprintf(
					fpouts,
					"      BOND PARAM: idx_global, idx_local, bond_composition, type, Kb, Req, alpha\n");
			for (kk=mole_first_bond_idx[jj]; kk<mole_first_bond_idx[jj+1]; kk++)
			{
				fprintf(fpouts, "  Bond %d %d, %d-%d, %d, %lf, %lf, %lf\n", kk,
						kk -mole_first_bond_idx[jj], bond_idx[kk][0],
						bond_idx[kk][1], bond_type[kk], Kb[kk], Req[kk],
						alpha[kk]);
			}
			fprintf(
					fpouts,
					"      ANGLE PARAM: idx_global, idx_local, angle_composition, type, Ktheta, ThetaEq, angle_para_3, angle_para_4, angle_para_5\n");
			for (kk=mole_first_angle_idx[jj]; kk<mole_first_angle_idx[jj+1]; kk++)
			{
				fprintf(
						fpouts,
						"  Angle %d %d, %d-%d-%d, %d, %lf, %lf, %lf, %lf, %lf\n",
						kk, kk-mole_first_angle_idx[jj], angle_idx[kk][0],
						angle_idx[kk][1], angle_idx[kk][2], angle_type[kk],
						Ktheta[kk], Thetaeq[kk], agl_para_3[kk],
						agl_para_4[kk], agl_para_5[kk]);
			}
			fprintf(
					fpouts,
					"      DIHEDRAL PARAM: idx_global, idx_local, dihedral_composition, type, c1, c2, c3, c4\n");
			for (kk=mole_first_dih_idx[jj]; kk<mole_first_dih_idx[jj+1]; kk++)
			{
				fprintf(fpouts,
						"  Dih %d %d, %d-%d-%d-%d, %d, %lf, %lf, %lf, %lf\n",
						kk, kk-mole_first_dih_idx[jj], dih_idx[kk][0],
						dih_idx[kk][1], dih_idx[kk][2], dih_idx[kk][3],
						dih_type[kk], c1[kk], c2[kk], c3[kk], c4[kk]);
			}
			fprintf(
					fpouts,
					"      IMPROPER PARAM: idx_global, idx_local, improper_composition, type, Komega, Omega0\n");
			for (kk=mole_first_imp_idx[jj]; kk<mole_first_imp_idx[jj+1]; kk++)
			{
				fprintf(fpouts, "  Imp %d %d, %d-%d-%d-%d, %d, %lf, %lf\n", kk,
						kk-mole_first_imp_idx[jj], imp_idx[kk][0],
						imp_idx[kk][1], imp_idx[kk][2], imp_idx[kk][3],
						imp_type[kk], komega[kk], omega0[kk]);
			}
			fprintf(
					fpouts,
					"      NON-BONDED PARAM: idx_global, idx_local, non-bonded_composition\n");
			for (kk=mole_first_nbp_idx[jj]; kk<mole_first_nbp_idx[jj+1]; kk++)
			{
				fprintf(fpouts,"  Nbp %d %d", kk, kk-mole_first_nbp_idx[jj], nbp_idx[kk][0], nbp_idx[kk][1]);
			}
		}
	}

	fprintf(fpouts, "The initial energies etc.:\n");
	fprintf(fpouts, "utot=%lf\n", uinter+uintra+ukin); // utot
	fprintf(fpouts, "upot=%lf\n", uinter+uintra); // upot
	fprintf(fpouts, "ukin=%lf\n", ukin);
	fprintf(fpouts, "uinter=%lf\n", uinter);
	fprintf(fpouts, "uintra=%lf\n", uintra);
	fprintf(fpouts, "uvdw=%lf\n", uvdw);
	fprintf(fpouts, "uer_vdw=%lf\n", uvdw-unbp_vdw);
	fprintf(fpouts, "unbp_vdw=%lf\n", unbp_vdw);
	fprintf(fpouts, "ubond=%lf\n", ubond);
	fprintf(fpouts, "uangle=%lf\n", uangle);
	fprintf(fpouts, "udih=%lf\n", udih);
	fprintf(fpouts, "uimp=%lf\n", uimp);
	fprintf(fpouts, "uewald=%lf\n", uewald);
	fprintf(fpouts, "ureal=%lf\n", ureal);
	fprintf(fpouts, "ufourier=%lf\n", ufourier);
	fprintf(fpouts, "uself=%lf\n", uself);
	fprintf(fpouts, "uexcl=%lf\n", uexcl);
	fprintf(fpouts, "uvacuum=%lf\n", uvacuum);
	fprintf(fpouts, "uGz0=%lf\n", uGz0);

	fprintf(fpouts, "usflj=%lf\n", usflj);
	fprintf(fpouts, "uwolf=%lf\n", uwolf);
	fprintf(fpouts, "uwolf_real=%lf\n", uwolf_real);
	fprintf(fpouts, "uwolf_con=%lf\n", uwolf_con);
	fprintf(fpouts, "ucoulomb=%lf\n", ucoulomb);

	fprintf(fpouts, "tinst=%lf\n", tinst);

	fprintf(fpouts, "virial=%lf\n", virial_inter+virial_intra);
	fprintf(fpouts, "virial_inter=%lf\n", virial_inter);
	fprintf(fpouts, "virial_intra=%lf\n", virial_intra);
	fprintf(fpouts, "pideal=%lf\n", pideal=natom/(boxlx*boxly*boxlz) *tinst
			*kb_1e30);
	// pljlrc already calculated in initialization part
	pinst = pideal + (virial_inter+virial_intra)*virial_to_pressure /(boxlx
			*boxly*boxlz) + pljlrc;
	fprintf(fpouts, "pressure=%lf\n", pinst);

	fprintf(fpouts, "uljlrc=%lf\npljlrc=%lf\n", uljlrc, pljlrc);

	fflush(fpouts);

	return 0;
}
