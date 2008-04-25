/**
 * Project: mspms2
 * File: echo.c
 * 
 * Copyright (C) 2008    Yang Wang <ywangd@gmail.com>
 * Created @ 2007
 * Modified @ Apr 24, 2008
 * 
 * Description:
 * 
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include "mspms2.h"

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
	fprintf(fpouts, "Reduced Units:\n");
	fprintf(fpouts, "  Energy(%lf K)\n", epsilon_base);  
	fprintf(fpouts, "  Temperature(%lf K)\n", epsilon_base);
	fprintf(fpouts, "  Separation(%lf Angstrom)\n", sigma_base);
	fprintf(fpouts, "  Mass (%lf*E-27 kg/molecule)\n", mass_base);
	fprintf(fpouts, "  Pressure(%lf bar)\n", pressure_base);
	fprintf(fpouts, "  Time (%lf fs)\n", time_base);
	fprintf(fpouts, "  Virial(%lf K)\n", epsilon_base);
	fprintf(fpouts, "  Force (%lf K/Angstrom)\n", force_base);
	fprintf(fpouts, "  Velocity (%lf Angstrom/fs)\n", velocity_base);
	fprintf(fpouts, "==================================================\n");
	fprintf(fpouts, "System Name: %s\n", sysname);

	switch (what_simulation)
	{
	case MOLECULAR_DYNAMICS:
		strcpy(szSimulation, "MD");
		break;
	case HYBRID_MONTE_CARLO:
		strcpy(szSimulation, "HMC");
		break;
	default:
		strcpy(szSimulation, "Unknown");
		break;
	}
	switch (what_ensemble)
	{
	case NVE:
		strcpy(szEnsemble, "NVE");
		break;
	case NVT:
		strcpy(szEnsemble, "NVT");
		break;
	case NPT:
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
	/*
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
				ii, nmole_per_specie[ii], sample_natom_per_mole[ii]);
		fprintf(fpouts, "%d bonds, %d angles, %d dihedrals, ",
				sample_nbond_per_mole[ii], sample_nangle_per_mole[ii], sample_ndih_per_mole[ii]);
		fprintf(fpouts, "%d improper dihedrals, and %d non-bonded pairs.\n",
				sample_nimp_per_mole[ii], sample_nnbp_per_mole[ii]);
	}
	*/
	switch (iStart_option)
	{
	case NEW:
		strcpy(szOutput, "New Run");
		break;
	case CONTINUE:
		strcpy(szOutput, "Continue Run");
		break;
	case CONFIG_ONLY:
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
	if (iChargeType != ELECTROSTATIC_NONE)
	{
		fprintf(fpouts,
				"Electrostatic interaction calculations are enabled (%d).\n",
				iChargeType);
		if (iChargeType==ELECTROSTATIC_EWALD)
		{
			fprintf(
					fpouts,
					"    %dD Ewald summation is enabled for electrostatic interactions, ",
					fEwald_Dim);
			switch (fEwald_BC)
			{
			case EWALD_BC_TINFOIL:
				strcpy(szOutput, "Tinfoil");
				break;
			case EWALD_BC_VACUUM:
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
		else if (iChargeType==ELECTROSTATIC_WOLF)
		{
			fprintf(fpouts, "    Wolf method is enabled.\n");
			fprintf(fpouts, "    kappa=%lf\n", kappa);
		}
		else if (iChargeType==ELECTROSTATIC_SIMPLE_COULOMB)
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
	if (what_simulation == MOLECULAR_DYNAMICS)
	{
		if (what_ensemble == NPT)
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
		else if (what_ensemble == NVT)
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
	else if (what_simulation == HYBRID_MONTE_CARLO)
	{
		fprintf(fpouts, "HMC data section is required:\n");
		fprintf(fpouts, "    %d MD step(s) are needed for one HMC step.\n");
		fprintf(
				fpouts,
				"    The probability of MD move is %lf, with required acceptance ratio %lf. The delt (%lf) will be adjusted every %d trials.\n",
				pdisp, rreq_disp, delt, nstep_delt_adj_cycle);
		fprintf(
				fpouts,
				"    The probability of Volume Change move is %lf, with required acceptance ratio %lf. The delv (%lf) will be adjusted every %d trials.\n",
				pvolm, rreq_volm, delv, nstep_delv_adj_cycle);
		fprintf(fpouts,
				"    The probability of insertion/deletion move is %lf.\n",
				pmake);
	}

	if (iSF_type != SF_NONE)
	{
		fprintf(fpouts, "Solid-fluid interaction calculation is enabled.\n");
		if (iSF_type==SF_NANOTUBE_HYPERGEO)
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
		else if (iSF_type==SF_NANOTUBE_ATOM_EXPLICIT)
		{
			fprintf(fpouts,
					"Atom explicit Nanotube (or other absorbent) is enabled.\n");
			fprintf(fpouts, "There are totally %d atoms. ", solid_natom);
			if (fSolid_type==SOLID_UNIFORM)
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
			else if (fSolid_type==SOLID_HETERO)
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
		else if (iSF_type==SF_NANOTUBE_TASOS)
		{
			fprintf(fpouts,
					"Tasos's interpolation of Nanotube potential is enabled.\n ");
		}
		else if (iSF_type==SF_NANOTUBE_MY_INTERP)
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
	/*
	fprintf(fpouts, "System configuration:\n");
	for (ii=0; ii<nspecie; ii++)
	{
		fprintf(fpouts, "Speice %d\n", ii);
		for (jj=specie_first_mole_idx[ii]; jj<specie_last_mole_idx[ii]; jj++)
		{
			fprintf(fpouts, " Mole %d, index %d in specie %d.\n", jj, jj
					-specie_first_mole_idx[ii], ii);
			fprintf(
					fpouts,
					"      ATOM PARAM: name, idx_global, idx_local, coordinates, soft-forces, rigid-forces\n");
			for (kk=mole_first_atom_idx[jj]; kk<mole_last_atom_idx[jj]; kk++)
			{
				iAtomPosInMole = kk - mole_first_atom_idx[jj];
				fprintf(fpouts,
						"  %s %d %d %lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n",
						atomname[kk], kk, iAtomPosInMole, xx[kk], yy[kk], zz[kk],
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
			for (kk=mole_first_bond_idx[jj]; kk<mole_last_bond_idx[jj]; kk++)
			{
				fprintf(fpouts, "  Bond %d %d, %d-%d, %d, %lf, %lf, %lf\n", kk,
						kk -mole_first_bond_idx[jj], bond_idx[kk][0],
						bond_idx[kk][1], bond_type[kk], Kb[kk], Req[kk],
						alpha[kk]);
			}
			fprintf(
					fpouts,
					"      ANGLE PARAM: idx_global, idx_local, angle_composition, type, Ktheta, ThetaEq, angle_para_3, angle_para_4, angle_para_5\n");
			for (kk=mole_first_angle_idx[jj]; kk<mole_last_angle_idx[jj]; kk++)
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
			for (kk=mole_first_dih_idx[jj]; kk<mole_last_dih_idx[jj]; kk++)
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
			for (kk=mole_first_imp_idx[jj]; kk<mole_last_imp_idx[jj]; kk++)
			{
				fprintf(fpouts, "  Imp %d %d, %d-%d-%d-%d, %d, %lf, %lf\n", kk,
						kk-mole_first_imp_idx[jj], imp_idx[kk][0],
						imp_idx[kk][1], imp_idx[kk][2], imp_idx[kk][3],
						imp_type[kk], komega[kk], omega0[kk]);
			}
			fprintf(
					fpouts,
					"      NON-BONDED PARAM: idx_global, idx_local, non-bonded_composition\n");
			for (kk=mole_first_nbp_idx[jj]; kk<mole_last_nbp_idx[jj]; kk++)
			{
				fprintf(fpouts,"  Nbp %d %d", kk, kk-mole_first_nbp_idx[jj], nbp_idx[kk][0], nbp_idx[kk][1]);
			}
		}
	}
	*/
	
	// Calculate energies and pressures
	upot = uinter + uintra;
	utot = upot + ukin;
	virial = virial_inter + virial_intra;
	// Add energy of thermostat. 
	// If nose hoover is not used, they will just be zero
	if (what_simulation==MOLECULAR_DYNAMICS)
	{
		if (what_ensemble==NVT)
		{
			utot = utot + unhts + unhtss;
		}
		else if (what_ensemble==NPT)
		{
			utot = utot + utsbs;
		}
	}
	// Calculate ideal pressure part, rho*K*T = rho(*)*T(*).
	pideal = natom/boxv*tinst;
	// Do not need to recalculate lrc here. It should be calculated
	// elsewhere when variables changed.
	pinst = pideal + virial/3.0/boxv;
	// Add long range corrections into total energy and pressure if needed.
	if (isLJlrcOn)
	{
		utot += uljlrc;
		pinst += pljlrc;
	}
	fprintf(fpouts, "The initial energies etc.:\n");
	fprintf(fpouts, "utot=%lf\n", utot); // utot
	fprintf(fpouts, "upot=%lf\n", upot); // upot
	fprintf(fpouts, "ukin=%lf\n", ukin);
	fprintf(fpouts, "uinter=%lf\n", uinter);
	fprintf(fpouts, "uintra=%lf\n", uintra);
	fprintf(fpouts, "uvdw=%lf\n", uvdw);
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
	fprintf(fpouts, "udftmcff=%lf (%lf ev/atom)\n", udftmcff, udftmcff/EV_TO_K/natom);
	fprintf(fpouts, "usg=%lf\n", usg);
	fprintf(fpouts, "ushift=%lf\n", ushift);
	fprintf(fpouts, "tinst=%lf\n", tinst);
	fprintf(fpouts, "virial=%lf\n", virial);
	fprintf(fpouts, "virial_inter=%lf\n", virial_inter);
	fprintf(fpouts, "virial_intra=%lf\n", virial_intra);
	fprintf(fpouts, "pideal=%lf\n", pideal);
	// pljlrc already calculated in initialization part
	fprintf(fpouts, "pressure=%lf\n", pinst);
	fprintf(fpouts, "uljlrc=%lf\npljlrc=%lf\n", uljlrc, pljlrc);
	fflush(fpouts);
	fflush(stderr);

	return 0;
}
