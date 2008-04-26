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

	fprintf(fpouts, "=================================================================================\n");
	fprintf(fpouts, "Energies are for all molecules.\n");
	fprintf(fpouts, "Only total energy (utot) contains long range corrections.\n");
	fprintf(fpouts, "Pressure is with long range corrections.\n");
	fprintf(fpouts, "=================================================================================\n");
	fprintf(fpouts, "This program uses reduced units:\n");
	fprintf(fpouts, "  Energy(%lf K)\n", epsilon_base);  
	fprintf(fpouts, "  Temperature(%lf K)\n", epsilon_base);
	fprintf(fpouts, "  Separation(%lf Angstrom)\n", sigma_base);
	fprintf(fpouts, "  Mass (%lf*E-27 kg/molecule)\n", mass_base);
	fprintf(fpouts, "  Pressure(%lf bar)\n", pressure_base);
	fprintf(fpouts, "  Time (%lf fs)\n", time_base);
	fprintf(fpouts, "  Virial(%lf K)\n", epsilon_base);
	fprintf(fpouts, "  Force (%lf K/Angstrom)\n", force_base);
	fprintf(fpouts, "  Velocity (%lf Angstrom/fs)\n", velocity_base);
	fprintf(fpouts, "=================================================================================\n");
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
		break;
	}
	fprintf(fpouts, "Running %s simulation (%d) with %s ensemble (%d).\n",
			szSimulation, what_simulation, szEnsemble, what_ensemble);
	fprintf(fpouts, "Length: %lf, Width: %lf, Height: %lf. Volume: %lf.\n", boxlx, boxly, boxlz, boxv);
	fprintf(fpouts, "nspecie=%d, nmole=%d, natom=%d\n", nspecie, nmole, natom);
	fprintf(fpouts, "nbond=%d, nangle=%d, ndih=%d, nimp=%d, nnbp=%d\n",
			nbond, nangle, ndih, nimp, nnbp);
	fprintf(fpouts, "nconstraints=%d.\n", nconstraint);
	switch (iStart_option)
	{
	case NEW:
		strcpy(szOutput, "New run");
		break;
	case CONTINUE:
		strcpy(szOutput, "Continue previoius runs");
		break;
	case CONFIG_ONLY:
		strcpy(szOutput, "New run from old configuration");
		break;
	default:
		strcpy(szOutput, "Unknown");
		break;
	}
	fprintf(fpouts, "iStart_option=%d (%s)\n", iStart_option, szOutput);
	fprintf(fpouts, "nstep=%d\n", nstep);
	fprintf(fpouts, "nstep_eq=%d\n", nstep_eq);
	fprintf(fpouts, "nstep_ave=%d\n", nstep_ave);
	fprintf(fpouts, "nstep_print=%d\n", nstep_print);
	fprintf(fpouts, "nstep_save=%d\n", nstep_save);
	fprintf(fpouts, "nstep_ss=%d\n", nstep_ss);
	fprintf(fpouts, "nstep_trj=%d\n", nstep_trj);
	fprintf(fpouts, "Random seeds: %d (ij), %d (jk)\n", ij, jk);
	fprintf(fpouts, "Required T(*): %lf.\n", treq);
	fprintf(fpouts, "Required P(*): %lf.\n", preq);
	fprintf(fpouts, "delt=%lf (%d inner steps) --> delts=%lf\n", delt, nstep_inner, delts);
	fprintf(fpouts, "rcutoff(LJ)=%lf; rcutoffelec=%lf.\n", rcutoff, rcutoffelec);
	
	fprintf(fpouts, "1,4 LJ interaction modifier (f0): %lf.\n", f0);
	if (f0==0.0)
	{
		fprintf(fpouts,
				"--->No 1,4 LJ interaction will be calculated, e.g. TraPPE.\n");
	}
	else if (f0==0.5)
	{
		fprintf(fpouts,
				"--->The 1,4 LJ interactions are reduced by half, e.g. OPLS potential.\n");
	}
	else
	{
		fprintf(fpouts, "--->The 1,4 LJ interactions are reduced.\n");
	}
	if (isLJswitchOn)
	{
		fprintf(
				fpouts, "Switch potential is enabled (%d) with Cuton=%lf.\n", isLJswitchOn, rcuton);
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
			fprintf(fpouts, "with %s boundary condition (%d).\n", szOutput, fEwald_BC);
			fprintf(fpouts, "    kappa=%lf, ", kappa);
			fprintf(fpouts, "kmaxx=%d, kmaxy=%d, kmaxz=%d, kmaxsq=%d\n", KMAXX, KMAXY, KMAXZ, KSQMAX);
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
	}
	else
	{
		fprintf(fpouts,
				"Electrostatic interaction calculations are disabled (%d).\n", iChargeType);
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
			fprintf(fpouts, "    Initialized variables: vts=%lf, vbs=%lf, rts=%lf, utsbs=%lf\n", vts, vbs, rts, utsbs);
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

	PSAMPLE_MOLECULE pSampleMole;
	int ii1, ii2, ii3, ii4;
	fprintf(fpouts, "=================================================================================\n");
	fprintf(fpouts, "Sample molecule configurations:\n");
	for (ii=0; ii<nspecie; ii++)
	{
		fprintf(fpouts, "---------------------------------------------------------------------------------\n");
		pSampleMole = sample_mole + ii;
		fprintf(fpouts, "Speice %d (%s): %d moles\n", ii, pSampleMole->mole_name, nmole_per_specie[ii]);
		fprintf(fpouts, "--->Sample mole has %d atoms.\n",pSampleMole->natom);
		for (jj=0;jj<pSampleMole->natom;jj++)
		{
			fprintf(fpouts, "%s  %lf, %lf, %lf; %lf, %lf, %lf; %d, %d\n", 
					pSampleMole->atom_name[jj], pSampleMole->xx[jj], pSampleMole->yy[jj], pSampleMole->zz[jj],
					pSampleMole->sigma[jj], pSampleMole->epsilon[jj], pSampleMole->charge[jj],
					pSampleMole->ghost_type[jj], pSampleMole->interp_type[jj]);
		}
		fprintf(fpouts, "--->Sample mole has %d bonds.\n", pSampleMole->nbond);
		for (jj=0;jj<pSampleMole->nbond;jj++)
		{ 
			ii1 = pSampleMole->bnd_idx[jj][0];
			ii2 = pSampleMole->bnd_idx[jj][1];
			fprintf(fpouts, "%d) %d-%d (%s-%s); %d; %lf, %lf, %lf\n",
					jj, ii1, ii2, pSampleMole->atom_name[ii1], pSampleMole->atom_name[ii2],
					pSampleMole->bnd_type[jj], pSampleMole->bnd_para_0[jj],
					pSampleMole->bnd_para_1[jj], pSampleMole->bnd_para_2[jj]);
		}
		fprintf(fpouts, "--->Sample mole has %d angles.\n", pSampleMole->nangle);
		for (jj=0;jj<pSampleMole->nangle;jj++)
		{
			ii1 = pSampleMole->agl_idx[jj][0];
			ii2 = pSampleMole->agl_idx[jj][1];
			ii3 = pSampleMole->agl_idx[jj][2];
			fprintf(fpouts, "%d) %d-%d-%d (%s-%s-%s); %d; %lf, %lf, %lf, %lf, %lf\n",
				jj, ii1, ii2, ii3, pSampleMole->atom_name[ii1], pSampleMole->atom_name[ii2],
				pSampleMole->atom_name[ii3], pSampleMole->agl_type[jj], pSampleMole->agl_para_1[jj], pSampleMole->agl_para_2[jj],
			       pSampleMole->agl_para_3[jj], pSampleMole->agl_para_4[jj], pSampleMole->agl_para_5[jj]);
		}
		fprintf(fpouts, "--->Sample mole has %d dihedrals.\n", pSampleMole->ndih);
		for (jj=0;jj<pSampleMole->ndih;jj++)
		{
		    ii1 = pSampleMole->dih_idx[jj][0];
		    ii2 = pSampleMole->dih_idx[jj][1];
		    ii3 = pSampleMole->dih_idx[jj][2];
		    ii4 = pSampleMole->dih_idx[jj][3];
		    fprintf(fpouts, "%d) %d-%d-%d-%d (%s-%s-%s-%s); %d; %lf, %lf, %lf, %lf\n",
			    jj, ii1, ii2, ii3, ii4, pSampleMole->atom_name[ii1], pSampleMole->atom_name[ii2],
			    pSampleMole->atom_name[ii3], pSampleMole->atom_name[ii4], pSampleMole->dih_type[jj], pSampleMole->dih_para_1[jj],
			    pSampleMole->dih_para_2[jj], pSampleMole->dih_para_3[jj], pSampleMole->dih_para_4[jj]);
		}
		fprintf(fpouts, "--->Sample mole has %d impropers.\n", pSampleMole->nimp);
		for (jj=0;jj<pSampleMole->nimp;jj++)
		{
		    ii1 = pSampleMole->imp_idx[jj][0];
		    ii2 = pSampleMole->imp_idx[jj][1];
		    ii3 = pSampleMole->imp_idx[jj][2];
		    ii4 = pSampleMole->imp_idx[jj][3];
		    fprintf(fpouts, "%d) %d-%d-%d-%d (%s-%s-%s-%s); %d; %lf, %lf\n",
			    jj, ii1, ii2, ii3, ii4, pSampleMole->atom_name[ii1], pSampleMole->atom_name[ii2],
			    pSampleMole->atom_name[ii3], pSampleMole->atom_name[ii4], pSampleMole->imp_type[jj], pSampleMole->imp_para_0[jj],
			    pSampleMole->imp_para_1[jj]);
		}
		fprintf(fpouts, "--->Sample mole has %d non-bonded pairs.\n", pSampleMole->nnbp);
		for (jj=0;jj<pSampleMole->nnbp;jj++)
		{
		    ii1 = pSampleMole->nbp_idx[jj][0];
		    ii2 = pSampleMole->nbp_idx[jj][1];
		    fprintf(fpouts, "%d) %d-%d (%s-%s)\n",
			    jj, ii1, ii2, pSampleMole->atom_name[ii1], pSampleMole->atom_name[ii2]);
		} 
	}
	fprintf(fpouts, "=================================================================================\n");
	fprintf(fpouts, "Initial soft forces and hard forces (unit is K/Angstrom):\n");
	fprintf(fpouts, "---------------------------------------------------------------------------------\n");
	for (ii=0;ii<natom;ii++)
	{
		fprintf(fpouts, "%d: (soft) %10.5lf, %10.5lf, %10.5lf, (hard) %12.5lf, %12.5lf, %12.5lf\n",
				ii, fxl[ii]*force_base, fyl[ii]*force_base, fzl[ii]*force_base, 
				fxs[ii]*force_base, fys[ii]*force_base, fzs[ii]*force_base);
	}
	fprintf(fpouts, "---------------------------------------------------------------------------------\n");	
	fprintf(fpouts, "Initial velocities (unit is Angstrom/fs):\n");
	fprintf(fpouts, "---------------------------------------------------------------------------------\n");
	for (ii=0;ii<natom;ii++)
	{
		fprintf(fpouts, "%d: vx=%10.7lf, vy=%10.7lf, vz=%10.7lf\n", 
				ii, vx[ii]*velocity_base, vy[ii]*velocity_base, vz[ii]*velocity_base);
	}
	fprintf(fpouts, "=================================================================================\n");
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
	fprintf(fpouts, "The initial energies and pressures:\n");
	fprintf(fpouts, "---------------------------------------------------------------------------------\n");	
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
	fprintf(fpouts, "pressure=%lf\n", pinst);
	fprintf(fpouts, "uljlrc=%lf\npljlrc=%lf\n", uljlrc, pljlrc);
	fflush(fpouts);
	fflush(stderr);
	
	return 0;
}
