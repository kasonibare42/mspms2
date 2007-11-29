/***************************************************************************
 *            box.h
 *
 *  Sun Jul 17 2005
 *  Created by: Yang Wang
 ****************************************************************************/
#ifndef _BOX_H
typedef class _simbox
{
    public:
	int		id;
	int		nspecie;
	int		nnanotube;
	int		ngroove;
	int		nmole;
	int		natom;
	double		rfree; // 1.0/(3*natom -3)
	double		length, width, height, volume;
	double		weight;
	double		rho_atom, rho_mole, rho_real;
	double 		temperature;
	double		delt; // time step
	double		delv; // volume move
	double		disp; // max displacement
	double		cut, cutsq, r_cut, r_cutsq, r_cuttri;
	double		b0; // constant for reaction field

	// following variables only for checking molecular positions on nanotubes
	int		nMole_inside_general;
	int		nMole_first_layer_general;
	int		nMole_groove_general;
	int		nMole_gas_phase;
	int		nMole_to_inside;
	int		nMole_inside_to_gas;
	int		nMole_noninside_to_gas;

	PSIMSPECIE	specie_list;
	PSIMNANOTUBE	nanotube_list;
	PSIMGROOVE	groove_list;
	// counter[0] total moves
	// counter[1] moves in eq or taking data run
	// counter[2] canonical moves
	// counter[3] accepted canonical moves
	// counter[4] volume change moves
	// counter[5] accepted volume change moves
	// counter[6] insertions
	// counter[7] accepted insertions
	// counter[8] deletions
	// counter[9] accepted deletions
	// counter[10] number of average cycles
	long		counter[_MAX_COUNTER_BOX];
	// 0 utot
	// 1 utotcs
	// 2 ukin
	// 3 upot
	// 4 upotcs
	// 5 ulj
	// 6 uljcs
	// 7 ucl
	// 8 usflj
	// 9 usfljcs
	// 10 ubond
	// 11 uangle
	// 12 udih
	// 13 unplj
	// 14 unpljcs
	// 15 unpcl
	// 16 pressure
	// 17 plj
	// 18 pcl
	// 19 psflj
	// 20 pbond
	// 21 pangle
	// 22 pdih
	// 23 pnplj
	// 24 pnpcl
	// 25 temperature
	// 26 nmole
	// 27 volume
	// 28 rho_mole
	// 29 pideal
	ACCUMU		accumu[_MAX_ACCUMU_BOX];

	// all energies are for every molecule
	double		utot, utotcs;
	double		ukin, upot, upotcs; 
	double		ulj, uljcs, uljlrc;
	double		ucl, uclcs, ucllrc; // long range for cl interaction
	double		usflj, usfljcs;
	double		ubond, uangle, udih;
	double		unplj, unpljcs;
	double		unpcl;

	double 		*uljlrc_term;
	double		*pljlrc_term;
	double		uljs_1, uljs_2; // shift energy for LJ interactions (including atom explicit sf LJ)
	double		usfljs[_MAX_NANOTUBE][_POLYNOMIAL_ORDER]; // shift energy for smooth nanotubes

	double		pressure, pideal; // total pressure and ideal part of pressure
	double		plj, pcl, psflj, pbond, pangle, pdih, pnplj, pnpcl;
	double		pljlrc;

	_simbox();
	~_simbox();
	PSIMMOLE get_mole(int index);


}SIMBOX, *PSIMBOX;
#define _BOX_H   
#endif /* BOX_H */
