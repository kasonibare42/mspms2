#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include "vars.h"

/// Echo the input data and initialized variables
int echo()
{
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
	fprintf(fpouts, "%s\n", sysname);
	fprintf(fpouts, "natom=%d\n", natom);
	fprintf(fpouts, "nmole = %d\n", nmole);
	fprintf(fpouts, "nconstraint=%d\n", nconstraint);

	fprintf(fpouts, "%d atoms.\n", natom);
	fprintf(fpouts, "%d bonds.\n", nbond);
	fprintf(fpouts, "%d angles.\n", nangle);
	fprintf(fpouts, "%d dihedrals.\n", ndih);
	fprintf(fpouts, "%d impropers.\n", nimp);
	fprintf(fpouts, "%d nonbonded pairs.\n", nnbp);
	fprintf(fpouts, "%lf kappa\n", kappa);
	fprintf(fpouts, "%d %d %d %d KMAX etc.\n", KMAXX, KMAXY, KMAXZ, KSQMAX);
	fprintf(fpouts, "preq=%lf\nQts=%lf\nQbs=%lf\n", preq, Qts, Qbs);
	// nvt
	// charge on
	// sf on
	// LJ switch on
	// copy in.mspms here?

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
	fprintf(fpouts, "pideal=%lf\n", pideal=natom/(boxlx*boxly*boxlz)*tinst
			*kb_1e30);
	// pljlrc already calculated in initialization part
	pinst = pideal + (virial_inter+virial_intra)*virial_to_pressure/(boxlx
			*boxly*boxlz) + pljlrc;
	fprintf(fpouts, "pressure=%lf\n", pinst);

	fprintf(fpouts, "uljlrc=%lf\npljlrc=%lf\n", uljlrc, pljlrc);

	fflush(fpouts);

	return 0;
}
