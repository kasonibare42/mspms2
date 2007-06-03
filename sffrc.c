#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "vars.h"

int sffrc()
{
    int ii;
    double usflj_tasos;
    double tasos_force[3];

    usflj = 0.0;

    if (sf_type==nanotube_polynomial)
    {
    }
    else if (sf_type==nanotube_atom_explicit)
    {}
    else if (sf_type==nanotube_tasos)
    {
	for (ii=0;ii<natom;ii++)
	{
	    // call Tasos's code for force calculations
	    cforce_atom_(&tasostype[ii],&xx[ii],&yy[ii],&zz[ii],&usflj_tasos,tasos_force);

	    // the tasos energy has unit of K and the force has unit of K/Angstrom
	    // change them to J/mol and J/mol/Angstrom
	    usflj_tasos *= Rgas;
	    tasos_force[0] *= Rgas;
	    tasos_force[1] *= Rgas;
	    tasos_force[2] *= Rgas;

	    usflj += usflj_tasos;
	    fxl[ii] += tasos_force[0];
	    fyl[ii] += tasos_force[1];
	    fzl[ii] += tasos_force[2];
	} // natom loop
    }
}

