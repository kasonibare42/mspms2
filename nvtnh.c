#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vars.h"

int nvtnh()
{
    int ii;
    double mvsq;

    mvsq = 2.0*ukin;
    // compute the driving force for the thermostat
    Gts = (mvsq - NRT)/Qts;

    // advance the thermostat velocity 1/4 time step
    vts = vts + dt_outer4*Gts;

    // advance the thermostat position 1/2 timestep
    rts = rts + dt_outer2*vts;

    // advance velocities of atoms @ 1/2 time step
    AA = exp(-dt_outer2*vts);

    for (ii=0;ii<natom;ii++)
    {
	vx[ii] *= AA;
	vy[ii] *= AA;
	vz[ii] *= AA;
    }

    // compute new driving force for the thermostat and advance velocities @ 1/4 time step
    mvsq = mvsq*AA*AA;
    Gts = (mvsq - NRT)/Qts;
    vts = vts + dt_outer4*Gts;

    ukin = 0.5*mvsq;

    ukin_nhts = 0.5*Gts*vts*vts;
    upot_nhts = treq*Rgas*rts*nfree;
}
