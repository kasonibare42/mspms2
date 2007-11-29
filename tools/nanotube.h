#ifndef _NANOTUBE_H
typedef struct _tubeatom 
{
    double xx, yy, zz;
} TUBEATOM, *PTUBEATOM;

typedef class _simnanotube
{
    public:
	int		id;
	int		type; // nanotube type, from 0 to 19
	int		ntubeatom;
	double		length;
	double		radius;
	double		weight;
	double		xx, yy; // coordinates of the z axis for smooth tubes
	double		sigma, epsilon, charge; // parameters for carbon atoms in the tube
	PTUBEATOM	tubeatom_list;

	_simnanotube();
	~_simnanotube();

} SIMNANOTUBE, *PSIMNANOTUBE;

typedef class _simgroove
{
    public:
	int		id;
	double		xx, yy;
	double		radius;
} SIMGROOVE, *PSIMGROOVE;

#define _NANOTUBE_H
#endif
