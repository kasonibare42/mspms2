#include <stdlib.h>
#include "nanotube.h"

_simnanotube::_simnanotube()
{
    tubeatom_list = NULL;
}

_simnanotube::~_simnanotube()
{
    if (tubeatom_list)
    {
	delete []tubeatom_list;
	tubeatom_list = NULL;
    }
}
