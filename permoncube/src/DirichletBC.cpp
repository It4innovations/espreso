/*
 * DirichletBC.cpp
 *
 *  Created on: Oct 17, 2013
 *      Author: mila
 */

#include "DirichletBC.h"

CDirichletBC::CDirichletBC(): ind(NULL), val(NULL), n(0)
{
}

CDirichletBC::~CDirichletBC() {
	if (ind) delete [] ind;
	if (val) delete [] val;
}

