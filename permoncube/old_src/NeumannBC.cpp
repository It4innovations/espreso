/*
 * NeumannBC.cpp
 *
 *  Created on: Oct 17, 2013
 *      Author: mila
 */

#include "NeumannBC.h"

CNeumannBC::CNeumannBC():
ind(NULL),
val(NULL),
n(0)
{
}

CNeumannBC::~CNeumannBC() {
	if (ind) delete [] ind;
	if (val) delete [] val;
}

