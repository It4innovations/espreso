
#include "ContactBC.h"

CContactBC::CContactBC(): ind(NULL), val(NULL), n(0)
{
}

CContactBC::~CContactBC() {
	if (ind) {delete [] ind; }
	if (val) {delete [] val; }
}

