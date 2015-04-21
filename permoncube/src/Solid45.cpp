
#include "Solid45.h"

CSolid45::CSolid45():
indSubdomain(0),
ordinalNumber(0)
{
	//int inod_glob[8];
 memset(inod_glob,0, 8 * sizeof(longInt));// ZeroMemm2(inod_glob, 8 * sizeof(int));
}

CSolid45::~CSolid45() {
	// TODO Auto-generated destructor stub
}

