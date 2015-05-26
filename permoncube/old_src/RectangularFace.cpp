
#include "RectangularFace.h"

CRectangularFace::CRectangularFace():
iFaceSub(0)
{
	//int inod_glob[4];
	memset(inod_glob, 0, 4 * sizeof(longInt));// ZeroMemm2(inod_glob, 4 * sizeof(int));
}

CRectangularFace::~CRectangularFace() {
	// TODO Auto-generated destructor stub
}

