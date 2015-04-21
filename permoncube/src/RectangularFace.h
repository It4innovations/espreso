
#ifndef RECTANGULARFACE_H_
#define RECTANGULARFACE_H_

#include "utility.h"

class CRectangularFace {
public:
	  longInt inod_glob[4];
	  int iFaceSub;
	  int iFace;
	  int iElem;

public:
	CRectangularFace();
	virtual ~CRectangularFace();
};

#endif /* RECTANGULARFACE_H_ */
