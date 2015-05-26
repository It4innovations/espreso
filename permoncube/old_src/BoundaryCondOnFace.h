

#ifndef BOUNDARYCONDONFACE_H_
#define BOUNDARYCONDONFACE_H_

#include "utility.h"

class CBoundaryCondOnFace{
public:
	  int iFace;
	  double val_xyz[3];
public:
	  CBoundaryCondOnFace();
	virtual ~CBoundaryCondOnFace();
};

#endif /* BOUNDARYCONDONFACE_H_ */
