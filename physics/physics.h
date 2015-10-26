
#ifndef PHYSICS_PHYSICS_H_
#define PHYSICS_PHYSICS_H_

#include "essolver.h"
#include "esmesh.h"
#include "esbem.h"

namespace physics {

enum MatrixComposer {
	FEM,
	BEM,
	ELMER
};

class Physics {

public:
	virtual void init() = 0;
	virtual void pre_solve_update() = 0;
	virtual void post_solve_update() = 0;
	virtual void solve() = 0;
	virtual void finalize() = 0;

	virtual ~Physics() {};

};

}



#endif /* PHYSICS_PHYSICS_H_ */
