
#ifndef PHYSICS_LINEAR_LINEAR_H_
#define PHYSICS_LINEAR_LINEAR_H_

#include "../physics.h"

namespace physics {

template <MatrixComposer TMatrixComposer>
class Linear: public Physics {

protected:
	void K(SparseSolver &K);

};


}


#endif /* PHYSICS_LINEAR_LINEAR_H_ */
