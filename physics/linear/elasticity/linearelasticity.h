#ifndef LINEARELASTICITY_H_
#define LINEARELASTICITY_H_

#include "../linear.h"

namespace physics {

template <MatrixComposer TMatrixComposer>
class LinearElasticity: public Linear<TMatrixComposer> {

public:
	LinearElasticity(): Linear<TMatrixComposer>() {};

};

}

#endif /* LINEARELASTICITY_H_ */
