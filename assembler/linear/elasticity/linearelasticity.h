#ifndef ASSEMBLER_LINEAR_ELASTICITY_LINEARELASTICITY_H_
#define ASSEMBLER_LINEAR_ELASTICITY_LINEARELASTICITY_H_

#include "../linear.h"

namespace espreso {

template <class TInput>
class LinearElasticity: public Linear<TInput> {

public:
	LinearElasticity(TInput &input): Linear<TInput>(input) {};

protected:
	void inertia(std::vector<double> &inertia);
	void C(DenseMatrix &C, eslocal material);
	double CP() { return 1; }
	double rho() { return 7.85e-9; }
};

}

#include "linearelasticity.hpp"

#endif /* ASSEMBLER_LINEAR_ELASTICITY_LINEARELASTICITY_H_ */
