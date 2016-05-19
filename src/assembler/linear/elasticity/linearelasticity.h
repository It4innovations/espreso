#ifndef ASSEMBLER_LINEAR_ELASTICITY_LINEARELASTICITY_H_
#define ASSEMBLER_LINEAR_ELASTICITY_LINEARELASTICITY_H_

#include "../linear.h"

namespace espreso {

template <class TInput>
class LinearElasticity: public Linear<TInput> {

public:
	LinearElasticity(TInput &input): Linear<TInput>(input) {};

	virtual ~LinearElasticity() {};

protected:
	virtual void inertia(std::vector<double> &inertia, const Material &material);
	virtual void C(DenseMatrix &C, const Material &material);
	virtual double CP() { return 1; }
};

}

#include "linearelasticity.hpp"

#endif /* ASSEMBLER_LINEAR_ELASTICITY_LINEARELASTICITY_H_ */
