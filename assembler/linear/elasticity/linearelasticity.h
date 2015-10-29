#ifndef ASSEMBLER_LINEAR_ELASTICITY_LINEARELASTICITY_H_
#define ASSEMBLER_LINEAR_ELASTICITY_LINEARELASTICITY_H_

#include "../linear.h"

namespace assembler {

template <MatrixComposer TMatrixComposer>
class LinearElasticity: public Linear<TMatrixComposer> {

public:
	LinearElasticity(const mesh::Mesh &mesh): Linear<TMatrixComposer>(mesh) {};

protected:
	size_t DOFs() { return 3; }
	void inertia(std::vector<double> &inertia)
	{
		inertia.resize(3, 0);
		inertia[2] = 9810.0 * 7.85e-9;
	}
	void C(DenseMatrix &C);
	double CP() { return 1; }
	double rho() { return 7.85e-9; }

};

}

#include "linearelasticity.hpp"

#endif /* ASSEMBLER_LINEAR_ELASTICITY_LINEARELASTICITY_H_ */
