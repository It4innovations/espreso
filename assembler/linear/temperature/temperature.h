
#ifndef ASSEMBLER_LINEAR_TEMPERATURE_TEMPERATURE_H_
#define ASSEMBLER_LINEAR_TEMPERATURE_TEMPERATURE_H_

#include "../linear.h"

namespace assembler {

template <MatrixComposer TMatrixComposer>
class Temperature: public Linear<TMatrixComposer> {

public:
	Temperature(const mesh::Mesh &mesh): Linear<TMatrixComposer>(mesh) {};

protected:
	size_t DOFs() { return 1; }
	void inertia(std::vector<double> &inertia)
	{
		inertia.resize(1, 0);
	}
	void C(DenseMatrix &C);
	double CP() { return 1; }
	double rho() { return 7.85e-9; }
};

}

#include "temperature.hpp"




#endif /* ASSEMBLER_LINEAR_TEMPERATURE_TEMPERATURE_H_ */
