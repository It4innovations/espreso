
#ifndef ASSEMBLER_LINEAR_TEMPERATURE_TEMPERATURE_H_
#define ASSEMBLER_LINEAR_TEMPERATURE_TEMPERATURE_H_

#include "../linear.h"

namespace espreso {

template <class TInput>
class Temperature: public Linear<TInput> {

public:
	Temperature(TInput &input): Linear<TInput>(input) {};

protected:
	void inertia(std::vector<double> &inertia, const Material &material);
	void C(DenseMatrix &C, const Material &material);
	double CP() { return 1; }
	double rho() { return 7.85e-9; }
};

}

#include "temperature.hpp"




#endif /* ASSEMBLER_LINEAR_TEMPERATURE_TEMPERATURE_H_ */
