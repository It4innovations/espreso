
#ifndef ASSEMBLER_LINEAR_TEMPERATURE_TEMPERATURE_H_
#define ASSEMBLER_LINEAR_TEMPERATURE_TEMPERATURE_H_

#include "../linear.h"

namespace espreso {

template <class TInput>
class Temperature: public Linear<TInput> {

public:
	Temperature(TInput &input): Linear<TInput>(input) {};

	virtual ~Temperature() {};

protected:
	virtual void inertia(std::vector<double> &inertia, const Material &material);
	virtual void C(DenseMatrix &C, const Material &material);
	virtual double CP() { return 1; }
};

}

#include "temperature.hpp"




#endif /* ASSEMBLER_LINEAR_TEMPERATURE_TEMPERATURE_H_ */
