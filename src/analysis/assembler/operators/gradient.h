
#ifndef SRC_ANALYSIS_ASSEMBLER_OPERATORS_GRADIENT_H_
#define SRC_ANALYSIS_ASSEMBLER_OPERATORS_GRADIENT_H_

#include "analysis/assembler/operator.h"
#include "analysis/assembler/parameter.h"
#include "analysis/assembler/math.hpp"
#include "esinfo/meshinfo.h"
#include "mesh/store/elementstore.h"

namespace espreso {

struct OutputGradient: public ActionOperator {
	OutputGradient(int interval, const ParameterData &dND, const ParameterData &temperature, NamedData *gradient)
	: dND(dND, interval),
	  temp(temperature, interval),
	  gradient(gradient->data.data() + info::mesh->dimension * info::mesh->elements->eintervals[interval].begin)
	{

	}

	InputParameterIterator dND, temp;
	double* gradient;

	void operator++()
	{
		++dND; ++temp;
		gradient += info::mesh->dimension;
	}

	void move(int n)
	{
		dND += n; temp += n;
		gradient += n * info::mesh->dimension;
	}
};

template<size_t nodes, size_t gps>
struct OutputGradient2D: public OutputGradient {
	using OutputGradient::OutputGradient;

	void operator()()
	{
		gradient[0] = gradient[1] = 0;
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			ADDM2NMN1<nodes>(1. / gps, dND.data + 2 * nodes * gpindex, temp.data, gradient);
		}
	}
};

template<size_t nodes, size_t gps>
struct OutputGradient3D: public OutputGradient {
	using OutputGradient::OutputGradient;

	void operator()()
	{
		gradient[0] = gradient[1] = gradient[2] = 0;
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			ADDM3NMN1<nodes>(1. / gps, dND.data + 3 * nodes * gpindex, temp.data, gradient);
		}
	}
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_OPERATORS_GRADIENT_H_ */
