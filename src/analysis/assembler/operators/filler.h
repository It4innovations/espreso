
#ifndef SRC_PHYSICS_ASSEMBLER_OPERATORS_FILLER_H_
#define SRC_PHYSICS_ASSEMBLER_OPERATORS_FILLER_H_

#include "analysis/assembler/operator.h"
#include "analysis/assembler/parameter.h"

namespace espreso {

struct MatrixFiller: public ActionOperator {
	MatrixFiller(int interval, const ParameterData &stiffness, double *data, const esint *position)
	: ActionOperator(interval, false, true),
	  stiffness(stiffness, interval),
	  data(data), position(position) {}

	InputParameterIterator stiffness;
	double *data;
	const esint *position;

	void operator++()
	{
		// increment by operator()
	}

	void reset()
	{

	}
};

template<size_t nodes, size_t gps, size_t dimension>
struct MatrixUpperFiller: public MatrixFiller {
	using MatrixFiller::MatrixFiller;

	void operator()()
	{
		for (size_t r = 0; r < nodes * dimension; ++r, stiffness.data += nodes * dimension) {
			for (size_t c = r; c < nodes * dimension; ++c) {
				data[*position++] += *(stiffness.data + c);
			}
		}
	}
};

template<size_t nodes, size_t gps, size_t dimension>
struct MatrixLowerFiller: public MatrixFiller {
	using MatrixFiller::MatrixFiller;

	void operator()()
	{
		for (size_t r = 0; r < nodes * dimension; ++r, stiffness.data += nodes * dimension) {
			for (size_t c = 0; c <= r; ++c) {
				data[*position++] += *(stiffness.data + c);
			}
		}
	}
};

template<size_t nodes, size_t gps, size_t dimension>
struct MatrixFullFiller: public MatrixFiller {
	using MatrixFiller::MatrixFiller;

	void operator()()
	{
		for (size_t r = 0; r < nodes * dimension; ++r) {
			for (size_t c = 0; c < nodes * dimension; ++c) {
				printf("%lu %lu : [ %d ] = %f\n", r, c, *position, *(stiffness.data));
				data[*position++] += *stiffness.data++;
			}
		}
		printf("====\n");
	}
};

}


#endif /* SRC_PHYSICS_ASSEMBLER_OPERATORS_FILLER_H_ */
