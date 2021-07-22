
#ifndef SRC_PHYSICS_ASSEMBLER_OPERATORS_FILLER_H_
#define SRC_PHYSICS_ASSEMBLER_OPERATORS_FILLER_H_

#include "analysis/assembler/operator.h"
#include "analysis/assembler/parameter.h"

namespace espreso {

struct MatrixFiller: public ActionOperator {
	MatrixFiller(int interval, const ParameterData &local, double *global, const esint *position)
	: ActionOperator(interval, false, true),
	  local(local, interval),
	  global(global), position(position) {}

	InputParameterIterator local;
	double *global;
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
		for (size_t r = 0; r < nodes * dimension; ++r, local.data += nodes * dimension) {
			for (size_t c = r; c < nodes * dimension; ++c) {
				global[*position++] += *(local.data + c);
			}
		}
	}
};

template<size_t nodes, size_t gps, size_t dimension>
struct MatrixLowerFiller: public MatrixFiller {
	using MatrixFiller::MatrixFiller;

	void operator()()
	{
		for (size_t r = 0; r < nodes * dimension; ++r, local.data += nodes * dimension) {
			for (size_t c = 0; c <= r; ++c) {
				global[*position++] += *(local.data + c);
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
				global[*position++] += *local.data++;
			}
		}
	}
};

template<size_t nodes, size_t gps, size_t dimension>
struct VectorFiller: public ActionOperator {
	VectorFiller(int interval, const ParameterData &rhs, double *global, const esint *position)
	: ActionOperator(interval, false, true),
	  rhs(rhs, interval),
	  global(global), position(position) {}

	InputParameterIterator rhs;
	double *global;
	const esint *position;

	void operator++()
	{
		// increment by operator()
	}

	void operator()()
	{
		for (size_t r = 0; r < nodes * dimension; ++r) {
			global[*position++] += *rhs.data++;
		}
	}

	void reset()
	{

	}
};

}


#endif /* SRC_PHYSICS_ASSEMBLER_OPERATORS_FILLER_H_ */
