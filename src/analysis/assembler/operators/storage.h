
#ifndef SRC_ANALYSIS_ASSEMBLER_OPERATORS_STORAGE_H_
#define SRC_ANALYSIS_ASSEMBLER_OPERATORS_STORAGE_H_

#include "analysis/assembler/operator.h"
#include "math/simd/simd.h"

#include <vector>

namespace espreso {

template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype, class Physics, class Parameter>
struct StorageStore: ActionOperator, Physics {
	const char* name() const { return "StorageStore"; }

	StorageStore(size_t interval, ActionOperator *op, size_t elements, size_t size, const Parameter &parameter)
	: parameter(parameter),
	  size(size),
	  inc(op->isconst ? 0 : size / sizeof(double) / SIMD::size),
	  data(op->isconst ? size / sizeof(double) : elements * size / sizeof(double) / SIMD::size + size / sizeof(double)),
	  iterator(data.data())
	{
		isconst = op->isconst;
		action = Action::ASSEMBLE;
		op->action = Action::ASSEMBLE;
	}

	Parameter parameter;
	size_t size, inc;
	std::vector<double> data;
	double* iterator;

	void move(int n)
	{
		iterator += n * inc;
	}

	void simd(typename Physics::Element &element)
	{
		memcpy(iterator, parameter(element), size);
		move(SIMD::size);
	}
};

template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype, class Physics, class Parameter>
struct StorageLoad: ActionOperator, Physics {
	const char* name() const { return "StorageLoad"; }

	StorageLoad(size_t interval, StorageStore<nodes, gps, ndim, edim, etype, Physics, Parameter> *store, const Parameter &parameter)
	: parameter(parameter),
	  size(store->size),
	  inc(store->inc),
	  iterator(store->iterator)
	{
		isconst = store->isconst;
		action = ActionOperator::REASSEMBLE | Action::SOLUTION;
	}

	Parameter parameter;
	const size_t size, inc;
	const double* iterator;

	void move(int n)
	{
		iterator += n * inc;
	}

	void simd(typename Physics::Element &element)
	{
		memcpy(parameter(element), iterator, size);
		move(SIMD::size);
	}
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_OPERATORS_STORAGE_H_ */
