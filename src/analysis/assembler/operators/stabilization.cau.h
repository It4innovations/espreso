
#ifndef SRC_ANALYSIS_ASSEMBLER_OPERATORS_STABILIZATION_CAU_H_
#define SRC_ANALYSIS_ASSEMBLER_OPERATORS_STABILIZATION_CAU_H_

#include "analysis/assembler/operator.h"
#include "math/simd/simd.h"

namespace espreso {

struct Stabilization: ActionOperator {
	const char* name() const { return "Stabilization"; }

	Stabilization(size_t interval)
	{
		isconst = false;
		action = Action::ASSEMBLE | Action::REASSEMBLE;
	}

	void move(int n)
	{

	}
};

template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype, class Physics> struct StabilizationCAU;

template <size_t nodes, size_t gps, size_t edim, size_t etype, class Physics>
struct StabilizationCAU<nodes, gps, 2, edim, etype, Physics>: Stabilization, Physics {
	using Stabilization::Stabilization;

	void simd(typename Physics::Element &element)
	{

	}
};

template <size_t nodes, size_t gps, size_t edim, size_t etype, class Physics>
struct StabilizationCAU<nodes, gps, 3, edim, etype, Physics>: Stabilization, Physics {
	using Stabilization::Stabilization;

	void simd(typename Physics::Element &element)
	{

	}
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_OPERATORS_STABILIZATION_CAU_H_ */
