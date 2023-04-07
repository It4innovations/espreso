
#ifndef SRC_ANALYSIS_ASSEMBLER_OPERATORS_DISPLACEMENT_H_
#define SRC_ANALYSIS_ASSEMBLER_OPERATORS_DISPLACEMENT_H_

#include "analysis/assembler/operator.h"
#include "basis/containers/serializededata.h"

namespace espreso {

struct Displacement: ActionOperator {
	const char* name() const { return "Displacement"; }

	serializededata<esint, esint>::const_iterator procNodes;
	const double * const source;

	Displacement(size_t interval, serializededata<esint, esint>::const_iterator procNodes, const double * const source)
	: procNodes(procNodes), source(source)
	{
		isconst = false;
		action = Action::SOLUTION;
	}

	void move(int n)
	{
		procNodes += n;
	}
};

template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype, class Physics>
struct DisplacementToElementNodes: Displacement, Physics {
	using Displacement::Displacement;

	void simd(typename Physics::Element &element)
	{
		for (size_t s = 0; s < SIMD::size; ++s, ++procNodes) {
			for (size_t n = 0; n < nodes; ++n) {
				for (size_t d = 0; d < ndim; ++d) {
					element.displacement[n][d][s] = source[ndim * procNodes->at(n) + d];
				}
			}
		}
	}

	void peel(typename Physics::Element &element, size_t size)
	{
		for (size_t s = 0; s < size; ++s, ++procNodes) {
			for (size_t n = 0; n < nodes; ++n) {
				for (size_t d = 0; d < ndim; ++d) {
					element.displacement[n][d][s] = source[ndim * procNodes->at(n) + d];
				}
			}
		}
	}
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_OPERATORS_DISPLACEMENT_H_ */
