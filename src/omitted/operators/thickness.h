
#ifndef SRC_ANALYSIS_ASSEMBLER_OPERATORS_THICKNESS_H_
#define SRC_ANALYSIS_ASSEMBLER_OPERATORS_THICKNESS_H_

#include <analysis/assembler/subkernel/operator.h>
#include "basis/containers/serializededata.h"

namespace espreso {

struct Thickness: ActionOperator {
	const char* name() const { return "Thickness"; }

	serializededata<esint, esint>::const_iterator procNodes;
	double *target;

	Thickness(size_t interval, serializededata<esint, esint>::const_iterator procNodes, double *target)
	: procNodes(procNodes), target(target)
	{
		isconst = false;
		action = Action::ASSEMBLE;
	}

	Thickness(size_t region, size_t interval, serializededata<esint, esint>::const_iterator procNodes, double *target)
	: procNodes(procNodes), target(target)
	{
		isconst = false;
		action = Action::ASSEMBLE;
	}

	void move(int n)
	{
		procNodes += n;
	}
};

template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype, class Physics>
struct ThicknessToNodes: Thickness, Physics {
	using Thickness::Thickness;

	void simd(typename Physics::Element &element)
	{
		peel(element, SIMD::size);
	}

	void peel(typename Physics::Element &element, size_t size)
	{
		for (size_t s = 0; s < size; ++s, ++procNodes) {
			for (size_t n = 0; n < nodes; ++n) {
				target[procNodes->at(n)] = element.ecf.thickness[n][s];
			}
		}
	}
};

template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype, class Physics>
struct ThicknessToElementNodes: Thickness, Physics {
	using Thickness::Thickness;

	void simd(typename Physics::Element &element)
	{
		peel(element, SIMD::size);
	}

	void peel(typename Physics::Element &element, size_t size)
	{
		for (size_t s = 0; s < size; ++s, ++procNodes) {
			for (size_t n = 0; n < nodes; ++n) {
				element.thickness[n][s] = target[procNodes->at(n)];
			}
		}
	}
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_OPERATORS_THICKNESS_H_ */
