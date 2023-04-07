
#ifndef SRC_ANALYSIS_ASSEMBLER_OPERATORS_INITIALTEMPERATURE_H_
#define SRC_ANALYSIS_ASSEMBLER_OPERATORS_INITIALTEMPERATURE_H_

#include "analysis/assembler/operator.h"
#include "basis/containers/serializededata.h"

namespace espreso {

struct InitialTemperature: ActionOperator {
	const char* name() const { return "InitialTemperature"; }

	serializededata<esint, esint>::const_iterator procNodes;
	double *target;

	InitialTemperature(size_t interval, serializededata<esint, esint>::const_iterator procNodes, double *target)
	: procNodes(procNodes), target(target)
	{
		isconst = false;
		action = Action::ASSEMBLE | Action::REASSEMBLE | Action::SOLUTION;
	}

	InitialTemperature(size_t region, size_t interval, serializededata<esint, esint>::const_iterator procNodes, double *target)
	: procNodes(procNodes), target(target)
	{
		isconst = false;
		action = Action::ASSEMBLE | Action::REASSEMBLE;
	}

	void move(int n)
	{
		procNodes += n;
	}
};

template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype, class Physics>
struct InitialTemperatureToNodes: InitialTemperature, Physics {
	using InitialTemperature::InitialTemperature;

	void simd(typename Physics::Element &element)
	{
		peel(element, SIMD::size);
	}

	void peel(typename Physics::Element &element, size_t size)
	{
		for (size_t s = 0; s < size; ++s, ++procNodes) {
			for (size_t n = 0; n < nodes; ++n) {
				target[procNodes->at(n)] = element.temp[n][s];
			}
		}
	}
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_OPERATORS_INITIALTEMPERATURE_H_ */
