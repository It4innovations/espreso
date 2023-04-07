
#ifndef SRC_ANALYSIS_ASSEMBLER_OPERATORS_TEMPERATURE_H_
#define SRC_ANALYSIS_ASSEMBLER_OPERATORS_TEMPERATURE_H_

#include "analysis/assembler/operator.h"
#include "basis/containers/serializededata.h"

namespace espreso {

struct Temperature: ActionOperator {
	const char* name() const { return "Temperature"; }

	serializededata<esint, esint>::const_iterator procNodes;
	const double * const source;

	Temperature(size_t interval, serializededata<esint, esint>::const_iterator procNodes, const double * const source)
	: procNodes(procNodes), source(source)
	{
		isconst = false;
		action = Action::ASSEMBLE | Action::REASSEMBLE | Action::SOLUTION;
	}

	void move(int n)
	{
		procNodes += n;
	}
};

template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype, class Physics>
struct TemperatureToElementNodes: Temperature, Physics {
	using Temperature::Temperature;

	void simd(typename Physics::Element &element)
	{
		for (size_t s = 0; s < SIMD::size; ++s, ++procNodes) {
			for (size_t n = 0; n < nodes; ++n) {
				element.temp[n][s] = source[procNodes->at(n)];
			}
		}
	}

	void peel(typename Physics::Element &element, size_t size)
	{
		for (size_t s = 0; s < size; ++s, ++procNodes) {
			for (size_t n = 0; n < nodes; ++n) {
				element.temp[n][s] = source[procNodes->at(n)];
			}
		}
	}
};

template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype, class Physics>
struct TemperatureToElementNodesAndGPs: Temperature, Physics {
	using Temperature::Temperature;

	void simd(typename Physics::Element &element)
	{
		for (size_t s = 0; s < SIMD::size; ++s, ++procNodes) {
			for (size_t n = 0; n < nodes; ++n) {
				element.temp[n][s] = source[procNodes->at(n)];
			}
		}
		for (size_t gp = 0; gp < gps; ++gp) {
			element.gptemp[gp] = zeros();
			for (size_t n = 0; n < nodes; ++n) {
				element.gptemp[gp] = element.gptemp[gp] + load1(element.N[gp][n]) * element.temp[n];
			}
		}
	}

	void peel(typename Physics::Element &element, size_t size)
	{
		for (size_t s = 0; s < size; ++s, ++procNodes) {
			for (size_t n = 0; n < nodes; ++n) {
				element.temp[n][s] = source[procNodes->at(n)];
			}
		}
		for (size_t gp = 0; gp < gps; ++gp) {
			element.gptemp[gp] = zeros();
			for (size_t n = 0; n < nodes; ++n) {
				element.gptemp[gp] = element.gptemp[gp] + load1(element.N[gp][n]) * element.temp[n];
			}
		}
	}
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_OPERATORS_TEMPERATURE_H_ */
