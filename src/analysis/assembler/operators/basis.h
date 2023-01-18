
#ifndef SRC_ANALYSIS_ASSEMBLER_OPERATORS_BASIS_H_
#define SRC_ANALYSIS_ASSEMBLER_OPERATORS_BASIS_H_

#include "analysis/assembler/operator.h"
#include "analysis/assembler/parameter.h"
#include "analysis/assembler/math.hpp"

namespace espreso {

template <size_t nodes, size_t gps, size_t ndim, size_t edim, class Physics>
struct SetBaseFunctions: ActionOperator, Physics {
	InputParameterIterator w, N, dN;

	SetBaseFunctions(
			int interval,
			const ParameterData &w,
			const ParameterData &N,
			const ParameterData &dN)
	: w(w, interval),
	  N(N, interval, 0),
	  dN(dN, interval, 0)
	{

	}

	void sisd(typename Physics::Element &element)
	{
		for (size_t gp = 0; gp < gps; ++gp) {
			element.w[gp] = w[gp];
			for (size_t node = 0; node < nodes; ++node) {
				element.N[gps * node + gp] = N[gps * node + gp];
				for (size_t d = 0; d < edim; ++d) {
					element.dN[edim * gps * node + gp * edim + d] = dN[edim * gps * node + gp * edim + d];
				}
			}
		}
	}

	void simd(typename Physics::Element &element)
	{
		for (size_t s = 0; s < SIMD::size; ++s) {
			for (size_t gp = 0; gp < gps; ++gp) {
				element.w[SIMD::size * gp + s] = w[gp];
				for (size_t node = 0; node < nodes; ++node) {
					element.N[SIMD::size * (gps * node + gp) + s] = N[gps * node + gp];
					for (size_t d = 0; d < edim; ++d) {
						element.dN[SIMD::size * (edim * gps * node + gp * edim + d) + s] = dN[edim * gps * node + gp * edim + d];
					}
				}
			}
		}
	}
};


}

#endif /* SRC_ANALYSIS_ASSEMBLER_OPERATORS_BASIS_H_ */
