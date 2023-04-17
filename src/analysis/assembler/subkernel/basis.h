
#ifndef SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_BASIS_H_
#define SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_BASIS_H_

#include "subkernels.h"
#include "analysis/assembler/operators/basis.h"

#include "mesh/element.h"

namespace espreso {

template <Element::CODE code, size_t nodes, size_t gps, size_t edim, class Physics>
struct BasisElementKernel: BasisKernel, Physics {
	BasisElementKernel(const BasisKernel &base): BasisKernel(base) { }

	void simd(typename Physics::Element &element)
	{
		for (size_t gp = 0; gp < gps; ++gp) {
			element.w[gp] = GaussPoints<code, nodes, gps, edim>::w[gp];
			for (size_t n = 0; n < nodes; ++n) {
				element.N[gp][n] = GaussPoints<code, nodes, gps, edim>::N[gp * nodes + n];
				for (size_t d = 0; d < edim; ++d) {
					element.dN[gp][n][d] = GaussPoints<code, nodes, gps, edim>::dN[gp * edim * nodes + d * nodes + n];
				}
			}
		}
	}
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_BASIS_H_ */
