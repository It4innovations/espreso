
#ifndef SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_HEATTRANSFER_HEATSOURCE_H_
#define SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_HEATTRANSFER_HEATSOURCE_H_

#include "analysis/assembler/subkernel/boundarycondition.h"

namespace espreso {

template <size_t nodes, size_t gps, size_t ndim, class Physics> struct HeatSourceKernel: BoundaryCondition, Physics {
	HeatSourceKernel(const BoundaryCondition &base): BoundaryCondition(base) {}

	void simd(typename Physics::Element &element)
	{
		printf("source\n");
		double * __restrict__ out = rhs;
		for (size_t n = 0; n < nodes; ++n) {
			SIMD heat = load(out + n * SIMD::size);
			for (size_t gp = 0; gp < gps; ++gp) {
				heat = heat + element.ecf.heatSource[gp] * element.det[gp] * load1(element.w[gp]) * load1(element.N[gp][n]);
			}
			store(out + n * SIMD::size, heat);
		}
		rhs += nodes * SIMD::size;
	}
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_HEATTRANSFER_HEATSOURCE_H_ */
