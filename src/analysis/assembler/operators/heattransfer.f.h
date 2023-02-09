
#ifndef SRC_ANALYSIS_ASSEMBLER_OPERATORS_HEATTRANSFER_F_H_
#define SRC_ANALYSIS_ASSEMBLER_OPERATORS_HEATTRANSFER_F_H_

#include "analysis/assembler/operator.h"

namespace espreso {

template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype, class Physics>
struct HeatSource: ActionOperator, Physics {
	OutputParameterIterator rhs;

	HeatSource(size_t interval, ParameterData &rhs)
	: rhs(rhs, interval)
	{
		isconst = false;
		action = Action::ASSEMBLE;
	}

	void move(int n)
	{
		rhs += n;
	}

	void simd(typename Physics::Element &element)
	{
		double * __restrict__ out = rhs.data;
		for (size_t n = 0; n < nodes; ++n) {
			SIMD heat = load(out + n * SIMD::size);
			for (size_t gp = 0; gp < gps; ++gp) {
				heat = heat + element.ecf.heatSource[gp] * element.det[gp] * load1(element.w[gp]) * load1(element.N[gp][n]);
			}
			store(out + n * SIMD::size, heat);
		}
		move(SIMD::size);
	}
};

template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype, class Physics>
struct BoundaryHeat: ActionOperator, Physics {
	OutputParameterIterator rhs;
	const double area;

	BoundaryHeat(size_t region, size_t interval, ParameterData &rhs)
	: rhs(rhs, interval), area(1.0 / info::mesh->boundaryRegions[region]->area)
	{
		isconst = false;
		action = Action::ASSEMBLE;
	}

	void move(int n)
	{
		rhs += n;
	}

	void simd(typename Physics::Element &element)
	{
		double * __restrict__ out = rhs.data;
		for (size_t n = 0; n < nodes; ++n) {
			SIMD heat = load(out + n * SIMD::size);
			for (size_t gp = 0; gp < gps; ++gp) {
				SIMD q;
				q = q + element.ecf.heatFlow[gp] * load1(area);
				q = q + element.ecf.heatFlux[gp];
				q = q + element.ecf.htc[gp] * element.ecf.extTemp[gp];

				heat = heat + q * element.det[gp] * load1(element.w[gp]) * load1(element.N[gp][n]);
			}
			store(out + n * SIMD::size, heat);
		}
		move(SIMD::size);
	}
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_OPERATORS_HEATTRANSFER_F_H_ */
