
#ifndef SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_HEATTRANSFER_EXTERNALHEAT_H_
#define SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_HEATTRANSFER_EXTERNALHEAT_H_

#include "analysis/assembler/subkernel/boundarycondition.h"

namespace espreso {

struct ExternalHeat: BoundaryCondition {

	ExternalHeat()
	: area(0)
	{
		isconst = false;
		action = Assembler::ASSEMBLE | Assembler::REASSEMBLE;
	}

	void activate(double area, double *rhs)
	{
		this->area = area;
		this->rhs = rhs;
		this->isactive = true;
	}

	double area;
};

template <size_t nodes, size_t gps, size_t ndim, class Physics> struct ExternalHeatKernel;

template <size_t nodes, size_t gps, class Physics>
struct ExternalHeatKernel<nodes, gps, 2, Physics>: ExternalHeat, Physics {
	ExternalHeatKernel(const ExternalHeat &base): ExternalHeat(base) { area = 1 / base.area; }

	void simd(typename Physics::Element &element)
	{
		double * __restrict__ out = rhs;
		for (size_t n = 0; n < nodes; ++n) {
			SIMD heat = load(out + n * SIMD::size);
			for (size_t gp = 0; gp < gps; ++gp) {
				SIMD q;
				q = q + element.ecf.heatFlow[gp] * load1(area);
				q = q + element.ecf.heatFlux[gp];
				q = q + element.ecf.htc[gp] * element.ecf.extTemp[gp];

				heat = heat + q * element.ecf.thickness[gp] * element.det[gp] * load1(element.w[gp]) * load1(element.N[gp][n]);
			}
			store(out + n * SIMD::size, heat);
		}
		rhs += nodes * SIMD::size;
	}
};

template <size_t nodes, size_t gps, class Physics>
struct ExternalHeatKernel<nodes, gps, 3, Physics>: ExternalHeat, Physics {
	ExternalHeatKernel(const ExternalHeat &base): ExternalHeat(base) {}

	void simd(typename Physics::Element &element)
	{
		double * __restrict__ out = rhs;
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
		rhs += nodes * SIMD::size;
	}
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_HEATTRANSFER_EXTERNALHEAT_H_ */
