
#ifndef SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_HEATTRANSFER_OP_EXTERNALHEAT_H_
#define SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_HEATTRANSFER_OP_EXTERNALHEAT_H_

#include "analysis/assembler/general/element.h"
#include "analysis/assembler/general/boundarycondition.h"

namespace espreso {

struct ExternalHeat: BoundaryCondition {

	ExternalHeat()
	: area(0)
	{
		isconst = false;
		action = SubKernel::ASSEMBLE | SubKernel::REASSEMBLE | SubKernel::ITERATION;
	}

	void activate()
	{
		this->isactive = true;
	}

	double area;
};

template <size_t nodes> struct ExternalHeatKernel;

template <size_t nodes>
struct ExternalHeatKernel: ExternalHeat {
	ExternalHeatKernel(const ExternalHeat &base): ExternalHeat(base) { area = 1 / base.area; }

	template <typename Element>
	void simd(Element &element, size_t gp)
	{
		for (size_t n = 0; n < nodes; ++n) {
			SIMD q;
			q = q + element.ecf.heatFlow * load1(area);
			q = q + element.ecf.heatFlux;
			q = q + element.ecf.htc * element.ecf.extTemp;

			element.f[n] = element.f[n] + q * element.thickness.gp * element.det * load1(element.w[gp]) * load1(element.N[gp][n]);

		}
	}
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_HEATTRANSFER_OP_EXTERNALHEAT_H_ */
