
#ifndef SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_HEATTRANSFER_OP_HEATSOURCE_H_
#define SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_HEATTRANSFER_OP_HEATSOURCE_H_

#include "analysis/assembler/general/subkernel.h"
#include "analysis/assembler/general/boundarycondition.h"

namespace espreso {

template <size_t nodes> struct HeatSourceKernel: BoundaryCondition {
    HeatSourceKernel(const BoundaryCondition &base): BoundaryCondition(base) {}

    template <typename Element>
    void simd(Element &element, size_t gp)
    {
        for (size_t n = 0; n < nodes; ++n) {
            element.f[n] = element.f[n] + element.ecf.heatSource * element.det * load1(element.w[gp]) * load1(element.N[gp][n]);
        }
    }
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_HEATTRANSFER_OP_HEATSOURCE_H_ */
