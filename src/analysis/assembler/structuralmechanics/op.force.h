
#ifndef SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_FORCE_H_
#define SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_FORCE_H_

#include "analysis/assembler/general/boundarycondition.h"

namespace espreso {

struct Force: public BoundaryCondition {
    Force()
    {

    }

    void activate(ECFExpressionVector *force)
    {
        BoundaryCondition::activate(force);
    }
};

template <size_t ndim> struct ForceKernel: public Force {
    ForceKernel(const Force &base): Force(base) {}

    template <typename Element>
    void simd(Element &element)
    {
        for (size_t d = 0; d < ndim; ++d) {
            element.f[d] = element.f[d] + element.ecf.force[d];
        }
    }
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_FORCE_H_ */
