
#ifndef SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_SIGMA_H_
#define SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_SIGMA_H_

#include "analysis/assembler/general/subkernel.h"
#include "config/ecf/physics/structuralmechanics.h"

namespace espreso {

struct Sigma: SubKernel {
    const char* name() const { return "Sigma"; }

    Sigma()
    : behaviour(StructuralMechanicsGlobalSettings::ELEMENT_BEHAVIOUR::PLANE_STRAIN)
    {
        isconst = false;
        action = SubKernel::SOLUTION;
    }

    void activate(StructuralMechanicsGlobalSettings::ELEMENT_BEHAVIOUR behaviour)
    {
        this->behaviour = behaviour;
        this->isactive = 1;
    }

    StructuralMechanicsGlobalSettings::ELEMENT_BEHAVIOUR behaviour;
};

template <size_t nodes, size_t ndim> struct SigmaKernel;

template <size_t nodes>
struct SigmaKernel<nodes, 2>: Sigma {
    SigmaKernel(const Sigma &base): Sigma(base) {}

    template <typename Element>
    void reset(Element &element)
    {

    }

    template <typename Element>
    void simd(Element &element, size_t gp)
    {

    }
};

template <size_t nodes>
struct SigmaKernel<nodes, 3>: Sigma {
    SigmaKernel(const Sigma &base): Sigma(base) {}

    template <typename Element>
    void reset(Element &element)
    {
        element.sigma[0] = element.sigma[1] = element.sigma[2] = element.sigma[3] = element.sigma[4] = element.sigma[5] = zeros();
    }

    template <typename Element>
    void simd(Element &element, size_t gp)
    {
        multAB<6, 6, 1>(element.sigma, element.vC4, element.smallStrainTensor);
    }
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_SIGMA_H_ */
