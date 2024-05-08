
#ifndef SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_SMALLSTRAINTENSOR_H_
#define SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_SMALLSTRAINTENSOR_H_

#include "element.h"
#include "analysis/assembler/general/subkernel.h"

namespace espreso {

struct SmallStrainTensor: SubKernel {
    const char* name() const { return "SmallStrainTensor"; }

    SmallStrainTensor()
    {
        isconst = false;
        action = SubKernel::ITERATION | SubKernel::SOLUTION;
    }

    void activate()
    {
        this->isactive = 1;
    }
};

template <size_t nodes, size_t ndim> struct SmallStrainTensorKernel;

template <size_t nodes>
struct SmallStrainTensorKernel<nodes, 2>: SmallStrainTensor {
    SmallStrainTensorKernel(const SmallStrainTensor &base): SmallStrainTensor(base) {}

        template <typename Element>
        void simd(Element &element, size_t gp)
        {
            element.smallStrainTensor[0] = zeros();
            element.smallStrainTensor[1] = zeros();
            element.smallStrainTensor[2] = zeros();
            for (size_t n = 0; n < nodes; ++n) {
                element.smallStrainTensor[0] = element.smallStrainTensor[0] + element.dND[n][0] * element.displacement[n][0];
                element.smallStrainTensor[1] = element.smallStrainTensor[1] + element.dND[n][1] * element.displacement[n][1];

                element.smallStrainTensor[2] = element.smallStrainTensor[2] + element.dND[n][1] * element.displacement[n][0];
                element.smallStrainTensor[2] = element.smallStrainTensor[2] + element.dND[n][0] * element.displacement[n][1];
            }
        }
};

template <size_t nodes>
struct SmallStrainTensorKernel<nodes, 3>: SmallStrainTensor {
    SmallStrainTensorKernel(const SmallStrainTensor &base): SmallStrainTensor(base) {}

    template <typename Element>
    void simd(Element &element, size_t gp)
    {
        element.smallStrainTensor[0] = zeros();
        element.smallStrainTensor[1] = zeros();
        element.smallStrainTensor[2] = zeros();
        element.smallStrainTensor[3] = zeros();
        element.smallStrainTensor[4] = zeros();
        element.smallStrainTensor[5] = zeros();
        for (size_t n = 0; n < nodes; ++n) {
            element.smallStrainTensor[0] = element.smallStrainTensor[0] + element.dND[n][0] * element.displacement[n][0];
            element.smallStrainTensor[1] = element.smallStrainTensor[1] + element.dND[n][1] * element.displacement[n][1];
            element.smallStrainTensor[2] = element.smallStrainTensor[2] + element.dND[n][2] * element.displacement[n][2];

            element.smallStrainTensor[3] = element.smallStrainTensor[3] + element.dND[n][1] * element.displacement[n][0];
            element.smallStrainTensor[3] = element.smallStrainTensor[3] + element.dND[n][0] * element.displacement[n][1];

            element.smallStrainTensor[4] = element.smallStrainTensor[4] + element.dND[n][2] * element.displacement[n][1];
            element.smallStrainTensor[4] = element.smallStrainTensor[4] + element.dND[n][1] * element.displacement[n][2];

            element.smallStrainTensor[5] = element.smallStrainTensor[5] + element.dND[n][2] * element.displacement[n][0];
            element.smallStrainTensor[5] = element.smallStrainTensor[5] + element.dND[n][0] * element.displacement[n][2];
        }
    }
};



}

#endif /* SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_SMALLSTRAINTENSOR_H_ */
