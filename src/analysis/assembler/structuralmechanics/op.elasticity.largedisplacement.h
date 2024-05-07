
#ifndef SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_ELASTICITY_LARGEDISPLACEMENT_H_
#define SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_ELASTICITY_LARGEDISPLACEMENT_H_

#include "analysis/assembler/general/subkernel.h"

namespace espreso {

struct ElasticityLargeDisplacement: SubKernel {
    const char* name() const { return "ElasticityLargeDisplacement"; }

    ElasticityLargeDisplacement()
    {
        isconst = false;
        action = SubKernel::ASSEMBLE | SubKernel::REASSEMBLE;
    }

    void activate()
    {
        this->isactive = 1;
    }
};


template <size_t ndim> struct ElasticityLargeDisplacementKernel;

template <>
struct ElasticityLargeDisplacementKernel<2>: ElasticityLargeDisplacement {
    ElasticityLargeDisplacementKernel(const ElasticityLargeDisplacement &base): ElasticityLargeDisplacement(base) {}

    template <typename Element>
    void simd(Element &element, size_t gp)
    {

    }
};

template <>
struct ElasticityLargeDisplacementKernel<3>: ElasticityLargeDisplacement {
    ElasticityLargeDisplacementKernel(const ElasticityLargeDisplacement &base): ElasticityLargeDisplacement(base) {}

    template <typename Element>
    void simd(Element &element, size_t gp)
    {
        // 0 1 2  0 3 6
        // 3 4 5  1 4 7
        // 6 7 8  2 5 8

        // 0 3 4
        //   1 5
        //     6

        SIMD C05 = load1(0.5), C2 = load1(2);
        SIMD eVec0 = -C05 + C05 * (element.F[0] * element.F[0] + element.F[1] * element.F[3] + element.F[2] * element.F[6]);
        SIMD eVec1 = -C05 + C05 * (element.F[3] * element.F[1] + element.F[4] * element.F[4] + element.F[5] * element.F[7]);
        SIMD eVec2 = -C05 + C05 * (element.F[6] * element.F[2] + element.F[7] * element.F[5] + element.F[8] * element.F[8]);

        SIMD eVec3 =   C2 * C05 * (element.F[0] * element.F[3] + element.F[1] * element.F[4] + element.F[2] * element.F[5]);
        SIMD eVec4 =   C2 * C05 * (element.F[3] * element.F[6] + element.F[4] * element.F[7] + element.F[5] * element.F[8]);
        SIMD eVec5 =   C2 * C05 * (element.F[0] * element.F[6] + element.F[1] * element.F[7] + element.F[2] * element.F[8]);

        element.sVec[0] = element.elasticity[ 0] * eVec0 + element.elasticity[ 1] * eVec1 + element.elasticity[ 2] * eVec2 + element.elasticity[ 3] * eVec3 + element.elasticity[ 4] * eVec4 + element.elasticity[ 5] * eVec5;
        element.sVec[1] = element.elasticity[ 6] * eVec0 + element.elasticity[ 7] * eVec1 + element.elasticity[ 8] * eVec2 + element.elasticity[ 9] * eVec3 + element.elasticity[10] * eVec4 + element.elasticity[11] * eVec5;
        element.sVec[2] = element.elasticity[12] * eVec0 + element.elasticity[13] * eVec1 + element.elasticity[14] * eVec2 + element.elasticity[15] * eVec3 + element.elasticity[16] * eVec4 + element.elasticity[17] * eVec5;
        element.sVec[3] = element.elasticity[18] * eVec0 + element.elasticity[19] * eVec1 + element.elasticity[20] * eVec2 + element.elasticity[21] * eVec3 + element.elasticity[22] * eVec4 + element.elasticity[23] * eVec5;
        element.sVec[4] = element.elasticity[24] * eVec0 + element.elasticity[25] * eVec1 + element.elasticity[26] * eVec2 + element.elasticity[27] * eVec3 + element.elasticity[28] * eVec4 + element.elasticity[29] * eVec5;
        element.sVec[5] = element.elasticity[30] * eVec0 + element.elasticity[31] * eVec1 + element.elasticity[32] * eVec2 + element.elasticity[33] * eVec3 + element.elasticity[34] * eVec4 + element.elasticity[35] * eVec5;
    }
};

}



#endif /* SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_ELASTICITY_LARGEDISPLACEMENT_H_ */
