
#ifndef SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_INTEGRATION_DISPLACED_H_
#define SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_INTEGRATION_DISPLACED_H_

#include "analysis/assembler/general/subkernel.h"

namespace espreso {

struct IntegrationDisplaced: SubKernel {
    const char* name() const { return "IntegrationDisplaced"; }

    IntegrationDisplaced()
    {
        isconst = false;
        action = SubKernel::ASSEMBLE | SubKernel::REASSEMBLE;
    }

    void activate()
    {
        this->isactive = 1;
    }
};


template <size_t nodes, size_t ndim> struct IntegrationDisplacedKernel;

template <size_t nodes>
struct IntegrationDisplacedKernel<nodes, 2>: IntegrationDisplaced {
    IntegrationDisplacedKernel(const IntegrationDisplaced &base): IntegrationDisplaced(base) {}

    template <typename Element>
    void simd(Element &element, size_t gp)
    {
        SIMD JC0 = zeros(), JC1 = zeros(), JC2 = zeros(), JC3 = zeros();

        for (size_t n = 0; n < nodes; ++n) {
            SIMD coordsX = element.coords.node[n][0] + element.displacement[n][0];
            SIMD coordsY = element.coords.node[n][1] + element.displacement[n][1];
            SIMD dNX = load1(element.dN[gp][n][0]);
            SIMD dNY = load1(element.dN[gp][n][1]);

            JC0 = JC0 + dNX * coordsX;
            JC1 = JC1 + dNX * coordsY;
            JC2 = JC2 + dNY * coordsX;
            JC3 = JC3 + dNY * coordsY;
        }

        for (int r = 0; r < 2; ++r) {
            element.F[2 * r + 0] = element.invJ[0 * 2 + r] * JC0 + element.invJ[1 * 2 + r] * JC2;
            element.F[2 * r + 1] = element.invJ[0 * 2 + r] * JC1 + element.invJ[1 * 2 + r] * JC3;
        }
    }
};

template <size_t nodes>
struct IntegrationDisplacedKernel<nodes, 3>: IntegrationDisplaced {
    IntegrationDisplacedKernel(const IntegrationDisplaced &base): IntegrationDisplaced(base) {}

    template <typename Element>
    void simd(Element &element, size_t gp)
    {
        SIMD JC0 = zeros(), JC1 = zeros(), JC2 = zeros(), JC3 = zeros(), JC4 = zeros();
        SIMD JC5 = zeros(), JC6 = zeros(), JC7 = zeros(), JC8 = zeros();

        for (size_t n = 0; n < nodes; ++n) {
            SIMD coordsX = element.coords.node[n][0] + element.displacement[n][0];
            SIMD coordsY = element.coords.node[n][1] + element.displacement[n][1];
            SIMD coordsZ = element.coords.node[n][2] + element.displacement[n][2];
            SIMD dNX = load1(element.dN[gp][n][0]);
            SIMD dNY = load1(element.dN[gp][n][1]);
            SIMD dNZ = load1(element.dN[gp][n][2]);

            JC0 = JC0 + dNX * coordsX;
            JC1 = JC1 + dNX * coordsY;
            JC2 = JC2 + dNX * coordsZ;
            JC3 = JC3 + dNY * coordsX;
            JC4 = JC4 + dNY * coordsY;
            JC5 = JC5 + dNY * coordsZ;
            JC6 = JC6 + dNZ * coordsX;
            JC7 = JC7 + dNZ * coordsY;
            JC8 = JC8 + dNZ * coordsZ;
        }
        for (int r = 0; r < 3; ++r) {
            element.F[0 * 3 + r] = element.invJ[0 * 3 + r] * JC0 + element.invJ[1 * 3 + r] * JC3 + element.invJ[2 * 3 + r] * JC6;
            element.F[1 * 3 + r] = element.invJ[0 * 3 + r] * JC1 + element.invJ[1 * 3 + r] * JC4 + element.invJ[2 * 3 + r] * JC7;
            element.F[2 * 3 + r] = element.invJ[0 * 3 + r] * JC2 + element.invJ[1 * 3 + r] * JC5 + element.invJ[2 * 3 + r] * JC8;
        }
    }
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_INTEGRATION_DISPLACED_H_ */
