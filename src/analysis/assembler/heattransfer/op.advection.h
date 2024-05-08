
#ifndef SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_HEATTRANSFER_OP_ADVECTION_H_
#define SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_HEATTRANSFER_OP_ADVECTION_H_

#include "analysis/assembler/general/boundarycondition.h"

namespace espreso {

struct Advection: BoundaryCondition {
    double sigma;

    Advection()
    : sigma(0)
    {
        isconst = false;
        action = SubKernel::ASSEMBLE | SubKernel::REASSEMBLE;
    }

    void activate(ECFExpressionVector *expression, double sigma)
    {
        BoundaryCondition::activate(expression);
        this->sigma = sigma;
    }
};

template <size_t nodes, size_t ndim> struct AdvectionKernel;

template <size_t nodes> struct AdvectionKernel<nodes, 2>: Advection {
    AdvectionKernel(const Advection &base): Advection(base) {}

    template <typename Element>
    void simd(Element &element, size_t gp)
    {
        SIMD ux = element.ecf.advection[0] * element.ecf.density * element.ecf.heatCapacity;
        SIMD uy = element.ecf.advection[1] * element.ecf.density * element.ecf.heatCapacity;

        SIMD usq = ux * ux + uy * uy;
        SIMD besq;
        for (size_t n = 0; n < nodes; ++n) {
            SIMD be = ux * element.dND[n][0] + uy * element.dND[n][1];
            besq = besq + be * be;
        }

        SIMD unorm = sqrt(usq), benorm = sqrt(besq);
        SIMD urnorm = positive_guarded_recip(unorm), bernorm = positive_guarded_recip(benorm);

        SIMD C05 = load1(.5);
        SIMD C1 = load1(1);
        SIMD C2 = load1(2);
        SIMD he = C2 * unorm * bernorm;
        SIMD rhe = benorm * C05 * urnorm;
        SIMD Pe = C2 * element.conductivity[0] * rhe * urnorm;
        SIMD tau = max(SIMD(), C1 - Pe);
        SIMD advection = C05 * he * tau * urnorm;

        SIMD stabilization = load1(sigma) * he * unorm;
        element.conductivity[0] = element.conductivity[0] + stabilization;
        element.conductivity[3] = element.conductivity[3] + stabilization;

        for (size_t n = 0; n < nodes; ++n) {
            SIMD scale = element.det * load1(element.w[gp]) * (load1(element.N[gp][n]) + advection * (ux * element.dND[n][0] + uy * element.dND[n][1]));
            for (size_t m = 0; m < nodes; ++m) {
                element.K[n * nodes + m] = element.K[n * nodes + m] + scale * (ux * element.dND[m][0] + uy * element.dND[m][1]);
            }
        }
    }
};

template <size_t nodes> struct AdvectionKernel<nodes, 3>: Advection {
    AdvectionKernel(const Advection &base): Advection(base) {}

    template <typename Element>
    void simd(Element &element, size_t gp)
    {
        SIMD ux = element.ecf.advection[0] * element.ecf.density * element.ecf.heatCapacity;
        SIMD uy = element.ecf.advection[1] * element.ecf.density * element.ecf.heatCapacity;
        SIMD uz = element.ecf.advection[2] * element.ecf.density * element.ecf.heatCapacity;

        SIMD usq = ux * ux + uy * uy + uz * uz;
        SIMD besq;
        for (size_t n = 0; n < nodes; ++n) {
            SIMD be = ux * element.dND[n][0] + uy * element.dND[n][1] + uz * element.dND[n][2];
            besq = besq + be * be;
        }

        SIMD unorm = sqrt(usq), benorm = sqrt(besq);
        SIMD urnorm = positive_guarded_recip(unorm), bernorm = positive_guarded_recip(benorm);

        SIMD C05 = load1(.5);
        SIMD C1 = load1(1);
        SIMD C2 = load1(2);
        SIMD he = C2 * unorm * bernorm;
        SIMD rhe = benorm * C05 * urnorm;
        SIMD Pe = C2 * element.conductivity[0] * rhe * urnorm;
        SIMD tau = max(SIMD(), C1 - Pe);
        SIMD advection = C05 * he * tau * urnorm;

        SIMD stabilization = load1(sigma) * he * unorm;
        element.conductivity[0] = element.conductivity[0] + stabilization;
        element.conductivity[4] = element.conductivity[4] + stabilization;
        element.conductivity[8] = element.conductivity[8] + stabilization;

        for (size_t n = 0; n < nodes; ++n) {
            SIMD scale = element.det * load1(element.w[gp]) * (load1(element.N[gp][n]) + advection * (ux * element.dND[n][0] + uy * element.dND[n][1] + uz * element.dND[n][2]));
            for (size_t m = 0; m < nodes; ++m) {
                element.K[n * nodes + m] = element.K[n * nodes + m] + scale * (ux * element.dND[m][0] + uy * element.dND[m][1] + uz * element.dND[m][2]);
            }
        }
    }
};

}



#endif /* SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_HEATTRANSFER_OP_ADVECTION_H_ */
