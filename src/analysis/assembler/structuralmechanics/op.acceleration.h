
#ifndef SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_ACCELERATION_H_
#define SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_ACCELERATION_H_

#include "op.elementcondition.h"

namespace espreso {

template <size_t nodes, size_t ndim> struct AccelerationKernel;

template <size_t nodes>
struct AccelerationKernel<nodes, 2>: ElementCondition {
    AccelerationKernel(const ElementCondition &base): ElementCondition(base) {}

    template <typename Element>
    void simd(Element &element, size_t gp)
    {
        switch (behaviour) {
        case StructuralMechanicsGlobalSettings::ELEMENT_BEHAVIOUR::PLANE_STRAIN:
        case StructuralMechanicsGlobalSettings::ELEMENT_BEHAVIOUR::PLANE_STRESS:
        case StructuralMechanicsGlobalSettings::ELEMENT_BEHAVIOUR::PLANE_STRESS_WITH_THICKNESS:
        {
            for (size_t n = 0; n < nodes; ++n) {
                SIMD scale = element.thickness.gp * element.det * load1(element.w[gp]) * element.ecf.density * load1(element.N[gp][n]);
                element.f[0 * nodes + n] = element.f[0 * nodes + n] + scale * element.ecf.acceleration[0];
                element.f[1 * nodes + n] = element.f[1 * nodes + n] + scale * element.ecf.acceleration[1];
            }
        } break;
        case StructuralMechanicsGlobalSettings::ELEMENT_BEHAVIOUR::AXISYMMETRIC:
        {
            for (size_t n = 0; n < nodes; ++n) {
                SIMD scale = element.det * load1(element.w[gp]) * load1(2 * M_PI) * element.coords.gp[0] * element.ecf.density * load1(element.N[gp][n]);
                element.f[1 * nodes + n] = element.f[1 * nodes + n] + scale * element.ecf.acceleration[1];
            }
        } break;
        }
    }
};

template <size_t nodes>
struct AccelerationKernel<nodes, 3>: ElementCondition {
    AccelerationKernel(const ElementCondition &base): ElementCondition(base) {}

    template <typename Element>
    void simd(Element &element, size_t gp)
    {
        for (size_t n = 0; n < nodes; ++n) {
            SIMD scale = element.det * load1(element.w[gp]) * element.ecf.density * load1(element.N[gp][n]);
            element.f[0 * nodes + n] = element.f[0 * nodes + n] + scale * element.ecf.acceleration[0];
            element.f[1 * nodes + n] = element.f[1 * nodes + n] + scale * element.ecf.acceleration[1];
            element.f[2 * nodes + n] = element.f[2 * nodes + n] + scale * element.ecf.acceleration[2];
        }
    }
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_ACCELERATION_H_ */
