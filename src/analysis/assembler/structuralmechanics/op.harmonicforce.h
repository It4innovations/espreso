
#ifndef SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_HARMONICFORCE_H_
#define SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_HARMONICFORCE_H_

#include "analysis/assembler/general/boundarycondition.h"
#include "config/ecf/physics/structuralmechanics.h"

namespace espreso {

struct HarmonicForce: HarmonicBoundaryCondition {

    HarmonicForce()
    : behaviour(StructuralMechanicsGlobalSettings::ELEMENT_BEHAVIOUR::PLANE_STRAIN)
    {

    }

    void activate(ECFHarmonicExpressionVector *expression, StructuralMechanicsGlobalSettings::ELEMENT_BEHAVIOUR behaviour)
    {
        this->behaviour = behaviour;
        HarmonicBoundaryCondition::activate(expression);
    }

    StructuralMechanicsGlobalSettings::ELEMENT_BEHAVIOUR behaviour;
};

template <size_t ndim> struct HarmonicForceKernel;

template <>
struct HarmonicForceKernel<2>: HarmonicForce {
    HarmonicForceKernel(const HarmonicForce &base): HarmonicForce(base) {}

    template <typename Element>
    void simd(Element &element)
    {
        switch (behaviour) {
        case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRAIN:
        case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRESS:
    case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRESS_WITH_THICKNESS:
            element.f  [0] = element.f  [0] + element.ecf.harmonicForceMag[0] * element.ecf.harmonicForceCos[0];
            element.f  [1] = element.f  [1] + element.ecf.harmonicForceMag[1] * element.ecf.harmonicForceCos[1];
            element.imf[0] = element.imf[0] + element.ecf.harmonicForceMag[0] * element.ecf.harmonicForceSin[0];
            element.imf[1] = element.imf[1] + element.ecf.harmonicForceMag[1] * element.ecf.harmonicForceSin[1];
            break;
        case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::AXISYMMETRIC:
            element.f  [1] = element.f  [1] + element.ecf.harmonicForceMag[1] * element.ecf.harmonicForceCos[1];
            element.imf[1] = element.imf[1] + element.ecf.harmonicForceMag[1] * element.ecf.harmonicForceSin[1];
            break;
        }
    }
};

template <>
struct HarmonicForceKernel<3>: HarmonicForce {
    HarmonicForceKernel(const HarmonicForce &base): HarmonicForce(base) {}

    template <typename Element>
    void simd(Element &element)
    {
        element.f  [0] = element.f  [0] + element.ecf.harmonicForceMag[0] * element.ecf.harmonicForceCos[0];
        element.f  [1] = element.f  [1] + element.ecf.harmonicForceMag[1] * element.ecf.harmonicForceCos[1];
        element.f  [2] = element.f  [2] + element.ecf.harmonicForceMag[2] * element.ecf.harmonicForceCos[2];
        element.imf[0] = element.imf[0] + element.ecf.harmonicForceMag[0] * element.ecf.harmonicForceSin[0];
        element.imf[1] = element.imf[1] + element.ecf.harmonicForceMag[1] * element.ecf.harmonicForceSin[1];
        element.imf[2] = element.imf[2] + element.ecf.harmonicForceMag[2] * element.ecf.harmonicForceSin[2];
    }
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_HARMONICFORCE_H_ */
