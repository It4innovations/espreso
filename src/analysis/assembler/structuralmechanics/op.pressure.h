
#ifndef SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_PRESSURE_H_
#define SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_PRESSURE_H_

#include "analysis/assembler/general/boundarycondition.h"
#include "config/ecf/physics/structuralmechanics.h"

namespace espreso {

struct Pressure: public SubKernel {
    BoundaryCondition pressure;
    BoundaryCondition direction;

    Pressure()
    : behaviour(StructuralMechanicsGlobalSettings::ELEMENT_BEHAVIOUR::PLANE_STRAIN)
    {
        action = SubKernel::ASSEMBLE | SubKernel::REASSEMBLE | SubKernel::ITERATION;
    }

    void activate(ECFExpression &pressure, ECFExpressionVector &direction, StructuralMechanicsGlobalSettings::ELEMENT_BEHAVIOUR behaviour)
    {
        this->behaviour = behaviour;
        this->pressure.activate(&pressure);
        this->direction.activate(&direction);
        this->isconst = this->pressure.isconst && this->direction.isconst;
        this->isactive = 1;
        this->needCoordinates = this->pressure.needCoordinates && this->direction.needCoordinates;
    }

    StructuralMechanicsGlobalSettings::ELEMENT_BEHAVIOUR behaviour;
};

template <size_t nodes, size_t ndim> struct PressureKernel;

template <size_t nodes>
struct PressureKernel<nodes, 2>: Pressure {
    PressureKernel(const Pressure &base): Pressure(base) {}

    template <typename Element>
    void simd(Element &element, size_t gp)
    {
        switch (behaviour) {
        case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRAIN:
        case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRESS:
        case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRESS_WITH_THICKNESS:
            for (size_t n = 0; n < nodes; ++n) {
                SIMD pressure = element.det * load1(element.w[gp]) * load1(element.N[gp][n]) * element.ecf.pressure.pressure;
                element.f[0 * nodes + n] = element.f[0 * nodes + n] + pressure * element.ecf.pressure.direction[0];
                element.f[1 * nodes + n] = element.f[1 * nodes + n] + pressure * element.ecf.pressure.direction[1];
            }
            break;
        case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::AXISYMMETRIC:
            for (size_t n = 0; n < nodes; ++n) {
                SIMD pressure = element.det * load1(element.w[gp]) * load1(2 * M_PI) * element.coords.gp[0] * load1(element.N[gp][n]) * element.ecf.pressure.pressure;
                element.f[0 * nodes + n] = element.f[0 * nodes + n] + pressure * element.ecf.pressure.direction[0];
                element.f[1 * nodes + n] = element.f[1 * nodes + n] + pressure * element.ecf.pressure.direction[1];
            }
            break;
        }
    }
};

template <size_t nodes>
struct PressureKernel<nodes, 3>: Pressure {
    PressureKernel(const Pressure &base): Pressure(base) {}

    template <typename Element>
    void simd(Element &element, size_t gp)
    {
        for (size_t n = 0; n < nodes; ++n) {
            SIMD pressure = element.det * load1(element.w[gp]) * load1(element.N[gp][n]) * element.ecf.pressure.pressure;
            element.f[0 * nodes + n] = element.f[0 * nodes + n] + pressure * element.ecf.pressure.direction[0];
            element.f[1 * nodes + n] = element.f[1 * nodes + n] + pressure * element.ecf.pressure.direction[1];
            element.f[2 * nodes + n] = element.f[2 * nodes + n] + pressure * element.ecf.pressure.direction[2];
        }
    }
};

}



#endif /* SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_PRESSURE_H_ */
