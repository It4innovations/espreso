
#ifndef SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_NORMALPRESSURE_H_
#define SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_NORMALPRESSURE_H_

#include "analysis/assembler/general/boundarycondition.h"
#include "config/ecf/physics/structuralmechanics.h"

namespace espreso {

struct NormalPressupre: BoundaryCondition {

	NormalPressupre()
	: behaviour(StructuralMechanicsGlobalSettings::ELEMENT_BEHAVIOUR::PLANE_STRAIN)
	{

	}

	void activate(ECFExpression *expression, StructuralMechanicsGlobalSettings::ELEMENT_BEHAVIOUR behaviour)
	{
		this->behaviour = behaviour;
		BoundaryCondition::activate(expression);
	}

	StructuralMechanicsGlobalSettings::ELEMENT_BEHAVIOUR behaviour;
};

template <size_t nodes, size_t ndim> struct NormalPressureKernel;

template <size_t nodes>
struct NormalPressureKernel<nodes, 2>: NormalPressupre {
	NormalPressureKernel(const NormalPressupre &base): NormalPressupre(base) {}

	template <typename Element>
	void simd(Element &element, size_t gp)
	{
		switch (behaviour) {
		case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRAIN:
		case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRESS:
		case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRESS_WITH_THICKNESS:
			for (size_t n = 0; n < nodes; ++n) {
				SIMD pressure = element.det * load1(element.w[gp]) * load1(element.N[gp][n]) * element.ecf.normalPressure;
				element.f[0 * nodes + n] = element.f[0 * nodes + n] + pressure * element.normal[0];
				element.f[1 * nodes + n] = element.f[1 * nodes + n] + pressure * element.normal[1];
			}
			break;
		case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::AXISYMMETRIC:
			for (size_t n = 0; n < nodes; ++n) {
				SIMD pressure = element.det * load1(element.w[gp]) * load1(2 * M_PI) * element.coords.gp[0] * load1(element.N[gp][n]) * element.ecf.normalPressure;
				element.f[0 * nodes + n] = element.f[0 * nodes + n] + pressure * element.normal[0];
				element.f[1 * nodes + n] = element.f[1 * nodes + n] + pressure * element.normal[1];
			}
			break;
		}
	}
};

template <size_t nodes>
struct NormalPressureKernel<nodes, 3>: NormalPressupre {
	NormalPressureKernel(const NormalPressupre &base): NormalPressupre(base) {}

	template <typename Element>
	void simd(Element &element, size_t gp)
	{
		for (size_t n = 0; n < nodes; ++n) {
			SIMD pressure = element.det * load1(element.w[gp]) * load1(element.N[gp][n]) * element.ecf.normalPressure;
			element.f[0 * nodes + n] = element.f[0 * nodes + n] + pressure * element.normal[0];
			element.f[1 * nodes + n] = element.f[1 * nodes + n] + pressure * element.normal[1];
			element.f[2 * nodes + n] = element.f[2 * nodes + n] + pressure * element.normal[2];
		}
	}
};

}



#endif /* SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_NORMALPRESSURE_H_ */
