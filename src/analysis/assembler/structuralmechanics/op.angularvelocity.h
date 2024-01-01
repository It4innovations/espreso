
#ifndef SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_ANGULARVELOCITY_H_
#define SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_ANGULARVELOCITY_H_

#include "element.h"
#include "op.elementcondition.h"

namespace espreso {

template <size_t nodes, size_t ndim> struct AngularVelocityKernel;

template <size_t nodes>
struct AngularVelocityKernel<nodes, 2>: ElementCondition {
	AngularVelocityKernel(const ElementCondition &base): ElementCondition(base) {}

	template <typename Element>
	void simd(Element &element, size_t gp)
	{
		switch (behaviour) {
		case StructuralMechanicsGlobalSettings::ELEMENT_BEHAVIOUR::PLANE_STRAIN:
		case StructuralMechanicsGlobalSettings::ELEMENT_BEHAVIOUR::PLANE_STRESS:
		case StructuralMechanicsGlobalSettings::ELEMENT_BEHAVIOUR::PLANE_STRESS_WITH_THICKNESS:
		{
			for (size_t n = 0; n < nodes; ++n) {
				SIMD scale = element.thickness.gp * element.det * load1(element.w[gp]) * element.ecf.density * load1(element.N[gp][n]) * element.ecf.angularVelocity[2] * element.ecf.angularVelocity[2];
				element.f[0 * nodes + n] = element.f[0 * nodes + n] + scale * element.coords.gp[0];
				element.f[1 * nodes + n] = element.f[1 * nodes + n] + scale * element.coords.gp[1];
			}
		} break;
		case StructuralMechanicsGlobalSettings::ELEMENT_BEHAVIOUR::AXISYMMETRIC:
		{
			for (size_t n = 0; n < nodes; ++n) {
				SIMD scale = element.det * load1(element.w[gp]) * load1(2 * M_PI) * element.coords.gp[0] * element.ecf.density * load1(element.N[gp][n]) * element.ecf.angularVelocity[1] * element.ecf.angularVelocity[1];
				element.f[0 * nodes + n] = element.f[0 * nodes + n] + scale * element.coords.gp[0];
			}
		} break;
		}
	}
};

template <size_t nodes>
struct AngularVelocityKernel<nodes, 3>: ElementCondition {
	AngularVelocityKernel(const ElementCondition &base): ElementCondition(base) {}

	template <typename Element>
	void simd(Element &element, size_t gp)
	{
		for (size_t n = 0; n < nodes; ++n) {
			SIMD scale = element.det * load1(element.w[gp]) * element.ecf.density * load1(element.N[gp][n]);
			SIMD vx = element.ecf.angularVelocity[0] * element.ecf.angularVelocity[0];
			SIMD vy = element.ecf.angularVelocity[1] * element.ecf.angularVelocity[1];
			SIMD vz = element.ecf.angularVelocity[2] * element.ecf.angularVelocity[2];
			element.f[0 * nodes + n] = element.f[0 * nodes + n] + scale * element.coords.gp[0] * (vy + vz);
			element.f[1 * nodes + n] = element.f[1 * nodes + n] + scale * element.coords.gp[1] * (vx + vz);
			element.f[2 * nodes + n] = element.f[2 * nodes + n] + scale * element.coords.gp[2] * (vx + vy);
		}
	}
};

}



#endif /* SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_ANGULARVELOCITY_H_ */
