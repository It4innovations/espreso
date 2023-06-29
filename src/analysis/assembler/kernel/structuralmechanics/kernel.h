
#ifndef SRC_ANALYSIS_ASSEMBLER_KERNEL_STRUCTURALMECHANICS_KERNEL_H_
#define SRC_ANALYSIS_ASSEMBLER_KERNEL_STRUCTURALMECHANICS_KERNEL_H_

#include "analysis/assembler/module/assembler.h"
#include "analysis/assembler/subkernel/thickness.h"
#include "analysis/assembler/subkernel/structuralmechanics/acceleration.h"
#include "analysis/assembler/subkernel/structuralmechanics/angularvelocity.h"
#include "analysis/assembler/subkernel/structuralmechanics/normalpressure.h"
#include "analysis/assembler/kernel/structuralmechanics/subkernellist.h"

namespace espreso {

template <Element::CODE code, size_t nodes, size_t gps, size_t ndim, size_t edim, enum Behaviour behaviour, enum ElasticityModel ecfmodel, enum ElasticityModel model>
void compute(StructuralMechanicsSubKernelsList &subkernels, Assembler::Action action)
{
	typedef StructuralMechanicsElementDescriptor<nodes, gps, ndim, edim, behaviour, ecfmodel, model> Physics;
	typename Physics::Element element;

	BasisKernel<code, nodes, gps, edim> basis(subkernels.basis);
	CoordinatesKernel<nodes, gps, ndim, Physics> coordinates(subkernels.coordinates);
	ThicknessToNodes<nodes, ndim, Physics> thickness(subkernels.thickness);
	TemperatureKernel<nodes, gps, Physics> temperature(subkernels.temperature);
	IntegrationKernel<nodes, gps, ndim, edim, Physics> integration(subkernels.integration);
	DisplacementKernel<nodes, gps, ndim, Physics> displacement(subkernels.displacement);
	ElasticityKernel<gps, ndim, ecfmodel, Physics> elasticity(subkernels.elasticity);
	PlasticityKernel<nodes, gps, ndim, ecfmodel, Physics> plasticity(subkernels.plasticity, action);
	StructuralMechanicsCoordinateSystemKernel<gps, ndim, ecfmodel, model, Physics> coosystem(subkernels.coosystem);
	StructuralMechanicsStiffness<nodes, gps, ndim, behaviour, model, Physics> K(subkernels.K);
	AccelerationKernel<nodes, gps, behaviour, Physics> acceleration(subkernels.acceleration);
	AngularVelocityKernel<nodes, gps, behaviour, Physics> angularVelocity(subkernels.angularVelocity);
	StressKernel<nodes, gps, ndim, model, Physics> stress(subkernels.stress);

	std::vector<ExternalGPsExpression<gps, Physics>*> nonconst;
	for (size_t i = 0; i < subkernels.expressions.size(); ++i) {
		if (subkernels.expressions[i]->evaluator->isConst()) {
			dynamic_cast<ExternalGPsExpression<gps, Physics>*>(subkernels.expressions[i])->simd(element);
		} else {
			nonconst.push_back(dynamic_cast<ExternalGPsExpression<gps, Physics>*>(subkernels.expressions[i]));
		}
	}

	basis.simd(element);
	if (coosystem.isactive) {
		coosystem.simd(element);
	}

	thickness.setActiveness(action);
	temperature.setActiveness(action);
	elasticity.setActiveness(action);
	plasticity.setActiveness(action);
	coosystem.setActiveness(action);
	displacement.setActiveness(action);
	K.setActiveness(action);
	acceleration.setActiveness(action);
	angularVelocity.setActiveness(action);
	stress.setActiveness(action);

//	printf("sub-kernels: ");
	for (esint c = 0; c < subkernels.chunks; ++c) {
		coordinates.simd(element);
//		if (c == 0) printf("coordinates ");
		if (temperature.isactive) {
			temperature.simd(element);
//			if (c == 0) printf("temp ");
		}
		if (thickness.isactive) {
			thickness.simd(element);
//			if (c == 0) printf("thickness ");
		}
		integration.simd(element);
		if (displacement.isactive) {
			displacement.simd(element);
//			if (c == 0) printf("displacement ");
		}
//		if (c == 0) printf("integrate ");
		if (elasticity.isactive) {
			elasticity.simd(element);
//			if (c == 0) printf("elasticity ");
		}
		if (plasticity.isactive) {
			plasticity.simd(element);
//			if (c == 0) printf("plasticity ");
		}

		if (coosystem.isactive) {
			coosystem.simd(element);
//			if (c == 0) printf("coosystem ");
		}
		if (K.isactive) {
			K.simd(element);
//			if (c == 0) printf("K ");
		}
		if (acceleration.isactive) {
			acceleration.simd(element);
		}
		if (angularVelocity.isactive) {
			angularVelocity.simd(element);
		}
		if (stress.isactive) {
			stress.simd(element);
		}
	}
//	printf("\n");
}

}

#endif /* SRC_ANALYSIS_ASSEMBLER_KERNEL_STRUCTURALMECHANICS_KERNEL_H_ */
