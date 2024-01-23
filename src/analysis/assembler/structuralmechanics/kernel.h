
#ifndef SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_KERNEL_H_
#define SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_KERNEL_H_

#include "element.h"
#include "operators.h"
#include "mesh/element.h"

#include <iostream>

namespace espreso {

template <Element::CODE code, size_t nodes, size_t gps, size_t ndim, size_t edim>
void setElementKernel(StructuralMechanicsOperators &subkernels, SubKernel::Action action)
{
	typedef StructuralMechanicsElement<nodes, gps, ndim, edim> Element; Element element;

	if constexpr(ndim == 2) {
		if (subkernels.thickness.expression) {
			subkernels.expressions.node.push_back(new ExternalNodeExpression<ndim, Element>(
					subkernels.thickness.expression->evaluator,
					[] (Element &element, size_t &n, size_t &s, double value) { element.thickness.node[n][s] = value; }));
		}

		switch (subkernels.elasticity.coordinateSystem->type) {
		case CoordinateSystemConfiguration::TYPE::CARTESIAN:
			if (subkernels.elasticity.coordinateSystem->rotation.z.isset) {
				subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
					subkernels.elasticity.coordinateSystem->rotation.z.evaluator,
					[] (Element &element, size_t &gp, size_t &s, double value) { element.rotation.center[0][s] = value; }));
			}
			break;
		case CoordinateSystemConfiguration::TYPE::CYLINDRICAL:
			if (subkernels.elasticity.coordinateSystem->center.x.isset) {
				subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
					subkernels.elasticity.coordinateSystem->center.x.evaluator,
					[] (Element &element, size_t &gp, size_t &s, double value) { element.rotation.center[0][s] = value; }));
			}
			if (subkernels.elasticity.coordinateSystem->center.y.isset) {
				subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
					subkernels.elasticity.coordinateSystem->center.y.evaluator,
					[] (Element &element, size_t &gp, size_t &s, double value) { element.rotation.center[1][s] = value; }));
			}
			break;
		case CoordinateSystemConfiguration::TYPE::SPHERICAL:
			break;
		}

		subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
				subkernels.elasticity.configuration->young_modulus.get(0, 0).evaluator,
				[] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.youngModulus[0][s] = value; }));
		subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
				subkernels.elasticity.configuration->young_modulus.get(1, 1).evaluator,
				[] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.youngModulus[1][s] = value; }));

		subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
				subkernels.elasticity.configuration->poisson_ratio.get(0, 0).evaluator,
				[] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.poissonRatio[0][s] = value; }));
		subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
				subkernels.elasticity.configuration->poisson_ratio.get(1, 1).evaluator,
				[] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.poissonRatio[1][s] = value; }));

		subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
				subkernels.elasticity.configuration->shear_modulus.get(0, 0).evaluator,
				[] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.shearModulus[0][s] = value; }));
		subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
				subkernels.elasticity.configuration->shear_modulus.get(1, 1).evaluator,
				[] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.shearModulus[1][s] = value; }));

		if (subkernels.acceleration.expressionVector) {
			subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
					subkernels.acceleration.expressionVector->x.evaluator,
					[] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.acceleration[0][s] = value; }));
		}
		if (subkernels.acceleration.expressionVector) {
			subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
					subkernels.acceleration.expressionVector->y.evaluator,
					[] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.acceleration[1][s] = value; }));
		}

		if (subkernels.angularVelocity.expressionVector) {
			subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
					subkernels.angularVelocity.expressionVector->y.evaluator,
					[] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.angularVelocity[1][s] = value; }));
		}
		if (subkernels.angularVelocity.expressionVector) {
			subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
					subkernels.angularVelocity.expressionVector->z.evaluator,
					[] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.angularVelocity[2][s] = value; }));
		}
	}

	if constexpr(ndim == 3) {
		switch (subkernels.elasticity.coordinateSystem->type) {
		case CoordinateSystemConfiguration::TYPE::CARTESIAN:
			if (subkernels.elasticity.coordinateSystem->rotation.x.isset) {
				subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
					subkernels.elasticity.coordinateSystem->rotation.x.evaluator,
					[] (Element &element, size_t &gp, size_t &s, double value) { element.rotation.center[0][s] = value; }));
			}
			if (subkernels.elasticity.coordinateSystem->rotation.y.isset) {
				subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
					subkernels.elasticity.coordinateSystem->rotation.y.evaluator,
					[] (Element &element, size_t &gp, size_t &s, double value) { element.rotation.center[1][s] = value; }));
			}
			if (subkernels.elasticity.coordinateSystem->rotation.z.isset) {
				subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
					subkernels.elasticity.coordinateSystem->rotation.z.evaluator,
					[] (Element &element, size_t &gp, size_t &s, double value) { element.rotation.center[2][s] = value; }));
			}
			break;
		case CoordinateSystemConfiguration::TYPE::CYLINDRICAL:
			if (subkernels.elasticity.coordinateSystem->center.x.isset) {
				subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
					subkernels.elasticity.coordinateSystem->center.x.evaluator,
					[] (Element &element, size_t &gp, size_t &s, double value) { element.rotation.center[0][s] = value; }));
			}
			if (subkernels.elasticity.coordinateSystem->center.y.isset) {
				subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
					subkernels.elasticity.coordinateSystem->center.y.evaluator,
					[] (Element &element, size_t &gp, size_t &s, double value) { element.rotation.center[1][s] = value; }));
			}
			break;
		case CoordinateSystemConfiguration::TYPE::SPHERICAL:
			if (subkernels.elasticity.coordinateSystem->center.x.isset) {
				subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
					subkernels.elasticity.coordinateSystem->center.x.evaluator,
					[] (Element &element, size_t &gp, size_t &s, double value) { element.rotation.center[0][s] = value; }));
			}
			if (subkernels.elasticity.coordinateSystem->center.y.isset) {
				subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
					subkernels.elasticity.coordinateSystem->center.y.evaluator,
					[] (Element &element, size_t &gp, size_t &s, double value) { element.rotation.center[1][s] = value; }));
			}
			if (subkernels.elasticity.coordinateSystem->center.z.isset) {
				subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
					subkernels.elasticity.coordinateSystem->center.z.evaluator,
					[] (Element &element, size_t &gp, size_t &s, double value) { element.rotation.center[2][s] = value; }));
			}
			break;
		}

		subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
				subkernels.elasticity.configuration->young_modulus.get(0, 0).evaluator,
				[] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.youngModulus[0][s] = value; }));
		subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
				subkernels.elasticity.configuration->young_modulus.get(1, 1).evaluator,
				[] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.youngModulus[1][s] = value; }));
		subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
				subkernels.elasticity.configuration->young_modulus.get(2, 2).evaluator,
				[] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.youngModulus[2][s] = value; }));

		subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
				subkernels.elasticity.configuration->poisson_ratio.get(0, 0).evaluator,
				[] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.poissonRatio[0][s] = value; }));
		subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
				subkernels.elasticity.configuration->poisson_ratio.get(1, 1).evaluator,
				[] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.poissonRatio[1][s] = value; }));
		subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
				subkernels.elasticity.configuration->poisson_ratio.get(2, 2).evaluator,
				[] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.poissonRatio[2][s] = value; }));

		subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
				subkernels.elasticity.configuration->shear_modulus.get(0, 0).evaluator,
				[] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.shearModulus[0][s] = value; }));
		subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
				subkernels.elasticity.configuration->shear_modulus.get(1, 1).evaluator,
				[] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.shearModulus[1][s] = value; }));
		subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
				subkernels.elasticity.configuration->shear_modulus.get(2, 2).evaluator,
				[] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.shearModulus[2][s] = value; }));

		if (subkernels.plasticity.isactive) {
			subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
					subkernels.plasticity.configuration->initial_yield_stress.evaluator,
					[] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.initialYieldStress[s] = value; }));
			subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
					subkernels.plasticity.configuration->isotropic_hardening.evaluator,
					[] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.isotropicHardening[s] = value; }));
			subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
					subkernels.plasticity.configuration->kinematic_hardening.evaluator,
					[] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.kinematicHardening[s] = value; }));
			subkernels.plasticity.smallStrainTensorPlastic.resize(gps * (subkernels.chunks + 1) * SIMD::size * 6);
			subkernels.plasticity.xi.resize(gps * (subkernels.chunks + 1) * SIMD::size * 6 + SIMD::size);
		}

		if (subkernels.acceleration.expressionVector) {
			subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
					subkernels.acceleration.expressionVector->x.evaluator,
					[] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.acceleration[0][s] = value; }));
		}
		if (subkernels.acceleration.expressionVector) {
			subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
					subkernels.acceleration.expressionVector->y.evaluator,
					[] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.acceleration[1][s] = value; }));
		}
		if (subkernels.acceleration.expressionVector) {
			subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
					subkernels.acceleration.expressionVector->z.evaluator,
					[] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.acceleration[2][s] = value; }));
		}
		if (subkernels.angularVelocity.expressionVector) {
			subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
					subkernels.angularVelocity.expressionVector->x.evaluator,
					[] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.angularVelocity[0][s] = value; }));
		}
		if (subkernels.angularVelocity.expressionVector) {
			subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
					subkernels.angularVelocity.expressionVector->y.evaluator,
					[] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.angularVelocity[1][s] = value; }));
		}
		if (subkernels.angularVelocity.expressionVector) {
			subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
					subkernels.angularVelocity.expressionVector->z.evaluator,
					[] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.angularVelocity[2][s] = value; }));
		}
	}

	if (subkernels.material.configuration->density.isset) {
		subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
				subkernels.material.configuration->density.evaluator,
				[] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.density[s] = value; }));
	}
	if (subkernels.material.configuration->heat_capacity.isset) {
		subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
				subkernels.material.configuration->heat_capacity.evaluator,
				[] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.heatCapacity[s] = value; }));
	}

	BasisKernel<code, nodes, gps, edim> basis(subkernels.basis);
	CoordinatesKernel<nodes, ndim> coordinates(subkernels.coordinates);
	IntegrationKernel<nodes, ndim, edim> integration(subkernels.integration);
	ThicknessToNodes<nodes, ndim> thickness(subkernels.thickness);

	SIMD volume;
	basis.simd(element);
	for (size_t c = 0; c < subkernels.chunks; ++c) {
		coordinates.simd(element);
		thickness.simd(element);

		for (size_t gp = 0; gp < gps; ++gp) {
			integration.simd(element, gp);
			volume = volume + element.det * load1(element.w[gp]);
		}
	}

	subkernels.esize = sizeof(Element);
	for (size_t s = 0; s < SIMD::size; ++s) {
		subkernels.volume += volume[s];
	}
}

template <Element::CODE code, size_t nodes, size_t gps, size_t ndim, size_t edim>
void runElementKernel(StructuralMechanicsOperators &subkernels, SubKernel::Action action)
{
	typedef StructuralMechanicsElement<nodes, gps, ndim, edim> Element; Element element;

	BasisKernel<code, nodes, gps, edim> basis(subkernels.basis);
	CoordinatesKernel<nodes, ndim> coordinates(subkernels.coordinates);
	CoordinatesToGPsKernel<nodes, ndim> coordinatesToGPs(subkernels.coordinates);
	ThicknessToGp<nodes, ndim> thickness(subkernels.thickness);
	TemperatureKernel<nodes> temperature(subkernels.temperature);
	TemperatureToGPsKernel<nodes> temperatureToGPs(subkernels.temperature);
	IntegrationKernel<nodes, ndim, edim> integration(subkernels.integration);
	DisplacementKernel<nodes, ndim> displacement(subkernels.displacement);
	SmallStrainTensorKernel<nodes, ndim> smallStrainTensor(subkernels.smallStrainTensor);
	ElasticityKernel<ndim> elasticity(subkernels.elasticity);
	PlasticityKernel<nodes, ndim> plasticity(subkernels.plasticity, action);
	MatrixElasticityKernel<nodes, ndim> K(subkernels.K);
	AccelerationKernel<nodes, ndim> acceleration(subkernels.acceleration);
	AngularVelocityKernel<nodes, ndim> angularVelocity(subkernels.angularVelocity);
	SigmaKernel<nodes, ndim> sigma(subkernels.sigma);
	StressKernel<nodes, gps, ndim> stress(subkernels.stress);
	MatricFillerKernel<nodes> outK(subkernels.Kfiller);
	RHSFillerKernel<nodes> outRHS(subkernels.RHSfiller);
	RHSFillerKernel<nodes> outNRHS(subkernels.RHSfiller);

	struct {
		std::vector<ExternalNodeExpression<ndim, Element>*> node;
		std::vector<ExternalGPsExpression<ndim, Element>*> gp;
	} nonconst;

	for (size_t i = 0; i < subkernels.expressions.node.size(); ++i) {
		ExternalNodeExpression<ndim, Element>* exp = dynamic_cast<ExternalNodeExpression<ndim, Element>*>(subkernels.expressions.node[i]);
		if (subkernels.expressions.node[i]->evaluator->isConst()) {
			for (size_t n = 0; n < nodes; ++n) {
				exp->simd(element, n);
			}
		} else {
			nonconst.node.push_back(exp);
		}
	}

	for (size_t i = 0; i < subkernels.expressions.gp.size(); ++i) {
		ExternalGPsExpression<ndim, Element>* exp = dynamic_cast<ExternalGPsExpression<ndim, Element>*>(subkernels.expressions.gp[i]);
		if (subkernels.expressions.gp[i]->evaluator->isConst()) {
			for (size_t gp = 0; gp < gps; ++gp) {
				exp->simd(element, gp);
			}
		} else {
			nonconst.gp.push_back(exp);
		}
	}

	// pre-processing of possible constant parameters from ecf
	basis.simd(element);
	elasticity.simd(element);
	thickness.simd(element, 0);

	coordinatesToGPs.setActiveness(action);
	thickness.setActiveness(action);
	temperature.setActiveness(action);
	elasticity.setActiveness(action);
	plasticity.setActiveness(action);
	displacement.setActiveness(action);
	smallStrainTensor.setActiveness(action);
	K.setActiveness(action);
	acceleration.setActiveness(action);
	angularVelocity.setActiveness(action);
	sigma.setActiveness(action);
	stress.setActiveness(action);

	outK.setActiveness(action);
	outRHS.setActiveness(action);
	outNRHS.setActiveness(action);

	for (size_t c = 0; c < subkernels.chunks; ++c) {
		if (sigma.isactive) {
			sigma.reset(element);
		}
		coordinates.simd(element);
		if (temperature.isactive) {
			temperature.simd(element);
		}
		if (displacement.isactive) {
			displacement.simd(element);
		}
		for (size_t i = 0; i < nonconst.node.size(); ++i) {
			for (size_t n = 0; n < nodes; ++n) {
				nonconst.node[i]->simd(element, n);
			}
		}

		for (size_t gp = 0; gp < gps; ++gp) {
			integration.simd(element, gp);

			if (coordinatesToGPs.isactive) {
				coordinatesToGPs.simd(element, gp);
			}
			if (temperatureToGPs.isactive) {
				temperatureToGPs.simd(element, gp);
			}
			if (thickness.isactive) {
				thickness.simd(element, gp);
			}
			for (size_t i = 0; i < nonconst.gp.size(); ++i) {
				for (size_t n = 0; n < nodes; ++n) {
					nonconst.gp[i]->simd(element, n);
				}
			}
			if (elasticity.isactive) {
				elasticity.simd(element, gp);
			}
			if (smallStrainTensor.isactive) {
				smallStrainTensor.simd(element, gp);
			}
			if (plasticity.isactive) {
				plasticity.simd(element, gp);
			}
			if (K.isactive) {
				K.simd(element, gp);
			}

			if (acceleration.isactive) {
				acceleration.simd(element, gp);
			}
			if (angularVelocity.isactive) {
				angularVelocity.simd(element, gp);
			}
			if (sigma.isactive) {
				sigma.simd(element, gp);
			}
		}

		if (outK.isactive) {
			outK.simd(element.K);
		}
		if (outRHS.isactive) {
			outRHS.simd(element.f);
		}
		if (outNRHS.isactive) {
			outNRHS.simd(element.nf);
		}
		if (stress.isactive) {
			stress.simd(element);
		}
	}
}

template <Element::CODE code, size_t nodes, size_t gps, size_t ndim, size_t edim>
void setBoundaryKernel(StructuralMechanicsBoundaryOperators &subkernels, SubKernel::Action action)
{
	typedef StructuralMechanicsBoundary<nodes, gps, ndim, edim> Element; Element element;

	if (subkernels.normalPressure.expression) {
		subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
				subkernels.normalPressure.expression->evaluator,
				[] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.normalPressure[s] = value; }));
	}

	BasisKernel<code, nodes, gps, edim> basis(subkernels.basis);
	CoordinatesKernel<nodes, ndim> coordinates(subkernels.coordinates);
	IntegrationKernel<nodes, ndim, edim> integration(subkernels.integration);

	SIMD surface;
	basis.simd(element);
	for (size_t c = 0; c < subkernels.chunks; ++c) {
		coordinates.simd(element);
		for (size_t gp = 0; gp < gps; ++gp) {
			integration.simd(element, gp);
			surface = surface + element.det * load1(element.w[gp]);
		}
	}

	subkernels.esize = sizeof(Element);
	for (size_t s = 0; s < SIMD::size; ++s) {
		subkernels.surface += surface[s];
	}
}

template <Element::CODE code, size_t nodes, size_t gps, size_t ndim, size_t edim>
void runBoundaryKernel(const StructuralMechanicsBoundaryOperators &subkernels, SubKernel::Action action)
{
	typedef StructuralMechanicsBoundary<nodes, gps, ndim, edim> Element; Element element;

	BasisKernel<code, nodes, gps, edim> basis(subkernels.basis);
	CoordinatesKernel<nodes, ndim> coordinates(subkernels.coordinates);
	CoordinatesToGPsKernel<nodes, ndim> coordinatesToGPs(subkernels.coordinates);
	ThicknessFromNodes<nodes, ndim> thickness(subkernels.thickness);
	ThicknessToGp<nodes, ndim> thicknessToGPs(subkernels.thickness);
	IntegrationKernel<nodes, ndim, edim> integration(subkernels.integration);
	NormalPressureKernel<nodes, ndim> normalPressure(subkernels.normalPressure);
	RHSFillerKernel<nodes> outRHS(subkernels.RHSfiller);

	std::vector<ExternalGPsExpression<ndim, Element>*> nonconst;
	for (size_t i = 0; i < subkernels.expressions.gp.size(); ++i) {
		ExternalGPsExpression<ndim, Element>* exp = dynamic_cast<ExternalGPsExpression<ndim, Element>*>(subkernels.expressions.gp[i]);
		if (subkernels.expressions.gp[i]->evaluator->isConst()) {
			exp->simd(element, 0);
		} else {
			nonconst.push_back(exp);
		}
	}

	basis.simd(element);
	thickness.setActiveness(action);

	outRHS.setActiveness(action);

	for (size_t c = 0; c < subkernels.chunks; ++c) {
		coordinates.simd(element);

		if (thickness.isactive) {
			thickness.simd(element);
		}

		for (size_t gp = 0; gp < gps; ++gp) {
			if (coordinatesToGPs.isactive) {
				coordinatesToGPs.simd(element, gp);
			}

			for (size_t i = 0; i < nonconst.size(); ++i) {
				nonconst[i]->simd(element, gp);
			}

			if (thicknessToGPs.isactive) {
				thicknessToGPs.simd(element, gp);
			}

			integration.simd(element, gp);

			if (normalPressure.isactive) {
				normalPressure.simd(element, gp);
			}
		}

		if (outRHS.isactive) {
			outRHS.simd(element.f);
		}
	}
}

template <size_t ndim>
void setDirichletKernel(StructuralMechanicsBoundaryOperators &subkernels, SubKernel::Action action)
{
	typedef StructuralMechanicsDirichlet<ndim> Element; Element element;
	if (subkernels.displacement.expression) {
		if (subkernels.displacement.expression->x.isset) {
			auto setter = [] (Element &element, size_t &n, size_t &s, double value) { element.displacement.node[0][s] = value; };
			switch (info::mesh->dimension) {
			case 2: subkernels.expressions.node.push_back(new ExternalNodeExpression<2, Element>(subkernels.displacement.expression->x.evaluator, setter)); break;
			case 3: subkernels.expressions.node.push_back(new ExternalNodeExpression<3, Element>(subkernels.displacement.expression->x.evaluator, setter)); break;
			}
		}
		if (subkernels.displacement.expression->y.isset) {
			auto setter = [] (Element &element, size_t &n, size_t &s, double value) { element.displacement.node[1][s] = value; };
			switch (info::mesh->dimension) {
			case 2: subkernels.expressions.node.push_back(new ExternalNodeExpression<2, Element>(subkernels.displacement.expression->y.evaluator, setter)); break;
			case 3: subkernels.expressions.node.push_back(new ExternalNodeExpression<3, Element>(subkernels.displacement.expression->y.evaluator, setter)); break;
			}
		}
		if (subkernels.displacement.expression->z.isset) {
			auto setter = [] (Element &element, size_t &n, size_t &s, double value) { element.displacement.node[2][s] = value; };
			subkernels.expressions.node.push_back(new ExternalNodeExpression<3, Element>(subkernels.displacement.expression->z.evaluator, setter));
		}
	}
}

template <size_t ndim>
void runDirichletKernel(const StructuralMechanicsBoundaryOperators &subkernels, SubKernel::Action action)
{
	typedef StructuralMechanicsDirichlet<ndim> Element; Element element;

	CoordinatesKernel<1, ndim> coordinates(subkernels.coordinates);
	VectorSetterKernel<1, Element> set(subkernels.dirichlet, [] (auto &element, size_t &n, size_t &d, size_t &s) { return element.displacement.node[d][s]; });

	std::vector<ExternalNodeExpression<ndim, Element>*> nonconst;
	for (size_t i = 0; i < subkernels.expressions.node.size(); ++i) {
		ExternalNodeExpression<ndim, Element>* exp = dynamic_cast<ExternalNodeExpression<ndim, Element>*>(subkernels.expressions.node[i]);
		if (subkernels.expressions.node[i]->evaluator->isConst()) {
			exp->simd(element, 0);
		} else {
			nonconst.push_back(exp);
		}
	}

	for (size_t c = 0; c < subkernels.chunks; ++c) {
		coordinates.simd(element);
		for (size_t i = 0; i < nonconst.size(); ++i) {
			nonconst[i]->simd(element, 0);
		}
		set.simd(element);
	}
}

}



#endif /* SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_KERNEL_H_ */
