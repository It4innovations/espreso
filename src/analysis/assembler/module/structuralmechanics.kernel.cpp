
#include "structuralmechanics.h"
#include "analysis/assembler/module/assembler.hpp"
#include "analysis/assembler/subkernel/structuralmechanics/acceleration.h"
#include "analysis/assembler/subkernel/structuralmechanics/angularvelocity.h"
#include "analysis/assembler/subkernel/structuralmechanics/normalpressure.h"

#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.hpp"
#include "esinfo/envinfo.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/boundaryregionstore.h"

#include "analysis/scheme/steadystate.h"
#include "math/physics/matrix_distributed.h"

#include <numeric>
#include <algorithm>

namespace espreso {

template <size_t gps, size_t ndim, enum ElasticityModel model, class Physics> struct SetElasticity;

template <size_t gps, size_t ndim, class Physics> struct SetElasticity<gps, ndim, ElasticityModel::ISOTROPIC, Physics> {
	static void analyze(StructuralMechanics::SubKernels &subkernels)
	{
		subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				subkernels.elasticity.configuration->young_modulus.get(0, 0).evaluator,
				[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.youngModulus[gp][s] = value; }));
		subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				subkernels.elasticity.configuration->poisson_ratio.get(0, 0).evaluator,
				[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.poissonRatio[gp][s] = value; }));
	}
};

template <size_t gps, class Physics> struct SetElasticity<gps, 2, ElasticityModel::ORTHOTROPIC, Physics> {
	static void analyze(StructuralMechanics::SubKernels &subkernels)
	{
		subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				subkernels.elasticity.configuration->young_modulus.get(0, 0).evaluator,
				[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.youngModulus[gp][0][s] = value; }));
		subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				subkernels.elasticity.configuration->young_modulus.get(1, 1).evaluator,
				[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.youngModulus[gp][1][s] = value; }));

		subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				subkernels.elasticity.configuration->poisson_ratio.get(0, 0).evaluator,
				[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.poissonRatio[gp][0][s] = value; }));
		subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				subkernels.elasticity.configuration->poisson_ratio.get(1, 1).evaluator,
				[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.poissonRatio[gp][1][s] = value; }));

		subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				subkernels.elasticity.configuration->shear_modulus.get(0, 0).evaluator,
				[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.shearModulus[gp][0][s] = value; }));
		subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				subkernels.elasticity.configuration->shear_modulus.get(1, 1).evaluator,
				[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.shearModulus[gp][1][s] = value; }));
	}
};

template <size_t gps, class Physics> struct SetElasticity<gps, 3, ElasticityModel::ORTHOTROPIC, Physics> {
	static void analyze(StructuralMechanics::SubKernels &subkernels)
	{
		subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				subkernels.elasticity.configuration->young_modulus.get(0, 0).evaluator,
				[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.youngModulus[gp][0][s] = value; }));
		subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				subkernels.elasticity.configuration->young_modulus.get(1, 1).evaluator,
				[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.youngModulus[gp][1][s] = value; }));
		subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				subkernels.elasticity.configuration->young_modulus.get(2, 2).evaluator,
				[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.youngModulus[gp][2][s] = value; }));

		subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				subkernels.elasticity.configuration->poisson_ratio.get(0, 0).evaluator,
				[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.poissonRatio[gp][0][s] = value; }));
		subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				subkernels.elasticity.configuration->poisson_ratio.get(1, 1).evaluator,
				[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.poissonRatio[gp][1][s] = value; }));
		subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				subkernels.elasticity.configuration->poisson_ratio.get(2, 2).evaluator,
				[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.poissonRatio[gp][2][s] = value; }));

		subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				subkernels.elasticity.configuration->shear_modulus.get(0, 0).evaluator,
				[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.shearModulus[gp][0][s] = value; }));
		subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				subkernels.elasticity.configuration->shear_modulus.get(1, 1).evaluator,
				[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.shearModulus[gp][1][s] = value; }));
		subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				subkernels.elasticity.configuration->shear_modulus.get(2, 2).evaluator,
				[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.shearModulus[gp][2][s] = value; }));
	}
};

template <size_t gps, class Physics> struct SetElasticity<gps, 2, ElasticityModel::ANISOTROPIC, Physics> {
	static void analyze(StructuralMechanics::SubKernels &subkernels)
	{
		// TODO
	}
};

template <size_t gps, class Physics> struct SetElasticity<gps, 3, ElasticityModel::ANISOTROPIC, Physics> {
	static void analyze(StructuralMechanics::SubKernels &subkernels)
	{
		// TODO
	}
};

template <size_t gps, size_t ndim, class Physics> struct SetPlasticity {
	static void analyze(StructuralMechanics::SubKernels &subkernels)
	{
		// TODO
	}
};

template <size_t gps, size_t ndim, class Physics> struct SetTranslation;

template <size_t gps, class Physics> struct SetTranslation<gps, 2, Physics> {
	static void analyze(StructuralMechanics::SubKernels &subkernels)
	{
		if (subkernels.coosystem.configuration) {
			switch (subkernels.coosystem.type) {
			case CoordinateSystemConfiguration::TYPE::CARTESIAN:
				subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
						subkernels.coosystem.configuration->rotation.z.evaluator,
						[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.center[gp][0][s] = value; }));
				break;
			case CoordinateSystemConfiguration::TYPE::CYLINDRICAL:
				subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
						subkernels.coosystem.configuration->center.x.evaluator,
						[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.center[gp][0][s] = value; }));
				subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
						subkernels.coosystem.configuration->center.y.evaluator,
						[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.center[gp][1][s] = value; }));
				break;
			}
		}
	}
};

template <size_t gps, class Physics> struct SetTranslation<gps, 3, Physics> {
	static void analyze(StructuralMechanics::SubKernels &subkernels)
	{
		if (subkernels.coosystem.configuration) {
			switch (subkernels.coosystem.type) {
			case CoordinateSystemConfiguration::TYPE::CARTESIAN:
				subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
						subkernels.coosystem.configuration->rotation.x.evaluator,
						[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.center[gp][0][s] = value; }));
				subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
						subkernels.coosystem.configuration->rotation.y.evaluator,
						[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.center[gp][1][s] = value; }));
				subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
						subkernels.coosystem.configuration->rotation.z.evaluator,
						[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.center[gp][2][s] = value; }));
				break;
			case CoordinateSystemConfiguration::TYPE::CYLINDRICAL:
			case CoordinateSystemConfiguration::TYPE::SPHERICAL:
				subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
						subkernels.coosystem.configuration->center.x.evaluator,
						[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.center[gp][0][s] = value; }));
				subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
						subkernels.coosystem.configuration->center.y.evaluator,
						[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.center[gp][1][s] = value; }));
				subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
						subkernels.coosystem.configuration->center.z.evaluator,
						[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.center[gp][2][s] = value; }));
				break;
			}
		}
	}
};

template <size_t gps, size_t ndim, class Physics> struct SetThickness {
	static void analyze(StructuralMechanics::SubKernels &subkernels)
	{

	}
};

template <size_t gps, class Physics> struct SetThickness<gps, 2, Physics> {
	static void analyze(StructuralMechanics::SubKernels &subkernels)
	{
		if (subkernels.thickness.isactive) {
			subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
					subkernels.thickness.expression->evaluator,
					[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.thickness[gp][s] = value; }));
		}
	}
};

template <size_t gps, class Physics> struct SetMaterial {
	static void analyze(StructuralMechanics::SubKernels &subkernels)
	{
		subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				subkernels.material.configuration->density.evaluator,
				[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.density[gp][s] = value; }));
		subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
					subkernels.material.configuration->heat_capacity.evaluator,
					[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.heatCapacity[gp][s] = value; }));
	}
};

template <size_t gps, size_t ndim, class Physics> struct SetAcceleration;

template <size_t gps, class Physics> struct SetAcceleration<gps, 2, Physics> {
	static void analyze(StructuralMechanics::SubKernels &subkernels)
	{
		if (subkernels.acceleration.isactive) {
			subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
					subkernels.acceleration.expressionVector->x.evaluator,
					[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.acceleration[gp][0][s] = value; }));
			subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
					subkernels.acceleration.expressionVector->y.evaluator,
					[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.acceleration[gp][1][s] = value; }));
		}
	}
};

template <size_t gps, class Physics> struct SetAcceleration<gps, 3, Physics> {
	static void analyze(StructuralMechanics::SubKernels &subkernels)
	{
		if (subkernels.acceleration.isactive) {
			subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
					subkernels.acceleration.expressionVector->x.evaluator,
					[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.acceleration[gp][0][s] = value; }));
			subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
					subkernels.acceleration.expressionVector->y.evaluator,
					[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.acceleration[gp][1][s] = value; }));
			subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
					subkernels.acceleration.expressionVector->z.evaluator,
					[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.acceleration[gp][2][s] = value; }));
		}
	}
};

template <size_t gps, size_t ndim, enum Behaviour behaviour, class Physics> struct SetAngularVelocity;

template <size_t gps, class Physics> struct SetAngularVelocity<gps, 2, Behaviour::PLANE, Physics> {
	static void analyze(StructuralMechanics::SubKernels &subkernels)
	{
		if (subkernels.angularVelocity.isactive) {
			subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
					subkernels.angularVelocity.expressionVector->z.evaluator,
					[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.angularVelocity[gp][0][s] = value; }));
		}
	}
};

template <size_t gps, class Physics> struct SetAngularVelocity<gps, 2, Behaviour::AXISYMMETRIC, Physics> {
	static void analyze(StructuralMechanics::SubKernels &subkernels)
	{
		if (subkernels.angularVelocity.isactive) {
			subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
					subkernels.angularVelocity.expressionVector->y.evaluator,
					[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.angularVelocity[gp][0][s] = value; }));
		}
	}
};

template <size_t gps, class Physics> struct SetAngularVelocity<gps, 3, Behaviour::VOLUME, Physics> {
	static void analyze(StructuralMechanics::SubKernels &subkernels)
	{
		if (subkernels.angularVelocity.isactive) {
			subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
					subkernels.angularVelocity.expressionVector->x.evaluator,
					[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.angularVelocity[gp][0][s] = value; }));
			subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
					subkernels.angularVelocity.expressionVector->y.evaluator,
					[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.angularVelocity[gp][1][s] = value; }));
			subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
					subkernels.angularVelocity.expressionVector->z.evaluator,
					[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.angularVelocity[gp][2][s] = value; }));
		}
	}
};

template <Element::CODE code, size_t nodes, size_t gps, size_t ndim, size_t edim, enum Behaviour behaviour, enum ElasticityModel ecfmodel, enum ElasticityModel model>
void preprocess(StructuralMechanics::SubKernels &subkernels)
{
	typedef StructuralMechanicsElementDescriptor<nodes, gps, ndim, edim, behaviour, ecfmodel, model> Physics;
	SetThickness<gps, ndim, Physics>::analyze(subkernels);
	SetMaterial<gps, Physics>::analyze(subkernels);
	SetElasticity<gps, ndim, ecfmodel, Physics>::analyze(subkernels);
	if (subkernels.coosystem.isactive) {
		SetTranslation<gps, ndim, Physics>::analyze(subkernels);
	}
	SetPlasticity<gps, ndim, Physics>::analyze(subkernels);
	if (subkernels.plasticity.isactive) {
		subkernels.plasticity.scale.resize(gps * subkernels.chunks * SIMD::size);
		subkernels.plasticity.eps.resize(gps * subkernels.chunks * SIMD::size * 6);
		subkernels.plasticity.xi.resize(gps * subkernels.chunks * SIMD::size);
	}
	SetAcceleration<gps, ndim, Physics>::analyze(subkernels);
	SetAngularVelocity<gps, ndim, behaviour, Physics>::analyze(subkernels);

	BasisKernel<code, nodes, gps, edim> basis(subkernels.basis);
	CoordinatesKernel<nodes, gps, ndim, Physics> coordinates(subkernels.coordinates);
	IntegrationKernel<nodes, gps, ndim, edim, Physics> integration(subkernels.integration);

	typename Physics::Element element;
	basis.simd(element);
	SIMD volume;
	for (esint c = 0; c < subkernels.chunks; ++c) {
		coordinates.simd(element);
		integration.simd(element);
		for (size_t gp = 0; gp < gps; ++gp) {
			volume = volume + element.det[gp] * load1(element.w[gp]);
		}
	}

	subkernels.esize = sizeof(typename Physics::Element);
	for (size_t s = 0; s < SIMD::size; ++s) {
		subkernels.volume += volume[s];
	}
}

template <Element::CODE code, size_t nodes, size_t gps, size_t ndim, size_t edim, enum Behaviour behaviour, enum ElasticityModel ecfmodel, enum ElasticityModel model>
void compute(StructuralMechanics::SubKernels &subkernels, Assembler::Action action)
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
	PlasticityKernel<gps, ndim, ecfmodel, Physics> plasticity(subkernels.plasticity);
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

template <Element::CODE code, size_t nodes, size_t gps, size_t ndim, size_t edim, enum Behaviour behaviour, enum ElasticityModel ecfmodel, enum ElasticityModel model>
void fill(StructuralMechanics::SubKernels &subkernels)
{
	typedef StructuralMechanicsElementDescriptor<nodes, gps, ndim, edim, behaviour, ecfmodel, model> Physics;
	typename Physics::Element element;

	MatricFillerKernel<nodes, Physics> K(subkernels.Kfiller);
	VectorFillerKernel<nodes, Physics> RHS(subkernels.RHSfiller);

	for (esint c = 0; c < subkernels.chunks; ++c) {
		if (K.isactive) {
			K.simd(element);
		}
		if (RHS.isactive) {
			RHS.simd(element);
		}
	}
}

template <Element::CODE code, size_t nodes, size_t gps, size_t ndim, size_t edim>
void preprocess(StructuralMechanics::BoundarySubKernels &subkernels)
{
	typedef StructuralMechanicsBoundaryDescriptor<nodes, gps, ndim, edim> Physics;
	typename Physics::Element element;
	if (subkernels.normalPressure.expression) {
		subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				subkernels.normalPressure.expression->evaluator,
				[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.normalPressure[gp][s] = value; }));
	}

	BasisKernel<code, nodes, gps, edim> basis(subkernels.basis);
	CoordinatesKernel<nodes, gps, ndim, Physics> coordinates(subkernels.coordinates);
	IntegrationKernel<nodes, gps, ndim, edim, Physics> integration(subkernels.integration);

	basis.simd(element);
	SIMD surface;
	for (esint c = 0; c < subkernels.chunks; ++c) {
		coordinates.simd(element);
		integration.simd(element);
		for (size_t gp = 0; gp < gps; ++gp) {
			surface = surface + element.det[gp] * load1(element.w[gp]);
		}
	}

	subkernels.esize = sizeof(typename Physics::Element);
	for (size_t s = 0; s < SIMD::size; ++s) {
		subkernels.surface += surface[s];
	}
}

template <Element::CODE code, size_t nodes, size_t gps, size_t ndim, size_t edim>
void compute(StructuralMechanics::BoundarySubKernels &subkernels, Assembler::Action action)
{
	typedef StructuralMechanicsBoundaryDescriptor<nodes, gps, ndim, edim> Physics;
	typename Physics::Element element;

	BasisKernel<code, nodes, gps, edim> basis(subkernels.basis);
	CoordinatesKernel<nodes, gps, ndim, Physics> coordinates(subkernels.coordinates);
	ThicknessToNodes<nodes, ndim, Physics> thickness(subkernels.thickness);
	IntegrationKernel<nodes, gps, ndim, edim, Physics> integration(subkernels.integration);
	NormalPressureKernel<nodes, gps, ndim, edim, Physics> normalPressure(subkernels.normalPressure);

	std::vector<ExternalGPsExpression<gps, Physics>*> nonconst;
	for (size_t i = 0; i < subkernels.expressions.size(); ++i) {
		if (subkernels.expressions[i]->evaluator->isConst()) {
			dynamic_cast<ExternalGPsExpression<gps, Physics>*>(subkernels.expressions[i])->simd(element);
		} else {
			nonconst.push_back(dynamic_cast<ExternalGPsExpression<gps, Physics>*>(subkernels.expressions[i]));
		}
	}

	basis.simd(element);
	thickness.setActiveness(action);
	normalPressure.setActiveness(action);

	for (esint c = 0; c < subkernels.chunks; ++c) {
		coordinates.simd(element);
//		if (c == 0) printf("coordinates ");
		if (thickness.isactive) {
			thickness.simd(element);
//			if (c == 0) printf("thickness ");
		}
		integration.simd(element);
		if (normalPressure.isactive) {
			normalPressure.simd(element);
		}
	}
}

template <Element::CODE code, size_t nodes, size_t gps, size_t ndim, size_t edim>
void fill(const StructuralMechanics::BoundarySubKernels &subkernels)
{
	typedef StructuralMechanicsBoundaryDescriptor<nodes, gps, ndim, edim> Physics;
	typename Physics::Element element;

	VectorFillerKernel<nodes, Physics> RHS(subkernels.RHSfiller);

	for (esint c = 0; c < subkernels.chunks; ++c) {
		if (RHS.isactive) {
			RHS.simd(element);
		}
	}
}

template <size_t ndim> void initDirichlet(StructuralMechanics::BoundarySubKernels &subkernels);

template <>
void initDirichlet<2>(StructuralMechanics::BoundarySubKernels &subkernels)
{
	typedef StructuralMechanicsBoundaryDescriptor<1, 1, 2, 0> Physics;
	if (subkernels.displacement.expression) {
		subkernels.expressions.push_back(new ExternalNodeExpression<1, Physics>(
				subkernels.displacement.expression->x.evaluator,
				[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.displacement[0][0][s] = value; }));
		subkernels.expressions.push_back(new ExternalNodeExpression<1, Physics>(
				subkernels.displacement.expression->y.evaluator,
				[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.displacement[0][1][s] = value; }));
	}
}

template <>
void initDirichlet<3>(StructuralMechanics::BoundarySubKernels &subkernels)
{
	typedef StructuralMechanicsBoundaryDescriptor<1, 1, 3, 0> Physics;
	if (subkernels.displacement.expression) {
		subkernels.expressions.push_back(new ExternalNodeExpression<1, Physics>(
				subkernels.displacement.expression->x.evaluator,
				[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.displacement[0][0][s] = value; }));
		subkernels.expressions.push_back(new ExternalNodeExpression<1, Physics>(
				subkernels.displacement.expression->y.evaluator,
				[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.displacement[0][1][s] = value; }));
		subkernels.expressions.push_back(new ExternalNodeExpression<1, Physics>(
				subkernels.displacement.expression->z.evaluator,
				[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.displacement[0][2][s] = value; }));
	}
}

template <size_t ndim>
void dirichlet(const StructuralMechanics::BoundarySubKernels &subkernels)
{
	typedef StructuralMechanicsBoundaryDescriptor<1, 1, ndim, 0> Physics;
	typename Physics::Element element;

	CoordinatesKernel<1, 1, ndim, Physics> coordinates(subkernels.coordinates);
	VectorSetterKernel<1, Physics> set(subkernels.dirichlet, [] (auto &element, size_t &n, size_t &d, size_t &s) { return element.displacement[n][d][s]; });

	std::vector<ExternalNodeExpression<1, Physics>*> nonconst;
	for (size_t i = 0; i < subkernels.expressions.size(); ++i) {
		if (subkernels.expressions[i]->evaluator->isConst()) {
			dynamic_cast<ExternalNodeExpression<1, Physics>*>(subkernels.expressions[i])->simd(element);
		} else {
			nonconst.push_back(dynamic_cast<ExternalNodeExpression<1, Physics>*>(subkernels.expressions[i]));
		}
	}

	for (esint c = 0; c < subkernels.chunks; ++c) {
		coordinates.simd(element);
		for (size_t i = 0; i < nonconst.size(); ++i) {
			nonconst[i]->simd(element);
		}
		set.simd(element);
	}
}

template <Element::CODE code, size_t nodes, size_t gps, size_t ndim, size_t edim, Behaviour behaviour, enum ElasticityModel ecfmodel, enum ElasticityModel model>
void runAction(StructuralMechanics::SubKernels &subkernels, Assembler::Action action)
{
	switch (action) {
	case Assembler::Action::PREPROCESS: preprocess<code, nodes, gps, ndim, edim, behaviour, ecfmodel, model>(subkernels); break;
	case Assembler::Action::FILL: fill<code, nodes, gps, ndim, edim, behaviour, ecfmodel, model>(subkernels); break;
	default: compute<code, nodes, gps, ndim, edim, behaviour, ecfmodel, model>(subkernels, action); break;
	}
}

template <Element::CODE code, size_t nodes, size_t gps, size_t ndim, size_t edim>
void runAction(StructuralMechanics::BoundarySubKernels &subkernels, Assembler::Action action)
{
	switch (action) {
	case Assembler::Action::PREPROCESS: preprocess<code, nodes, gps, ndim, edim>(subkernels); break;
	case Assembler::Action::FILL: fill<code, nodes, gps, ndim, edim>(subkernels); break;
	default: compute<code, nodes, gps, ndim, edim>(subkernels, action); break;
	}
}

template <Element::CODE code, size_t nodes, size_t gps, size_t ndim, size_t edim, enum Behaviour behaviour>
void runElasticity(StructuralMechanics::SubKernels &subkernels, Assembler::Action action)
{
	switch (subkernels.elasticity.configuration->model) {
	case LinearElasticPropertiesConfiguration::MODEL::ISOTROPIC:
		if (subkernels.coosystem.rotated || subkernels.plasticity.isactive) {
			runAction<code, nodes, gps, ndim, edim, behaviour, ElasticityModel::ISOTROPIC, ElasticityModel::SYMMETRIC>(subkernels, action);
		} else {
			runAction<code, nodes, gps, ndim, edim, behaviour, ElasticityModel::ISOTROPIC, ElasticityModel::ISOTROPIC>(subkernels, action);
		}
		break;
	case LinearElasticPropertiesConfiguration::MODEL::ORTHOTROPIC:
		if (subkernels.coosystem.rotated) {
			runAction<code, nodes, gps, ndim, edim, behaviour, ElasticityModel::ORTHOTROPIC, ElasticityModel::SYMMETRIC>(subkernels, action);
		} else {
			runAction<code, nodes, gps, ndim, edim, behaviour, ElasticityModel::ORTHOTROPIC, ElasticityModel::ORTHOTROPIC>(subkernels, action);
		}
		break;
	case LinearElasticPropertiesConfiguration::MODEL::ANISOTROPIC:
		runAction<code, nodes, gps, ndim, edim, behaviour, ElasticityModel::ANISOTROPIC, ElasticityModel::ANISOTROPIC>(subkernels, action);
		break;
	}
}

template <enum Behaviour behaviour>
static void runElement2D(StructuralMechanics::SubKernels &subkernels, Assembler::Action action)
{
	switch (subkernels.code) {
//	case static_cast<size_t>(Element::CODE::TRIANGLE3): runElasticity<Element::CODE::TRIANGLE3, 3, StructuralMechanicsGPC::TRIANGLE3, 2, 2, behaviour>(subkernels, action); break;
//	case static_cast<size_t>(Element::CODE::TRIANGLE6): runElasticity<Element::CODE::TRIANGLE6, 6, StructuralMechanicsGPC::TRIANGLE6, 2, 2, behaviour>(subkernels, action); break;
	case static_cast<size_t>(Element::CODE::SQUARE4  ): runElasticity<Element::CODE::SQUARE4  , 4, StructuralMechanicsGPC::SQUARE4  , 2, 2, behaviour>(subkernels, action); break;
//	case static_cast<size_t>(Element::CODE::SQUARE8  ): runElasticity<Element::CODE::SQUARE8  , 8, StructuralMechanicsGPC::SQUARE8  , 2, 2, behaviour>(subkernels, action); break;
	}
}

static void runElement3D(StructuralMechanics::SubKernels &subkernels, Assembler::Action action)
{
	switch (subkernels.code) {
//	case static_cast<size_t>(Element::CODE::TETRA4   ): runElasticity<Element::CODE::TETRA4   ,  4, StructuralMechanicsGPC::TETRA4    , 3, 3, Behaviour::VOLUME>(subkernels, action); break;
//	case static_cast<size_t>(Element::CODE::TETRA10  ): runElasticity<Element::CODE::TETRA10  , 10, StructuralMechanicsGPC::TETRA10   , 3, 3, Behaviour::VOLUME>(subkernels, action); break;
//	case static_cast<size_t>(Element::CODE::PYRAMID5 ): runElasticity<Element::CODE::PYRAMID5 ,  5, StructuralMechanicsGPC::PYRAMID5  , 3, 3, Behaviour::VOLUME>(subkernels, action); break;
//	case static_cast<size_t>(Element::CODE::PYRAMID13): runElasticity<Element::CODE::PYRAMID13, 13, StructuralMechanicsGPC::PYRAMID13 , 3, 3, Behaviour::VOLUME>(subkernels, action); break;
//	case static_cast<size_t>(Element::CODE::PRISMA6  ): runElasticity<Element::CODE::PRISMA6  ,  6, StructuralMechanicsGPC::PRISMA6   , 3, 3, Behaviour::VOLUME>(subkernels, action); break;
//	case static_cast<size_t>(Element::CODE::PRISMA15 ): runElasticity<Element::CODE::PRISMA15 , 15, StructuralMechanicsGPC::PRISMA15  , 3, 3, Behaviour::VOLUME>(subkernels, action); break;
	case static_cast<size_t>(Element::CODE::HEXA8    ): runElasticity<Element::CODE::HEXA8    ,  8, StructuralMechanicsGPC::HEXA8     , 3, 3, Behaviour::VOLUME>(subkernels, action); break;
//	case static_cast<size_t>(Element::CODE::HEXA20   ): runElasticity<Element::CODE::HEXA20   , 20, StructuralMechanicsGPC::HEXA20    , 3, 3, Behaviour::VOLUME>(subkernels, action); break;
	}
}

template <size_t ndim, size_t edim> void runBoundary(StructuralMechanics::BoundarySubKernels &subkernels, Assembler::Action action);

template <> void runBoundary<2, 1>(StructuralMechanics::BoundarySubKernels &subkernels, Assembler::Action action)
{
	switch (subkernels.code) {
	case static_cast<size_t>(Element::CODE::LINE2): runAction<Element::CODE::LINE2, 2, StructuralMechanicsGPC::LINE2, 2, 1>(subkernels, action); break;
	case static_cast<size_t>(Element::CODE::LINE3): runAction<Element::CODE::LINE3, 3, StructuralMechanicsGPC::LINE3, 2, 1>(subkernels, action); break;
	}
}

template <> void runBoundary<2, 0>(StructuralMechanics::BoundarySubKernels &subkernels, Assembler::Action action)
{
	switch (action) {
	case Assembler::Action::PREPROCESS: initDirichlet<2>(subkernels); break;
	case Assembler::Action::ASSEMBLE: case Assembler::Action::REASSEMBLE: dirichlet<2>(subkernels); break;
	}
}
template <> void runBoundary<3, 2>(StructuralMechanics::BoundarySubKernels &subkernels, Assembler::Action action)
{
	switch (subkernels.code) {
	case static_cast<size_t>(Element::CODE::TRIANGLE3): runAction<Element::CODE::TRIANGLE3, 3, StructuralMechanicsGPC::TRIANGLE3, 3, 2>(subkernels, action); break;
	case static_cast<size_t>(Element::CODE::TRIANGLE6): runAction<Element::CODE::TRIANGLE6, 6, StructuralMechanicsGPC::TRIANGLE6, 3, 2>(subkernels, action); break;
	case static_cast<size_t>(Element::CODE::SQUARE4  ): runAction<Element::CODE::SQUARE4  , 4, StructuralMechanicsGPC::SQUARE4  , 3, 2>(subkernels, action); break;
	case static_cast<size_t>(Element::CODE::SQUARE8  ): runAction<Element::CODE::SQUARE8  , 8, StructuralMechanicsGPC::SQUARE8  , 3, 2>(subkernels, action); break;
	}
}

template <> void runBoundary<3, 1>(StructuralMechanics::BoundarySubKernels &subkernels, Assembler::Action action)
{
	switch (subkernels.code) {
	case static_cast<size_t>(Element::CODE::LINE2): runAction<Element::CODE::LINE2, 2, StructuralMechanicsGPC::LINE2, 3, 1>(subkernels, action); break;
	case static_cast<size_t>(Element::CODE::LINE3): runAction<Element::CODE::LINE3, 3, StructuralMechanicsGPC::LINE3, 3, 1>(subkernels, action); break;
	}
}

template <> void runBoundary<3, 0>(StructuralMechanics::BoundarySubKernels &subkernels, Assembler::Action action)
{
	switch (action) {
	case Assembler::Action::PREPROCESS: initDirichlet<3>(subkernels); break;
	case Assembler::Action::ASSEMBLE: case Assembler::Action::REASSEMBLE: dirichlet<3>(subkernels); break;
	}
}

void StructuralMechanics::run(Action action, size_t interval)
{
	switch (info::mesh->dimension) {
	case 3: runElement3D(subkernels[interval], action); break;
	case 2:
		switch (settings.element_behaviour) {
		case StructuralMechanicsGlobalSettings::ELEMENT_BEHAVIOUR::PLANE_STRAIN:
		case StructuralMechanicsGlobalSettings::ELEMENT_BEHAVIOUR::PLANE_STRESS:
		case StructuralMechanicsGlobalSettings::ELEMENT_BEHAVIOUR::PLANE_STRESS_WITH_THICKNESS:
			runElement2D<Behaviour::PLANE>(subkernels[interval], action); break;
		case StructuralMechanicsGlobalSettings::ELEMENT_BEHAVIOUR::AXISYMMETRIC:
			runElement2D<Behaviour::AXISYMMETRIC>(subkernels[interval], action); break;
		}
	}
}

void StructuralMechanics::run(Action action, size_t region, size_t interval)
{
	switch (info::mesh->dimension) {
	case 2:
		switch (info::mesh->boundaryRegions[region]->dimension) {
		case 0: runBoundary<2, 0>(boundary[region][interval], action); break;
		case 1: runBoundary<2, 1>(boundary[region][interval], action); break;
		} break;
	case 3:
		switch (info::mesh->boundaryRegions[region]->dimension) {
		case 0: runBoundary<3, 0>(boundary[region][interval], action); break;
		case 1: runBoundary<3, 1>(boundary[region][interval], action); break;
		case 2: runBoundary<3, 2>(boundary[region][interval], action); break;
		} break;
	}
}

}
