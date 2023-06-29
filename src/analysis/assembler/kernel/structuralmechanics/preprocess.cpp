
#include "analysis/assembler/module/structuralmechanics.h"
#include "analysis/assembler/module/assembler.hpp"
#include "analysis/assembler/subkernel/structuralmechanics/acceleration.h"
#include "analysis/assembler/subkernel/structuralmechanics/angularvelocity.h"
#include "analysis/assembler/subkernel/structuralmechanics/normalpressure.h"

#include <numeric>
#include <algorithm>

namespace espreso {

template <size_t gps, size_t ndim, enum ElasticityModel model, class Physics> struct SetElasticity;

template <size_t gps, size_t ndim, class Physics> struct SetElasticity<gps, ndim, ElasticityModel::ISOTROPIC, Physics> {
	static void analyze(StructuralMechanicsSubKernelsList &subkernels)
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
	static void analyze(StructuralMechanicsSubKernelsList &subkernels)
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
	static void analyze(StructuralMechanicsSubKernelsList &subkernels)
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
	static void analyze(StructuralMechanicsSubKernelsList &subkernels)
	{
		// TODO
	}
};

template <size_t gps, class Physics> struct SetElasticity<gps, 3, ElasticityModel::ANISOTROPIC, Physics> {
	static void analyze(StructuralMechanicsSubKernelsList &subkernels)
	{
		// TODO
	}
};

template <size_t gps, size_t ndim, class Physics> struct SetPlasticity {
	static void analyze(StructuralMechanicsSubKernelsList &subkernels)
	{
		// TODO
	}
};

template <size_t gps, size_t ndim, class Physics> struct SetTranslation;

template <size_t gps, class Physics> struct SetTranslation<gps, 2, Physics> {
	static void analyze(StructuralMechanicsSubKernelsList &subkernels)
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
	static void analyze(StructuralMechanicsSubKernelsList &subkernels)
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
	static void analyze(StructuralMechanicsSubKernelsList &subkernels)
	{

	}
};

template <size_t gps, class Physics> struct SetThickness<gps, 2, Physics> {
	static void analyze(StructuralMechanicsSubKernelsList &subkernels)
	{
		if (subkernels.thickness.isactive) {
			subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
					subkernels.thickness.expression->evaluator,
					[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.thickness[gp][s] = value; }));
		}
	}
};

template <size_t gps, class Physics> struct SetMaterial {
	static void analyze(StructuralMechanicsSubKernelsList &subkernels)
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
	static void analyze(StructuralMechanicsSubKernelsList &subkernels)
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
	static void analyze(StructuralMechanicsSubKernelsList &subkernels)
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
	static void analyze(StructuralMechanicsSubKernelsList &subkernels)
	{
		if (subkernels.angularVelocity.isactive) {
			subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
					subkernels.angularVelocity.expressionVector->z.evaluator,
					[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.angularVelocity[gp][0][s] = value; }));
		}
	}
};

template <size_t gps, class Physics> struct SetAngularVelocity<gps, 2, Behaviour::AXISYMMETRIC, Physics> {
	static void analyze(StructuralMechanicsSubKernelsList &subkernels)
	{
		if (subkernels.angularVelocity.isactive) {
			subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
					subkernels.angularVelocity.expressionVector->y.evaluator,
					[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.angularVelocity[gp][0][s] = value; }));
		}
	}
};

template <size_t gps, class Physics> struct SetAngularVelocity<gps, 3, Behaviour::VOLUME, Physics> {
	static void analyze(StructuralMechanicsSubKernelsList &subkernels)
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
void preprocess(StructuralMechanicsSubKernelsList &subkernels)
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
		subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				subkernels.plasticity.configuration->initial_yield_stress.evaluator,
				[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.initialYieldStress[gp][s] = value; }));
		subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				subkernels.plasticity.configuration->isotropic_hardening.evaluator,
				[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.isotropicHardening[gp][s] = value; }));
		subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				subkernels.plasticity.configuration->kinematic_hardening.evaluator,
				[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.kinematicHardening[gp][s] = value; }));
		subkernels.plasticity.smallStrainTensorPlastic.resize(gps * (subkernels.chunks + 1) * SIMD::size * 6);
		subkernels.plasticity.xi.resize(gps * (subkernels.chunks + 1) * SIMD::size * 6 + SIMD::size);
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
void fill(StructuralMechanicsSubKernelsList &subkernels)
{
	typedef StructuralMechanicsElementDescriptor<nodes, gps, ndim, edim, behaviour, ecfmodel, model> Physics;
	typename Physics::Element element;

	MatricFillerKernel<nodes, Physics> K(subkernels.Kfiller);
	VectorFillerKernel<nodes, Physics> RHS(subkernels.RHSfiller), nRHS(subkernels.nRHSfiller);

	for (esint c = 0; c < subkernels.chunks; ++c) {
		if (K.isactive) {
			K.simd(element);
		}
		if (RHS.isactive) {
			RHS.simd(element);
		}
		if (nRHS.isactive) {
			nRHS.simd(element);
		}
	}
}

template <Element::CODE code, size_t nodes, size_t gps, size_t ndim, size_t edim, Behaviour behaviour, enum ElasticityModel ecfmodel, enum ElasticityModel model>
void runAction(StructuralMechanicsSubKernelsList &subkernels, Assembler::Action action)
{
	switch (action) {
	case Assembler::Action::PREPROCESS: preprocess<code, nodes, gps, ndim, edim, behaviour, ecfmodel, model>(subkernels); break;
	case Assembler::Action::FILL: fill<code, nodes, gps, ndim, edim, behaviour, ecfmodel, model>(subkernels); break;
	}
}

template <Element::CODE code, size_t nodes, size_t gps, size_t ndim, size_t edim>
void runAction(StructuralMechanicsBoundarySubKernelsList &subkernels, Assembler::Action action)
{
	switch (action) {
	case Assembler::Action::PREPROCESS: preprocess<code, nodes, gps, ndim, edim>(subkernels); break;
	case Assembler::Action::FILL: fill<code, nodes, gps, ndim, edim>(subkernels); break;
	}
}

template <Element::CODE code, size_t nodes, size_t gps, size_t ndim, size_t edim, enum Behaviour behaviour>
void runElasticity(StructuralMechanicsSubKernelsList &subkernels, Assembler::Action action)
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
static void runElement2D(StructuralMechanicsSubKernelsList &subkernels, Assembler::Action action)
{
	switch (subkernels.code) {
	case static_cast<size_t>(Element::CODE::TRIANGLE3): runElasticity<Element::CODE::TRIANGLE3, 3, StructuralMechanicsGPC::TRIANGLE3, 2, 2, behaviour>(subkernels, action); break;
	case static_cast<size_t>(Element::CODE::TRIANGLE6): runElasticity<Element::CODE::TRIANGLE6, 6, StructuralMechanicsGPC::TRIANGLE6, 2, 2, behaviour>(subkernels, action); break;
	case static_cast<size_t>(Element::CODE::SQUARE4  ): runElasticity<Element::CODE::SQUARE4  , 4, StructuralMechanicsGPC::SQUARE4  , 2, 2, behaviour>(subkernels, action); break;
	case static_cast<size_t>(Element::CODE::SQUARE8  ): runElasticity<Element::CODE::SQUARE8  , 8, StructuralMechanicsGPC::SQUARE8  , 2, 2, behaviour>(subkernels, action); break;
	}
}

static void runElement3D(StructuralMechanicsSubKernelsList &subkernels, Assembler::Action action)
{
	switch (subkernels.code) {
	case static_cast<size_t>(Element::CODE::TETRA4   ): runElasticity<Element::CODE::TETRA4   ,  4, StructuralMechanicsGPC::TETRA4    , 3, 3, Behaviour::VOLUME>(subkernels, action); break;
	case static_cast<size_t>(Element::CODE::TETRA10  ): runElasticity<Element::CODE::TETRA10  , 10, StructuralMechanicsGPC::TETRA10   , 3, 3, Behaviour::VOLUME>(subkernels, action); break;
	case static_cast<size_t>(Element::CODE::PYRAMID5 ): runElasticity<Element::CODE::PYRAMID5 ,  5, StructuralMechanicsGPC::PYRAMID5  , 3, 3, Behaviour::VOLUME>(subkernels, action); break;
	case static_cast<size_t>(Element::CODE::PYRAMID13): runElasticity<Element::CODE::PYRAMID13, 13, StructuralMechanicsGPC::PYRAMID13 , 3, 3, Behaviour::VOLUME>(subkernels, action); break;
	case static_cast<size_t>(Element::CODE::PRISMA6  ): runElasticity<Element::CODE::PRISMA6  ,  6, StructuralMechanicsGPC::PRISMA6   , 3, 3, Behaviour::VOLUME>(subkernels, action); break;
	case static_cast<size_t>(Element::CODE::PRISMA15 ): runElasticity<Element::CODE::PRISMA15 , 15, StructuralMechanicsGPC::PRISMA15  , 3, 3, Behaviour::VOLUME>(subkernels, action); break;
	case static_cast<size_t>(Element::CODE::HEXA8    ): runElasticity<Element::CODE::HEXA8    ,  8, StructuralMechanicsGPC::HEXA8     , 3, 3, Behaviour::VOLUME>(subkernels, action); break;
	case static_cast<size_t>(Element::CODE::HEXA20   ): runElasticity<Element::CODE::HEXA20   , 20, StructuralMechanicsGPC::HEXA20    , 3, 3, Behaviour::VOLUME>(subkernels, action); break;
	}
}

void StructuralMechanics::runPreprocess(Action action, size_t interval)
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

}
