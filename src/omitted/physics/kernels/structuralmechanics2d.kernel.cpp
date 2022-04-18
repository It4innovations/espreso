
#include "structuralmechanics2d.kernel.h"
#include "heattransfer2d.kernel.h"
#include "solverdataprovider/structuralmechanics2d.provider.h"
#include "esinfo/stepinfo.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.hpp"
#include "esinfo/meshinfo.h"

#include "physics/system/builder/builder.h"
#include "basis/containers/point.h"
#include "basis/evaluator/evaluator.h"
#include "math/matrix.dense.h"
#include "math/vector.dense.h"

#include <cmath>

using namespace espreso;

StructuralMechanics2DKernel::StructuralMechanics2DKernel(StructuralMechanics2DKernel *previous, PhysicsConfiguration &physics, StructuralMechanicsGlobalSettings &gsetting, StructuralMechanicsLoadStepConfiguration &configuration)
: KernelExecutor(new StructuralMechanics2DSolverDataProvider(configuration)),
  iterator(previous ? &previous->iterator : NULL, physics, gsetting, configuration, 2)
{
	boundaries.reserve(info::mesh->boundaryRegions.size());
	for (size_t i = 0; i < info::mesh->boundaryRegions.size(); ++i) {
		boundaries.emplace_back(info::mesh->boundaryRegions[i], iterator, configuration, 2);
	}
}

StructuralMechanics2DKernel::StructuralMechanics2DKernel(HeatTransfer2DKernel *previous, PhysicsConfiguration &physics, StructuralMechanicsGlobalSettings &gsetting, StructuralMechanicsLoadStepConfiguration &configuration)
: KernelExecutor(new StructuralMechanics2DSolverDataProvider(configuration)),
  iterator(previous ? &previous->iterator : NULL, physics, gsetting, configuration, 2)
{
	boundaries.reserve(info::mesh->boundaryRegions.size());
	for (size_t i = 0; i < info::mesh->boundaryRegions.size(); ++i) {
		boundaries.emplace_back(info::mesh->boundaryRegions[i], iterator, configuration, 2);
	}
}

StructuralMechanics2DKernel::~StructuralMechanics2DKernel()
{
	iterator.clear();
	for (size_t i = 0; i < info::mesh->boundaryRegions.size(); ++i) {
		boundaries[i].clear();
	}
	delete solverDataProvider;
}

void StructuralMechanics2DKernel::assembleMaterialMatrix(esint node, double *coordinates, const MaterialBaseConfiguration *mat, double time, double temp, MatrixDense &K) const
{
	double Ex = 0, Ey = 0, mi = 0;

	switch (mat->linear_elastic_properties.model) {

	case LinearElasticPropertiesConfiguration::MODEL::ISOTROPIC:
		Ex = Ey = mat->linear_elastic_properties.young_modulus.get(0, 0).evaluator->eval(Evaluator::Params().coords(2, coordinates).temp(&temp));
		mi = mat->linear_elastic_properties.poisson_ratio.get(0, 0).evaluator->eval(Evaluator::Params().coords(2, coordinates).temp(&temp));
		break;

	case LinearElasticPropertiesConfiguration::MODEL::ANISOTROPIC:
		eslog::error("Implement ANISOTROPIC MATERIAL.\n");
		break;

	case LinearElasticPropertiesConfiguration::MODEL::ORTHOTROPIC:
		Ex = mat->linear_elastic_properties.young_modulus.get(0, 0).evaluator->eval(Evaluator::Params().coords(2, coordinates).temp(&temp));
		Ey = mat->linear_elastic_properties.young_modulus.get(1, 1).evaluator->eval(Evaluator::Params().coords(2, coordinates).temp(&temp));
		mi = mat->linear_elastic_properties.poisson_ratio.get(0, 0).evaluator->eval(Evaluator::Params().coords(2, coordinates).temp(&temp));
		break;

	default:
		eslog::error("Linear elasticity 2D not supports set material model.\n");
	}

	switch (mat->linear_elastic_properties.model) {

	case LinearElasticPropertiesConfiguration::MODEL::ISOTROPIC:
	{

		switch (info::ecf->structural_mechanics_2d.element_behaviour) {

		case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRAIN:
		{
			double k = Ex * (1 - mi) / ((1 + mi) * (1 - 2 * mi));
			K(node, 0) = k * 1;
			K(node, 1) = k * 1;
			K(node, 2) = k * ((1 - 2 * mi) / (2 * (1 - mi)));
			K(node, 3) = k * (mi / (1 - mi));
			K(node, 4) = 0;
			K(node, 5) = 0;
			K(node, 6) = k * (mi / (1 - mi));
			K(node, 7) = 0;
			K(node, 8) = 0;
			return;
		}

		case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRESS:
		case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRESS_WITH_THICKNESS:
		{
			double k = Ex / (1 - mi * mi);
			K(node, 0) = k * 1;
			K(node, 1) = k * 1;
			K(node, 2) = k * ((1 -  mi) / 2);
			K(node, 3) = k * mi;
			K(node, 4) = 0;
			K(node, 5) = 0;
			K(node, 6) = k * mi;
			K(node, 7) = 0;
			K(node, 8) = 0;
			return;
		}

		case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::AXISYMMETRIC:
		{
			K.resize(K.nrows, 16);
			double k = Ex * (1 - mi) / ((1 + mi) * (1 - 2 * mi));
			K(node,  0) = k * 1;
			K(node,  1) = k * 1;
			K(node,  2) = k * ((1 - 2 * mi) / (2 * (1 - mi)));
			K(node,  3) = k * 1;

			K(node,  4) = k * (mi / (1 - mi));
			K(node,  5) = 0;
			K(node,  6) = k * (mi / (1 - mi));
			K(node,  7) = 0;
			K(node,  8) = k * (mi / (1 - mi));
			K(node,  9) = 0;

			K(node, 10) = k * (mi / (1 - mi));
			K(node, 11) = 0;
			K(node, 12) = 0;
			K(node, 13) = k * (mi / (1 - mi));
			K(node, 14) = k * (mi / (1 - mi));
			K(node, 15) = 0;
			return;
		}
		}
		break;
	}

	case LinearElasticPropertiesConfiguration::MODEL::ORTHOTROPIC:
	{
		eslog::error("IMPLEMENT: MATERIAL_MODEL::LINEAR_ELASTIC_ORTHOTROPIC.\n");
		return;
	}

	case LinearElasticPropertiesConfiguration::MODEL::ANISOTROPIC:
	{
		eslog::error("IMPLEMENT: MATERIAL_MODEL::LINEAR_ELASTIC_ANISOTROPIC.\n");
		return;
	}

	default:
		eslog::error("Structural mechanics 2D not supports set material model.\n");
	}
}

// source dX, dY

// target::
// dX   0
//  0  dY
// dY  dX
static void distribute3x2(double *target, double *source, size_t rows, size_t columns)
{
	memcpy(target                               , source          , sizeof(double) * columns);
	memcpy(target + 2 * rows * columns + columns, source          , sizeof(double) * columns);

	memcpy(target + 1 * rows * columns + columns, source + columns, sizeof(double) * columns);
	memcpy(target + 2 * rows * columns          , source + columns, sizeof(double) * columns);
}

// source dX, dY

// target::
// dX   0
//  0  dY
//  0  0
// dY  dX
static void distribute4x2(double *target, double *source, size_t rows, size_t columns)
{
	memcpy(target                               , source          , sizeof(double) * columns);
	memcpy(target + 3 * rows * columns + columns, source          , sizeof(double) * columns);

	memcpy(target + 1 * rows * columns + columns, source + columns, sizeof(double) * columns);
	memcpy(target + 3 * rows * columns          , source + columns, sizeof(double) * columns);
}

void StructuralMechanics2DKernel::processElement(const Builder &builder, const ElasticityElementIterator &iterator, InstanceFiller &filler) const
{
	esint size = iterator.element->nodes;
	filler.DOFs = 2 * size;

	const std::vector<MatrixDense> &N = *(iterator.element->N);
	const std::vector<MatrixDense> &NN = *(iterator.element->NN);
	const std::vector<MatrixDense> &dN = *(iterator.element->dN);
	const std::vector<double> &weighFactor = *(iterator.element->weighFactor);

	MatrixDense Ce(4, 4), XY(1, 2), coordinates(size, 2), J, invJ(2, 2), CeB, dND, B, precision, CEp, rhsT;
	MatrixDense K(size, 9), TE(size, 2), thickness(size, 1), inertia(size, 2), dens(size, 1);
	MatrixDense gpK(size, 9), gpTE(1, 2), gpThickness(1, 1), gpInertia(1, 2), gpDens(1, 1);
	double detJ, CP = 1, te;
	Point center;

	for (esint n = 0; n < size; n++) {
		Evaluator::Params params;
		params.coords(2, iterator.coordinates.data + 2 * n);
		params.temp(iterator.temperature.data + n);
		inertia(n, 0) = iterator.acceleration.data[2 * n + 0];
		inertia(n, 1) = iterator.acceleration.data[2 * n + 1];
		coordinates(n, 0) = iterator.coordinates.data[2 * n + 0];
		coordinates(n, 1) = iterator.coordinates.data[2 * n + 1];
		iterator.material->density.evaluator->evalVector(1, params, &dens(n, 0));

		thickness(n, 0) = iterator.thickness.data[n];
		switch (iterator.material->thermal_expansion.model) {
		case ThermalExpansionConfiguration::MODEL::ISOTROPIC:
			te = iterator.material->thermal_expansion.thermal_expansion.get(0, 0).evaluator->eval(params);
			TE(n, 0) = TE(n, 1) = (iterator.temperature.data[n] - iterator.initialTemperature.data[n]) * te;
			break;
		case ThermalExpansionConfiguration::MODEL::ORTHOTROPIC:
			te = iterator.material->thermal_expansion.thermal_expansion.get(0, 0).evaluator->eval(params);
			TE(n, 0) = (iterator.temperature.data[n] - iterator.initialTemperature.data[n]) * te;
			te = iterator.material->thermal_expansion.thermal_expansion.get(1, 1).evaluator->eval(params);
			TE(n, 1) = (iterator.temperature.data[n] - iterator.initialTemperature.data[n]) * te;
			break;
		}
		assembleMaterialMatrix(n, iterator.coordinates.data, iterator.material, step::time.current, iterator.temperature.data[n], K);
	}

	if (builder.matrices & (Builder::Request::K | Builder::Request::R)) {
		filler.Ke.resize(2 * size, 2 * size);
		filler.Ke.fill(0);
	}
	if (builder.matrices & Builder::Request::M) {
		filler.Me.resize(2 * size, 2 * size);
		filler.Me.fill(0);
	}
	if (builder.matrices & Builder::Request::R) {
		filler.Re.resize(2 * size);
		filler.Re.fill(0);
	}
	if (builder.matrices & Builder::Request::f) {
		filler.Fe.resize(2 * size);
		filler.Fe.fill(0);
	}

	for (size_t gp = 0; gp < N.size(); gp++) {
		J.multiply(dN[gp], coordinates);
		detJ = MATH::determinant2x2(J.vals);
		if (detJ < 0) { ++filler.invalid; detJ = -detJ; }
		MATH::Dense2x2inverse(J.vals, invJ.vals, detJ);

		gpThickness.multiply(N[gp], thickness);
		gpK.multiply(N[gp], K);
		dND.multiply(invJ, dN[gp]);
		gpDens.multiply(N[gp], dens);

		if (builder.matrices & Builder::Request::f) {
			gpTE.multiply(N[gp], TE);
			gpInertia.multiply(N[gp], inertia);
			XY.multiply(N[gp], coordinates);
		}

		if (builder.matrices & Builder::Request::M) {
			filler.Me.multiply(NN[gp], NN[gp], gpDens(0, 0) * detJ * weighFactor[gp] * CP, 1, true);
		}

		switch (info::ecf->structural_mechanics_2d.element_behaviour) {

		case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRESS:
		case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRAIN:
			gpThickness(0, 0) = 1;
		case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRESS_WITH_THICKNESS:

			Ce.resize(3, 3);
			Ce(0, 0) = gpK(0, 0);
			Ce(1, 1) = gpK(0, 1);
			Ce(2, 2) = gpK(0, 2);
			Ce(0, 1) = gpK(0, 3);
			Ce(0, 2) = gpK(0, 4);
			Ce(1, 2) = gpK(0, 5);
			Ce(1, 0) = gpK(0, 6);
			Ce(2, 0) = gpK(0, 7);
			Ce(2, 1) = gpK(0, 8);

			B.resize(Ce.nrows, 2 * size);
			B.fill(0);
			distribute3x2(B.vals, dND.vals, dND.nrows, dND.ncols);

			if (builder.matrices & (Builder::Request::K | Builder::Request::R)) {
				CeB.multiply(Ce, B);
				filler.Ke.multiply(B, CeB, detJ * weighFactor[gp] * gpThickness(0, 0), 1, true);
			}

			if (builder.matrices & Builder::Request::f) {
				precision.resize(Ce.nrows, 1);
				precision(0, 0) = gpTE(0, 1);
				precision(1, 0) = gpTE(0, 1);
				precision(2, 0) = 0;
				CEp.multiply(Ce, precision);
				rhsT.multiply(B, CEp, detJ * weighFactor[gp] * gpThickness(0, 0), 0, true, false);
				for (esint i = 0; i < 2 * size; i++) {
					filler.Fe[0][i] += gpDens(0, 0) * detJ * weighFactor[gp] * gpThickness(0, 0) * N[gp](0, i % size) * gpInertia(0, i / size);
					filler.Fe[0][i] += gpDens(0, 0) * detJ * weighFactor[gp] * gpThickness(0, 0) * N[gp](0, i % size) * XY(0, i / size) * pow(iterator.angularVelocity.data[2], 2);
					filler.Fe[0][i] += rhsT(i, 0);
				}
			}
			break;

		case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::AXISYMMETRIC:

			Ce.resize(4, 4);
			Ce(0, 0) = gpK(0,  0);
			Ce(1, 1) = gpK(0,  1);
			Ce(2, 2) = gpK(0,  2);
			Ce(3, 3) = gpK(0,  3);
			Ce(0, 1) = gpK(0,  4);
			Ce(0, 2) = gpK(0,  5);
			Ce(0, 3) = gpK(0,  6);
			Ce(1, 2) = gpK(0,  7);
			Ce(1, 3) = gpK(0,  8);
			Ce(2, 3) = gpK(0,  9);
			Ce(1, 0) = gpK(0, 10);
			Ce(2, 0) = gpK(0, 11);
			Ce(2, 1) = gpK(0, 12);
			Ce(3, 0) = gpK(0, 13);
			Ce(3, 1) = gpK(0, 14);
			Ce(3, 2) = gpK(0, 15);

			B.resize(Ce.nrows, 2 * size);
			B.fill(0);
			distribute4x2(B.vals, dND.vals, dND.nrows, dND.ncols);
			for(esint i = 0; i < N[gp].ncols; i++) {
				B(2, i) = N[gp](0, i) / XY(0, 0);
			}

			if (builder.matrices & (Builder::Request::K | Builder::Request::R)) {
				CeB.multiply(Ce, B);
				filler.Ke.multiply(B, CeB, detJ * weighFactor[gp] * 2 * M_PI * XY(0, 0), 1, true);
			}

			if (builder.matrices & Builder::Request::f) {
				precision.resize(Ce.nrows, 1);
				precision(0, 0) = gpTE(0, 0);
				precision(1, 0) = gpTE(0, 1);
				precision(2, 0) = precision(3, 0) = 0;
				CEp.multiply(Ce, precision);
				rhsT.multiply(B, CEp, detJ * weighFactor[gp] * 2 * M_PI * XY(0, 0), 0, true, false);
				for (esint i = 0; i < 2 * size; i++) {
					filler.Fe[0][i] += gpDens(0, 0) * detJ * weighFactor[gp] * 2 * M_PI * XY(0, 0) * N[gp](0, i % size) * gpInertia(0, i / size);
					filler.Fe[0][i] += rhsT(i, 0);
				}
				for (esint i = 0; i < size; i++) {
					filler.Fe[0][i] += gpDens(0, 0) * detJ * weighFactor[gp] * 2 * M_PI * XY(0, 0) * N[gp](0, i % size) * XY(0, 0) * pow(iterator.angularVelocity.data[1], 2);
					filler.Fe[0][size + i] += gpDens(0, 0) * detJ * weighFactor[gp] * 2 * M_PI * XY(0, 0) * N[gp](0, i % size) * XY(0, 1) * pow(iterator.angularVelocity.data[1], 2);
				}
			}
			break;
		}
	}

	if (iterator.designVariable.data) {
		double scale = iterator.minDesignVariable + pow(iterator.designVariable.data[0], iterator.penaltyFactor) * (1 - iterator.minDesignVariable);
		filler.Ke.scale(scale);
	}

	filler.insertK = builder.matrices & Builder::Request::K;
	filler.insertM = builder.matrices & Builder::Request::M;
	filler.insertC = false;
	filler.insertR = builder.matrices & Builder::Request::R;
	filler.insertF = builder.matrices & Builder::Request::f;
}

void StructuralMechanics2DKernel::processEdge(const Builder &builder, const ElasticityBoundaryIterator &iterator, InstanceFiller &filler) const
{
	esint size = iterator.element->nodes;
	filler.DOFs = 2 * iterator.element->nodes;
	if ((filler.insertF = (builder.matrices & Builder::Request::f))) {
		filler.Fe.resize(2 * size);
		filler.Fe.fill(0);
	}

	if (iterator.normalPressure.data == NULL) {
		return;
	}
	if (!(builder.matrices & (Builder::Request::K | Builder::Request::f))) {
		return;
	}

	const std::vector<MatrixDense> &N = *(iterator.element->N);
	const std::vector<MatrixDense> &dN = *(iterator.element->dN);
	const std::vector<double> &weighFactor = *(iterator.element->weighFactor);

	MatrixDense coordinates(size, 2), dND(1, 2), P(size, 1), normal(2, 1), matThickness(size, 1), XY(1, 2);
	MatrixDense gpP(1, 1), gpQ(1, 2), gpThickness(1, 1);

	for (esint n = 0; n < size; n++) {
		coordinates(n, 0) = iterator.coordinates.data[2 * n + 0];
		coordinates(n, 1) = iterator.coordinates.data[2 * n + 1];
		P(n, 0) = iterator.normalPressure.data[n];
		matThickness(n, 0) = 1; // iterator.thickness.data[n]; FIXME
	}

	for (size_t gp = 0; gp < N.size(); gp++) {
		dND.multiply(dN[gp], coordinates);
		double J = dND.norm();
		Point n(-dND(0, 1), dND(0, 0), 0);
		// e->rotateOutside(e->parentElements()[0], _mesh->coordinates(), n);
		normal(0, 0) = n.x / J;
		normal(1, 0) = n.y / J;
		gpP.multiply(N[gp], P);
		gpQ.multiply(normal, gpP);
		gpThickness.multiply(N[gp], matThickness);

		switch (info::ecf->structural_mechanics_2d.element_behaviour) {

		case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRESS:
		case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRAIN:
			gpThickness(0, 0) = 1;
		case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRESS_WITH_THICKNESS:
			for (esint i = 0; i < 2 * size; i++) {
				filler.Fe[0][i] += gpThickness(0, 0) * J * weighFactor[gp] * N[gp](0, i % size) * gpQ(0, i / size);
			}
			break;

		case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::AXISYMMETRIC:
			XY.multiply(N[gp], coordinates);
			for (esint i = 0; i < 2 * size; i++) {
				filler.Fe[0][i] += gpThickness(0, 0) * J * weighFactor[gp] * 2 * M_PI * XY(0, 0) * N[gp](0, i % size) * gpQ(0, i / size);
			}
			break;
		}
	}
}

void StructuralMechanics2DKernel::processNode(const Builder &builder, const ElasticityBoundaryIterator &iterator, InstanceFiller &filler) const
{
	filler.insertK = filler.insertF = false;
	filler.DOFs = 2;
	if (iterator.force.data == NULL) {
		return;
	}
	if (!(builder.matrices & Builder::Request::f)) {
		return;
	}

	if ((filler.insertF = (builder.matrices & Builder::Request::f))) {
		filler.Fe.resize(2);
		filler.Fe.fill(0);
	}

	if (iterator.force.data) {
		filler.Fe[0][0] = iterator.force.data[0];
		filler.Fe[0][1] = iterator.force.data[1];
	}
}

void StructuralMechanics2DKernel::elementSolution(ElasticityElementIterator &iterator)
{
	if (iterator.designVariable.data == NULL) {
		return;
	}

	esint size = iterator.element->nodes;

	const std::vector<MatrixDense> &N = *(iterator.element->N);
	const std::vector<MatrixDense> &dN = *(iterator.element->dN);
	const std::vector<double> &weighFactor = *(iterator.element->weighFactor);

	MatrixDense Ke(2 * size, 2 * size), Ce(4, 4), XY(1, 2), coordinates(size, 2), J, invJ(2, 2), CeB, dND, B;
	VectorDense disp(2 * size), Ku(2 * size);
	MatrixDense K(size, 9), thickness(size, 1);
	MatrixDense gpK(size, 9), gpThickness(1, 1);
	double detJ;

	for (esint n = 0; n < size; n++) {
		coordinates(n, 0) = iterator.coordinates.data[2 * n + 0];
		coordinates(n, 1) = iterator.coordinates.data[2 * n + 1];
		thickness(n, 0) = iterator.thickness.data[n];
		disp[n + 0 * size] = iterator.displacement.data[2 * n + 0];
		disp[n + 1 * size] = iterator.displacement.data[2 * n + 1];
		assembleMaterialMatrix(n, iterator.coordinates.data, iterator.material, step::time.current, iterator.temperature.data[n], K);
	}

	for (size_t gp = 0; gp < N.size(); gp++) {
		J.multiply(dN[gp], coordinates);
		detJ = MATH::determinant2x2(J.vals);
		if (detJ < 0) { detJ = -detJ; }
		MATH::Dense2x2inverse(J.vals, invJ.vals, detJ);

		gpThickness.multiply(N[gp], thickness);
		gpK.multiply(N[gp], K);
		dND.multiply(invJ, dN[gp]);

		switch (info::ecf->structural_mechanics_2d.element_behaviour) {

		case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRESS:
		case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRAIN:
			gpThickness(0, 0) = 1;
		case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRESS_WITH_THICKNESS:

			Ce.resize(3, 3);
			Ce(0, 0) = gpK(0, 0);
			Ce(1, 1) = gpK(0, 1);
			Ce(2, 2) = gpK(0, 2);
			Ce(0, 1) = gpK(0, 3);
			Ce(0, 2) = gpK(0, 4);
			Ce(1, 2) = gpK(0, 5);
			Ce(1, 0) = gpK(0, 6);
			Ce(2, 0) = gpK(0, 7);
			Ce(2, 1) = gpK(0, 8);

			B.resize(Ce.nrows, 2 * size);
			B.fill(0);
			distribute3x2(B.vals, dND.vals, dND.nrows, dND.ncols);

			CeB.multiply(Ce, B);
			Ke.multiply(B, CeB, detJ * weighFactor[gp] * gpThickness(0, 0), 1, true);
			break;

		case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::AXISYMMETRIC:

			Ce.resize(4, 4);
			Ce(0, 0) = gpK(0,  0);
			Ce(1, 1) = gpK(0,  1);
			Ce(2, 2) = gpK(0,  2);
			Ce(3, 3) = gpK(0,  3);
			Ce(0, 1) = gpK(0,  4);
			Ce(0, 2) = gpK(0,  5);
			Ce(0, 3) = gpK(0,  6);
			Ce(1, 2) = gpK(0,  7);
			Ce(1, 3) = gpK(0,  8);
			Ce(2, 3) = gpK(0,  9);
			Ce(1, 0) = gpK(0, 10);
			Ce(2, 0) = gpK(0, 11);
			Ce(2, 1) = gpK(0, 12);
			Ce(3, 0) = gpK(0, 13);
			Ce(3, 1) = gpK(0, 14);
			Ce(3, 2) = gpK(0, 15);

			B.resize(Ce.nrows, 2 * size);
			B.fill(0);
			distribute4x2(B.vals, dND.vals, dND.nrows, dND.ncols);
			for(esint i = 0; i < N[gp].ncols; i++) {
				B(2, i) = N[gp](0, i) / XY(0, 0);
			}

			CeB.multiply(Ce, B);
			Ke.multiply(B, CeB, detJ * weighFactor[gp] * 2 * M_PI * XY(0, 0), 1, true);
			break;
		}
	}

	double scale = iterator.minDesignVariable + pow(iterator.designVariable.data[0], iterator.penaltyFactor) * (1 - iterator.minDesignVariable);
	Ku.multiply(Ke, disp);
	double uku = Ku.dot(&disp);
	iterator.complianceDerivation.data[0] = -iterator.penaltyFactor * (1 - iterator.minDesignVariable) * pow(iterator.designVariable.data[0], iterator.penaltyFactor - 1) * uku;
	iterator.compliance.data[0] = scale * uku;
}





