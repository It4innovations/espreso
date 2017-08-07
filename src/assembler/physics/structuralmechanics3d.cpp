
#include "../../configuration/physics/structuralmechanics3d.h"
#include "structuralmechanics3d.h"

#include "../step.h"
#include "../instance.h"
#include "../constraints/equalityconstraints.h"

#include "../../mesh/settings/property.h"
#include "../../mesh/settings/evaluator.h"
#include "../../mesh/elements/element.h"
#include "../../mesh/structures/mesh.h"
#include "../../mesh/structures/material.h"
#include "../../mesh/structures/coordinates.h"
#include "../../mesh/structures/elementtypes.h"

#include "../../basis/matrices/denseMatrix.h"
#include "../../solver/specific/sparsesolvers.h"

using namespace espreso;

StructuralMechanics3D::StructuralMechanics3D(Mesh *mesh, Instance *instance, const StructuralMechanics3DConfiguration &configuration)
: Physics("STRUCTURAL MECHANICS 3D", mesh, instance), StructuralMechanics(configuration), _configuration(configuration)
{
	_equalityConstraints = new EqualityConstraints(*_instance, *_mesh, _mesh->nodes(), _mesh->faces(), pointDOFs(), pointDOFsOffsets());
}

void StructuralMechanics3D::prepare()
{
	_mesh->loadProperty(_configuration.displacement    , { "X", "Y", "Z" }     , { Property::DISPLACEMENT_X, Property::DISPLACEMENT_Y, Property::DISPLACEMENT_Z });
	_mesh->loadProperty(_configuration.acceleration    , { "X", "Y", "Z" }     , { Property::ACCELERATION_X, Property::ACCELERATION_Y, Property::ACCELERATION_Z });

	_mesh->loadMaterials(_configuration.materials, _configuration.material_set);

	for (size_t s = 1; s <= _configuration.physics_solver.load_steps; s++) {
		if (_configuration.physics_solver.load_steps_settings.find(s) != _configuration.physics_solver.load_steps_settings.end() &&
			_configuration.physics_solver.load_steps_settings.find(s)->second->espreso.regularization == REGULARIZATION::FIX_POINTS) {
			_mesh->computeFixPoints(8);
			break;
		}
	}

	StructuralMechanics::prepare();
}


void StructuralMechanics3D::analyticRegularization(size_t domain)
{
	if (_instance->K[domain].mtype != MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE) {
		ESINFO(ERROR) << "Cannot compute analytic regularization of not REAL_SYMMETRIC_POSITIVE_DEFINITE matrix. Set REGULARIZATION = NULL_PIVOTS";
	}

	ESTEST(MANDATORY) << "Too few FIX POINTS: " << _mesh->fixPoints(domain).size() << (_mesh->fixPoints(domain).size() > 3 ? TEST_PASSED : TEST_FAILED);

	size_t nodes = _mesh->coordinates().localSize(domain);
	_instance->N1[domain].rows = 3 * nodes;
	_instance->N1[domain].cols = 6;
	_instance->N1[domain].nnz = _instance->N1[domain].rows * _instance->N1[domain].cols;
	_instance->N1[domain].type = 'G';

	_instance->N1[domain].dense_values.reserve(_instance->N1[domain].nnz);

	for (size_t c = 0; c < 3; c++) {
		std::vector<double> kernel = { 0, 0, 0 };
		kernel[c] = 1;
		for (size_t i = 0; i < nodes; i++) {
			_instance->N1[domain].dense_values.insert(_instance->N1[domain].dense_values.end(), kernel.begin(), kernel.end());
		}
	}

	for (size_t i = 0; i < _mesh->coordinates().localSize(domain); i++) {
		const Point &p = _mesh->coordinates().get(i, domain);
		_instance->N1[domain].dense_values.push_back(-p.y);
		_instance->N1[domain].dense_values.push_back( p.x);
		_instance->N1[domain].dense_values.push_back(   0);
	}

	for (size_t i = 0; i < _mesh->coordinates().localSize(domain); i++) {
		const Point &p = _mesh->coordinates().get(i, domain);
		_instance->N1[domain].dense_values.push_back(-p.z);
		_instance->N1[domain].dense_values.push_back(   0);
		_instance->N1[domain].dense_values.push_back( p.x);
	}

	for (size_t i = 0; i < _mesh->coordinates().localSize(domain); i++) {
		const Point &p = _mesh->coordinates().get(i, domain);
		_instance->N1[domain].dense_values.push_back(   0);
		_instance->N1[domain].dense_values.push_back(-p.z);
		_instance->N1[domain].dense_values.push_back( p.y);
	}


	SparseMatrix Nt; // CSR matice s DOFY
	Nt.rows = 6;
	Nt.cols = _instance->K[domain].cols;
	Nt.nnz  = 9 * _mesh->fixPoints(domain).size();
	Nt.type = 'G';

	std::vector<eslocal> &ROWS = Nt.CSR_I_row_indices;
	std::vector<eslocal> &COLS = Nt.CSR_J_col_indices;
	std::vector<double>  &VALS = Nt.CSR_V_values;

	ROWS.reserve(Nt.rows + 1);
	COLS.reserve(Nt.nnz);
	VALS.reserve(Nt.nnz);

	ROWS.push_back(1);
	ROWS.push_back(ROWS.back() + _mesh->fixPoints(domain).size());
	ROWS.push_back(ROWS.back() + _mesh->fixPoints(domain).size());
	ROWS.push_back(ROWS.back() + _mesh->fixPoints(domain).size());
	ROWS.push_back(ROWS.back() + 2 * _mesh->fixPoints(domain).size());
	ROWS.push_back(ROWS.back() + 2 * _mesh->fixPoints(domain).size());
	ROWS.push_back(ROWS.back() + 2 * _mesh->fixPoints(domain).size());

	for (size_t c = 0; c < 3; c++) {
		for (size_t i = 0; i < _mesh->fixPoints(domain).size(); i++) {
			COLS.push_back(_mesh->fixPoints(domain)[i]->DOFIndex(domain, c) + 1);
		}
	}
	VALS.insert(VALS.end(), 3 * _mesh->fixPoints(domain).size(), 1);

	for (size_t i = 0; i < _mesh->fixPoints(domain).size(); i++) {
		const Point &p = _mesh->coordinates()[_mesh->fixPoints(domain)[i]->node(0)];
		COLS.push_back(_mesh->fixPoints(domain)[i]->DOFIndex(domain, 0) + 1);
		COLS.push_back(_mesh->fixPoints(domain)[i]->DOFIndex(domain, 1) + 1);
		VALS.push_back(-p.y);
		VALS.push_back( p.x);
	}

	for (size_t i = 0; i < _mesh->fixPoints(domain).size(); i++) {
		const Point &p = _mesh->coordinates()[_mesh->fixPoints(domain)[i]->node(0)];
		COLS.push_back(_mesh->fixPoints(domain)[i]->DOFIndex(domain, 0) + 1);
		COLS.push_back(_mesh->fixPoints(domain)[i]->DOFIndex(domain, 2) + 1);
		VALS.push_back(-p.z);
		VALS.push_back( p.x);
	}

	for (size_t i = 0; i < _mesh->fixPoints(domain).size(); i++) {
		const Point &p = _mesh->coordinates()[_mesh->fixPoints(domain)[i]->node(0)];
		COLS.push_back(_mesh->fixPoints(domain)[i]->DOFIndex(domain, 1) + 1);
		COLS.push_back(_mesh->fixPoints(domain)[i]->DOFIndex(domain, 2) + 1);
		VALS.push_back(-p.z);
		VALS.push_back( p.y);
	}

	SparseMatrix N;
	Nt.MatTranspose( N );

	_instance->RegMat[domain].MatMat(Nt, 'N', N);
	_instance->RegMat[domain].MatTranspose();
	_instance->RegMat[domain].RemoveLower();

	SparseSolverCPU NtN;
	NtN.ImportMatrix(_instance->RegMat[domain]);
	_instance->RegMat[domain].Clear();

	NtN.Factorization("Create RegMat");
	NtN.SolveMat_Sparse(Nt);
	NtN.Clear();

	_instance->RegMat[domain].MatMat(N, 'N', Nt);
	_instance->RegMat[domain].MatScale(_instance->K[domain].getDiagonalMaximum());
}

std::vector<std::pair<ElementType, Property> > StructuralMechanics3D::propertiesToStore() const
{
	return {};
}


void StructuralMechanics3D::assembleMaterialMatrix(const Step &step, const Element *e, eslocal node, double temp, DenseMatrix &K) const
{
	const Material* material = _mesh->materials()[e->param(Element::MATERIAL)];
	double Ex, Ey, Ez, miXY, miXZ, miYZ, Gx, Gy, Gz;

	switch (material->getModel(PHYSICS::STRUCTURAL_MECHANICS_2D)) {

	case MATERIAL_MODEL::LINEAR_ELASTIC_ISOTROPIC: {
		Ex = Ey = Ez = material->get(MATERIAL_PARAMETER::YOUNG_MODULUS_X)->evaluate(e->node(node), step.currentTime, temp);
		miXY = miXZ = miYZ = material->get(MATERIAL_PARAMETER::POISSON_RATIO_XY)->evaluate(e->node(node), step.currentTime, temp);

		double EE = Ex / ((1 + miXY) * (1 - 2 * miXY));

		K(node,  0) = EE * (1.0 - miXY);
		K(node,  1) = EE * (1.0 - miXY);
		K(node,  2) = EE * (1.0 - miXY);
		K(node,  3) = EE * (0.5 - miXY);
		K(node,  4) = EE * (0.5 - miXY);
		K(node,  5) = EE * (0.5 - miXY);

		K(node,  6) = EE * miXY;
		K(node,  7) = EE * miXY;
		K(node,  8) = 0;
		K(node,  9) = 0;
		K(node, 10) = 0;
		K(node, 11) = EE * miXY;
		K(node, 12) = 0;
		K(node, 13) = 0;
		K(node, 14) = 0;
		K(node, 15) = 0;
		K(node, 16) = 0;
		K(node, 17) = 0;
		K(node, 18) = 0;
		K(node, 19) = 0;
		K(node, 20) = 0;
		K(node, 21) = EE * miXY;
		K(node, 22) = EE * miXY;
		K(node, 23) = EE * miXY;
		K(node, 24) = 0;
		K(node, 25) = 0;
		K(node, 26) = 0;
		K(node, 27) = 0;
		K(node, 28) = 0;
		K(node, 29) = 0;
		K(node, 30) = 0;
		K(node, 31) = 0;
		K(node, 32) = 0;
		K(node, 33) = 0;
		K(node, 34) = 0;
		K(node, 35) = 0;
	} break;

	case MATERIAL_MODEL::LINEAR_ELASTIC_ANISOTROPIC: {
		K(node,  0) = material->get(MATERIAL_PARAMETER::D11)->evaluate(e->node(node), step.currentTime, temp);
		K(node,  1) = material->get(MATERIAL_PARAMETER::D22)->evaluate(e->node(node), step.currentTime, temp);
		K(node,  2) = material->get(MATERIAL_PARAMETER::D33)->evaluate(e->node(node), step.currentTime, temp);
		K(node,  3) = material->get(MATERIAL_PARAMETER::D44)->evaluate(e->node(node), step.currentTime, temp);
		K(node,  4) = material->get(MATERIAL_PARAMETER::D55)->evaluate(e->node(node), step.currentTime, temp);
		K(node,  5) = material->get(MATERIAL_PARAMETER::D66)->evaluate(e->node(node), step.currentTime, temp);

		K(node,  6) = material->get(MATERIAL_PARAMETER::D12)->evaluate(e->node(node), step.currentTime, temp);
		K(node,  7) = material->get(MATERIAL_PARAMETER::D13)->evaluate(e->node(node), step.currentTime, temp);
		K(node,  8) = material->get(MATERIAL_PARAMETER::D14)->evaluate(e->node(node), step.currentTime, temp);
		K(node,  9) = material->get(MATERIAL_PARAMETER::D15)->evaluate(e->node(node), step.currentTime, temp);
		K(node, 10) = material->get(MATERIAL_PARAMETER::D16)->evaluate(e->node(node), step.currentTime, temp);
		K(node, 11) = material->get(MATERIAL_PARAMETER::D23)->evaluate(e->node(node), step.currentTime, temp);
		K(node, 12) = material->get(MATERIAL_PARAMETER::D24)->evaluate(e->node(node), step.currentTime, temp);
		K(node, 13) = material->get(MATERIAL_PARAMETER::D25)->evaluate(e->node(node), step.currentTime, temp);
		K(node, 14) = material->get(MATERIAL_PARAMETER::D26)->evaluate(e->node(node), step.currentTime, temp);
		K(node, 15) = material->get(MATERIAL_PARAMETER::D34)->evaluate(e->node(node), step.currentTime, temp);
		K(node, 16) = material->get(MATERIAL_PARAMETER::D35)->evaluate(e->node(node), step.currentTime, temp);
		K(node, 17) = material->get(MATERIAL_PARAMETER::D36)->evaluate(e->node(node), step.currentTime, temp);
		K(node, 18) = material->get(MATERIAL_PARAMETER::D45)->evaluate(e->node(node), step.currentTime, temp);
		K(node, 19) = material->get(MATERIAL_PARAMETER::D46)->evaluate(e->node(node), step.currentTime, temp);
		K(node, 20) = material->get(MATERIAL_PARAMETER::D56)->evaluate(e->node(node), step.currentTime, temp);

		K(node, 21) = K(node,  6);
		K(node, 22) = K(node,  7);
		K(node, 23) = K(node, 11);
		K(node, 24) = K(node,  8);
		K(node, 25) = K(node, 12);
		K(node, 26) = K(node, 15);
		K(node, 27) = K(node,  9);
		K(node, 28) = K(node, 13);
		K(node, 29) = K(node, 16);
		K(node, 30) = K(node, 18);
		K(node, 31) = K(node, 10);
		K(node, 32) = K(node, 14);
		K(node, 33) = K(node, 17);
		K(node, 34) = K(node, 19);
		K(node, 35) = K(node, 20);
	} break;

	case MATERIAL_MODEL::LINEAR_ELASTIC_ORTHOTROPIC: {
		Ex = material->get(MATERIAL_PARAMETER::YOUNG_MODULUS_X)->evaluate(e->node(node), step.currentTime, temp);
		Ey = material->get(MATERIAL_PARAMETER::YOUNG_MODULUS_Y)->evaluate(e->node(node), step.currentTime, temp);
		Ez = material->get(MATERIAL_PARAMETER::YOUNG_MODULUS_Z)->evaluate(e->node(node), step.currentTime, temp);

		miXY = material->get(MATERIAL_PARAMETER::POISSON_RATIO_XY)->evaluate(e->node(node), step.currentTime, temp);
		miXZ = material->get(MATERIAL_PARAMETER::POISSON_RATIO_XZ)->evaluate(e->node(node), step.currentTime, temp);
		miYZ = material->get(MATERIAL_PARAMETER::POISSON_RATIO_YZ)->evaluate(e->node(node), step.currentTime, temp);

		Gx = material->get(MATERIAL_PARAMETER::SHEAR_MODULUS_XY)->evaluate(e->node(node), step.currentTime, temp);
		Gy = material->get(MATERIAL_PARAMETER::SHEAR_MODULUS_XZ)->evaluate(e->node(node), step.currentTime, temp);
		Gz = material->get(MATERIAL_PARAMETER::SHEAR_MODULUS_YZ)->evaluate(e->node(node), step.currentTime, temp);

		double miYX = miXY * Ey / Ex;
		double miZY = miYZ * Ez / Ey;
		double miZX = miXZ * Ex / Ez;

		double ksi = 1 - (miXY * miYX + miYZ * miZY + miXZ * miZX) - (miXY * miYZ * miZX + miYX * miZY * miXZ);

		double dxx = Ex * (1 - miYZ * miZY) / ksi;
		double dxy = Ey * (miXY + miXZ * miZY) / ksi;
		double dxz = Ez * (miXZ + miYZ * miXY)  /ksi;
		double dyy = Ey * (1 - miXZ * miZX) / ksi;
		double dyz = Ez * (miYZ + miYX * miXZ) / ksi;
		double dzz = Ez * (1 - miYX * miXY) / ksi;

		K(node,  0) = dxx;
		K(node,  1) = dyy;
		K(node,  2) = dzz;
		K(node,  3) = Gx;
		K(node,  4) = Gz;
		K(node,  5) = Gy;

		K(node,  6) = dxy;
		K(node,  7) = dxz;
		K(node,  8) = 0;
		K(node,  9) = 0;
		K(node, 10) = 0;
		K(node, 11) = dyz;
		K(node, 12) = 0;
		K(node, 13) = 0;
		K(node, 14) = 0;
		K(node, 15) = 0;
		K(node, 16) = 0;
		K(node, 17) = 0;
		K(node, 18) = 0;
		K(node, 19) = 0;
		K(node, 20) = 0;
		K(node, 21) = dxy;
		K(node, 22) = dxz;
		K(node, 23) = dyz;
		K(node, 24) = 0;
		K(node, 25) = 0;
		K(node, 26) = 0;
		K(node, 27) = 0;
		K(node, 28) = 0;
		K(node, 29) = 0;
		K(node, 30) = 0;
		K(node, 31) = 0;
		K(node, 32) = 0;
		K(node, 33) = 0;
		K(node, 34) = 0;
		K(node, 35) = 0;
	} break;

	default:
		ESINFO(ERROR) << "Structural mechanics 3D not supports set material model";
	}
}

void StructuralMechanics3D::processElement(const Step &step, Matrices matrices, const Element *e, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe, const std::vector<Solution*> &solution) const
{
	DenseMatrix Ce(6, 6), coordinates(e->nodes(), 3), J, invJ(3, 3), dND, B, epsilon, rhsT;
	DenseMatrix K(e->nodes(), 36), TE(e->nodes(), 3), inertia(e->nodes(), 3), dens(e->nodes(), 1);
	DenseMatrix gpK(e->nodes(), 36), gpTE(1, 3), gpInertia(1, 3), gpDens(1, 1);
	double detJ, temp, initTemp, CP = 1;

	const Material* material = _mesh->materials()[e->param(Element::MATERIAL)];

	for (size_t i = 0; i < e->nodes(); i++) {
		initTemp = e->getProperty(Property::INITIAL_TEMPERATURE, i, step.step, step.currentTime, 0, 0);
		temp = e->getProperty(Property::TEMPERATURE, i, step.step, step.currentTime, 0, initTemp);
		inertia(i, 0) = e->sumProperty(Property::ACCELERATION_X, i, step.step, step.currentTime, temp, 0);
		inertia(i, 1) = e->sumProperty(Property::ACCELERATION_Y, i, step.step, step.currentTime, temp, 0);
		inertia(i, 2) = e->sumProperty(Property::ACCELERATION_Z, i, step.step, step.currentTime, temp, 0);
		coordinates(i, 0) = _mesh->coordinates()[e->node(i)].x;
		coordinates(i, 1) = _mesh->coordinates()[e->node(i)].y;
		coordinates(i, 2) = _mesh->coordinates()[e->node(i)].z;
		dens(i, 0) = material->get(MATERIAL_PARAMETER::DENSITY)->evaluate(e->node(i), step.currentTime, temp);
		TE(i, 0) = (temp - initTemp) * material->get(MATERIAL_PARAMETER::THERMAL_EXPANSION_X)->evaluate(e->node(i), step.currentTime, temp);
		TE(i, 1) = (temp - initTemp) * material->get(MATERIAL_PARAMETER::THERMAL_EXPANSION_Y)->evaluate(e->node(i), step.currentTime, temp);
		TE(i, 2) = (temp - initTemp) * material->get(MATERIAL_PARAMETER::THERMAL_EXPANSION_Z)->evaluate(e->node(i), step.currentTime, temp);
		assembleMaterialMatrix(step, e, i, temp, K);
	}

	eslocal Ksize = pointDOFs().size() * e->nodes();

	Ke.resize(0, 0);
	Me.resize(0, 0);
	Re.resize(0, 0);
	fe.resize(0, 0);
	if (matrices & (Matrices::K | Matrices::R)) {
		Ke.resize(Ksize, Ksize);
		Ke = 0;
	}
	if (matrices & Matrices::M) {
		Me.resize(Ksize, Ksize);
		Me = 0;
	}
	if (matrices & Matrices::R) {
		Re.resize(Ksize, 1);
		Re = 0;
	}
	if (matrices & Matrices::f) {
		fe.resize(Ksize, 1);
		fe = 0;
	}

	for (size_t gp = 0; gp < e->gaussePoints(); gp++) {
		J.multiply(e->dN()[gp], coordinates);
		detJ = determinant3x3(J.values());
		inverse3x3(J.values(), invJ.values(), detJ);

		gpK.multiply(e->N()[gp], K);
		dND.multiply(invJ, e->dN()[gp]);
		gpDens.multiply(e->N()[gp], dens);

		if (matrices & Matrices::f) {
			gpTE.multiply(e->N()[gp], TE);
			gpInertia.multiply(e->N()[gp], inertia);
		}

		if (matrices & Matrices::M) {
			Me.multiply(e->N()[gp], e->N()[gp], gpDens(0, 0) * detJ * e->weighFactor()[gp] * CP, 1, true);
		}

		Ce.resize(6, 6);
		size_t k = 0;
		for (size_t i = 0; i < 6; i++) {
			Ce(i, i) = gpK(0, k++);
		}
		for (size_t i = 0; i < 6; i++) {
			for (size_t j = i + 1; j < 6; j++) {
				Ce(i, j) = gpK(0, k++);
			}
		}
		for (size_t i = 0; i < 6; i++) {
			for (size_t j = 0; j < i; j++) {
				Ce(i, j) = gpK(0, k++);
			}
		}

		B.resize(Ce.rows(), Ksize);
		distribute6x3(B.values(), dND.values(), dND.rows(), dND.columns());

		if (matrices & Matrices::K) {
			Ke.multiply(B, Ce * B, detJ * e->weighFactor()[gp], 1, true);
		}

		if (matrices & Matrices::f) {
			epsilon.resize(Ce.rows(), 1);
			epsilon(0, 0) = gpTE(0, 0);
			epsilon(1, 0) = gpTE(0, 1);
			epsilon(2, 0) = gpTE(0, 2);
			epsilon(3, 0) = epsilon(4, 0) = epsilon(5, 0) = 0;

			rhsT.multiply(B, Ce * epsilon, detJ * e->weighFactor()[gp], 0, true, false);
			for (eslocal i = 0; i < Ksize; i++) {
				fe(i, 0) += gpDens(0, 0) * detJ * e->weighFactor()[gp] * e->N()[gp](0, i % e->nodes()) * gpInertia(0, i / e->nodes());
				fe(i, 0) += rhsT(i, 0);
			}
		}
	}
}

void StructuralMechanics3D::processFace(const Step &step, Matrices matrices, const Element *e, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe, const std::vector<Solution*> &solution) const
{
	if (!e->hasProperty(Property::PRESSURE, step.step)) {
		Ke.resize(0, 0);
		Me.resize(0, 0);
		Re.resize(0, 0);
		fe.resize(0, 0);
		return;
	}
	if (!(matrices & (Matrices::K | Matrices::f))) {
		Ke.resize(0, 0);
		Me.resize(0, 0);
		Re.resize(0, 0);
		fe.resize(0, 0);
		return;
	}

	DenseMatrix coordinates(e->nodes(), 3), dND(1, 3), P(e->nodes(), 1), normal(1, 3);
	DenseMatrix gpP(1, 1), gpQ(1, 3);

	eslocal Ksize = pointDOFs().size() * e->nodes();
	Ke.resize(0, 0);
	Me.resize(0, 0);
	Re.resize(0, 0);
	fe.resize(0, 0);

	if (matrices & Matrices::f) {
		fe.resize(Ksize, 1);
		fe = 0;
	}

	for (size_t n = 0; n < e->nodes(); n++) {
		coordinates(n, 0) = _mesh->coordinates()[e->node(n)].x;
		coordinates(n, 1) = _mesh->coordinates()[e->node(n)].y;
		coordinates(n, 2) = _mesh->coordinates()[e->node(n)].z;
		P(n, 0) = e->getProperty(Property::PRESSURE, n, step.step, step.currentTime, 0, 0);
	}

	for (size_t gp = 0; gp < e->gaussePoints(); gp++) {
		dND.multiply(e->dN()[gp], coordinates);
		Point v2(dND(0, 0), dND(0, 1), dND(0, 2));
		Point v1(dND(1, 0), dND(1, 1), dND(1, 2));
		Point va = Point::cross(v1, v2);
		e->rotateOutside(e->parentElements()[0], _mesh->coordinates(), va);
		double J = va.norm();
		normal(0, 0) = va.x / va.norm();
		normal(0, 1) = va.y / va.norm();
		normal(0, 2) = va.z / va.norm();

		gpP.multiply(e->N()[gp], P);
		gpQ.multiply(normal, gpP, 1, 0, true);

		for (eslocal i = 0; i < Ksize; i++) {
			fe(i, 0) += J * e->weighFactor()[gp] * e->N()[gp](0, i % e->nodes()) * gpQ(0, i / e->nodes());
		}
	}
}

void StructuralMechanics3D::processEdge(const Step &step, Matrices matrices, const Element *e, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe, const std::vector<Solution*> &solution) const
{
	Ke.resize(0, 0);
	Me.resize(0, 0);
	Re.resize(0, 0);
	fe.resize(0, 0);
}

void StructuralMechanics3D::processNode(const Step &step, Matrices matrices, const Element *e, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe, const std::vector<Solution*> &solution) const
{
	if (
			e->hasProperty(Property::FORCE_X, step.step) ||
			e->hasProperty(Property::FORCE_Y, step.step) ||
			e->hasProperty(Property::FORCE_Z, step.step)) {

		Ke.resize(0, 0);
		Me.resize(0, 0);
		Re.resize(0, 0);
		fe.resize(pointDOFs().size(), 0);

		fe(0, 0) = e->sumProperty(Property::FORCE_X, 0, step.step, step.currentTime, 0, 0);
		fe(1, 0) = e->sumProperty(Property::FORCE_Y, 0, step.step, step.currentTime, 0, 0);
		fe(2, 0) = e->sumProperty(Property::FORCE_Z, 0, step.step, step.currentTime, 0, 0);
		return;
	}
	Ke.resize(0, 0);
	Me.resize(0, 0);
	Re.resize(0, 0);
	fe.resize(0, 0);
}

void StructuralMechanics3D::postProcessElement(const Step &step, const Element *e, std::vector<Solution*> &solution)
{

}

void StructuralMechanics3D::processSolution(const Step &step)
{
	if (!_configuration.post_process) {
		return;
	}
}





