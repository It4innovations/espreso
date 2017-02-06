
#include "../../old_physics/elasticity2d/assembler.h"

#include "../../../basis/matrices/denseMatrix.h"
#include "../../../basis/matrices/sparseVVPMatrix.h"
#include "../../../basis/matrices/sparseCSRMatrix.h"
#include "../../../configuration/physics/linearelasticity2d.h"
#include "../../../solver/generic/SparseMatrix.h"
#include "../../../solver/specific/sparsesolvers.h"

#include "../../../mesh/elements/element.h"
#include "../../../mesh/settings/evaluator.h"
#include "../../../mesh/structures/mesh.h"
#include "../../../mesh/structures/material.h"

#include "../../../mesh/elements/plane/square4.h"
#include "../../../mesh/elements/plane/square8.h"
#include "../../../mesh/elements/plane/triangle3.h"
#include "../../../mesh/elements/plane/triangle6.h"

#include "../../constraints/equalityconstraints.h"

#include "../../../output/vtk/vtk.h"

namespace espreso {

std::vector<Property> Elasticity2D::elementDOFs;
std::vector<Property> Elasticity2D::faceDOFs;
std::vector<Property> Elasticity2D::edgeDOFs;
std::vector<Property> Elasticity2D::pointDOFs = { Property::DISPLACEMENT_X, Property::DISPLACEMENT_Y };
std::vector<Property> Elasticity2D::midPointDOFs = { Property::DISPLACEMENT_X, Property::DISPLACEMENT_Y };

Elasticity2D::Elasticity2D(Mesh &mesh, Constraints &constraints, const LinearElasticity2DConfiguration &configuration)
: OldPhysics(
		mesh, constraints, configuration.espreso,
		MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE,
		elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs),
  _configuration(configuration) {};

void Elasticity2D::assembleStiffnessMatrices()
{
	ESINFO(PROGRESS2) << "Assemble matrices K, kernels, and RHS.";
	#pragma omp parallel for
	for (size_t p = 0; p < _mesh.parts(); p++) {
		composeSubdomain(p);
		K[p].mtype = mtype;
		ESINFO(PROGRESS2) << Info::plain() << ".";
	}
	ESINFO(PROGRESS2);
}

void Elasticity2D::prepareMeshStructures()
{
	Square4::setDOFs(elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs);
	Square8::setDOFs(elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs);
	Triangle3::setDOFs(elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs);
	Triangle6::setDOFs(elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs);

	matrixSize = _mesh.assignUniformDOFsIndicesToNodes(matrixSize, pointDOFs);
	_mesh.computeNodesDOFsCounters(pointDOFs);

	if (_solverConfiguration.regularization == REGULARIZATION::FIX_POINTS) {
		_mesh.computeFixPoints(4);
	}

	if (_solverConfiguration.method == ESPRESO_METHOD::HYBRID_FETI) {
		switch (_solverConfiguration.B0_type) {
		case B0_TYPE::CORNERS:
			_mesh.computePlaneCorners(1, true, false);
			break;
		case B0_TYPE::KERNELS:
			_mesh.computeEdgesSharedByDomains();
			break;
		default:
			break;
		}
	}

	_constraints.initMatrices(matrixSize);

	_mesh.loadNodeProperty(_configuration.displacement       , { "X", "Y" }, { Property::DISPLACEMENT_X, Property::DISPLACEMENT_Y });
	_mesh.loadNodeProperty(_configuration.temperature        , { }         , { Property::TEMPERATURE });
	_mesh.loadNodeProperty(_configuration.obstacle           , { }         , { Property::OBSTACLE });
	_mesh.loadNodeProperty(_configuration.normal_direction   , { }         , { Property::NORMAL_DIRECTION });
	_mesh.loadNodeProperty(_configuration.thickness          , { }         , { Property::THICKNESS });

	_mesh.loadProperty(_configuration.normal_presure     , { }         , { Property::PRESSURE });
	_mesh.loadProperty(_configuration.acceleration       , { "X", "Y" }, { Property::ACCELERATION_X, Property::ACCELERATION_Y });
	_mesh.loadProperty(_configuration.initial_temperature, { }         , { Property::INITIAL_TEMPERATURE });

	_mesh.loadMaterials(_configuration.materials, _configuration.material_set);
	_mesh.removeDuplicateRegions();
}

void Elasticity2D::saveMeshProperties(store::ResultStore &store)
{
	store.storeProperty("displacement", { Property::DISPLACEMENT_X, Property::DISPLACEMENT_Y }, store::ResultStore::ElementType::NODES);
	store.storeProperty("forces", { Property::FORCE_X, Property::FORCE_Y }, store::ResultStore::ElementType::NODES);
	store.storeProperty("obstacle", { Property::OBSTACLE }, store::ResultStore::ElementType::NODES);
	store.storeProperty("normal_direction", { Property::NORMAL_DIRECTION }, store::ResultStore::ElementType::NODES);
	if (_solverConfiguration.regularization == REGULARIZATION::FIX_POINTS) {
		store::VTK::fixPoints(store.configuration(), _mesh, "fixPoints");
	}
	if (_solverConfiguration.method == ESPRESO_METHOD::HYBRID_FETI) {
		switch (_solverConfiguration.B0_type) {
		case B0_TYPE::CORNERS:
		case B0_TYPE::COMBINED:
			store::VTK::mesh(store.configuration(), _mesh, "edges", store::ResultStore::ElementType::EDGES);
			store::VTK::corners(store.configuration(), _mesh, "corners");
			break;
		case B0_TYPE::KERNELS:
			store::VTK::mesh(store.configuration(), _mesh, "edges", store::ResultStore::ElementType::EDGES);
			break;
		default:
			ESINFO(GLOBAL_ERROR) << "Not implemented saving properties of B0";
		}
	}
}

void Elasticity2D::saveMeshResults(store::ResultStore &store, const std::vector<std::vector<double> > &results)
{
	store.storeValues("displacement", 2, results, store::ResultStore::ElementType::NODES);
}

void Elasticity2D::assembleB1()
{
	EqualityConstraints::insertDirichletToB1(_constraints, _mesh.nodes(), pointDOFs);
	EqualityConstraints::insertElementGluingToB1(_constraints, _mesh.nodes(), pointDOFs, K);
}

void Elasticity2D::assembleB0()
{
	if (_solverConfiguration.method == ESPRESO_METHOD::HYBRID_FETI) {
		switch (_solverConfiguration.B0_type) {
		case B0_TYPE::CORNERS:
			EqualityConstraints::insertDomainGluingToB0(_constraints, _mesh.corners(), pointDOFs);
			break;
		case B0_TYPE::KERNELS:
			EqualityConstraints::insertKernelsToB0(_constraints, _mesh.edges(), pointDOFs, R1);
			break;
		default:
			break;
		}
	}
}

static double determinant2x2(DenseMatrix &m)
{
	return fabs(
		m(0, 0) * m(1, 1) - m(1, 0) * m(0, 1)
	);
}

static void inverse(const DenseMatrix &m, DenseMatrix &inv, double det)
{
	inv.resize(m.rows(), m.columns());
	double detJx = 1 / det;
	inv(0, 0) =   detJx * m(1, 1);
	inv(0, 1) = - detJx * m(0, 1);
	inv(1, 0) = - detJx * m(1, 0);
	inv(1, 1) =   detJx * m(0, 0);
}

// B =
// dX   0
//  0  dY
// dY  dX
static void distribute(DenseMatrix &B, DenseMatrix &dND)
{
	eslocal columns = dND.rows() * dND.columns();
	const double *dNDx = dND.values();
	const double *dNDy = dND.values() + dND.columns();

	double *v = B.values();

	memcpy(&v[0], dNDx,                               sizeof(double) * dND.columns());
	memcpy(&v[2 * columns + dND.columns()],     dNDx, sizeof(double) * dND.columns());

	memcpy(&v[1 * columns + dND.columns()],     dNDy, sizeof(double) * dND.columns());
	memcpy(&v[2 * columns],                     dNDy, sizeof(double) * dND.columns());
}


// B =
// dX   0
//  0  dY
// N/x  0
// dY  dX
static void distribute(DenseMatrix &B, DenseMatrix &dND, const DenseMatrix &N, DenseMatrix &XY)
{
	eslocal columns = dND.rows() * dND.columns();
	const double *dNDx = dND.values();
	const double *dNDy = dND.values() + dND.columns();

	double *v = B.values();

	memcpy(&v[0], dNDx,                               sizeof(double) * dND.columns());
	memcpy(&v[3 * columns + dND.columns()],     dNDx, sizeof(double) * dND.columns());

	memcpy(&v[1 * columns + dND.columns()],     dNDy, sizeof(double) * dND.columns());
	memcpy(&v[3 * columns],                     dNDy, sizeof(double) * dND.columns());

	for(size_t i = 0; i < N.columns(); i++) {
		B(2, i) = N(0, i) / XY(0, 0);
	}
}

static void fillC(DenseMatrix &C, LinearElasticity2DConfiguration::ELEMENT_BEHAVIOUR behaviour, MATERIAL_MODEL model, DenseMatrix &dens, DenseMatrix &E, DenseMatrix &mi, DenseMatrix &T)
{
	switch (model) {

	case MATERIAL_MODEL::LINEAR_ELASTIC_ISOTROPIC:
	{

		switch (behaviour) {

		case LinearElasticity2DConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRAIN:
		{
			double k = E(0, 0) * (1 - mi(0, 0)) / ((1 + mi(0, 0)) * (1 - 2 * mi(0, 0)));
			C.resize(3, 3);
			C(0, 0) = C(1, 1) = k * 1;
			C(1, 0) = C(0, 1) = k * (mi(0, 0) / (1 - mi(0, 0)));
			C(2, 2) = k * ((1 - 2 * mi(0, 0)) / (2 * (1 - mi(0, 0))));
			return;
		}

		case LinearElasticity2DConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRESS:
		case LinearElasticity2DConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRESS_WITH_THICKNESS:
		{
			double k = E(0, 0) / (1 - mi(0, 0) * mi(0, 0));
			C.resize(3, 3);
			C(0, 0) = C(1, 1) = k * 1;
			C(1, 0) = C(0, 1) = k * mi(0, 0);
			C(2, 2) = k * ((1 -  mi(0, 0)) / 2);
			return;
		}

		case LinearElasticity2DConfiguration::ELEMENT_BEHAVIOUR::AXISYMMETRIC:
		{
			double k = E(0, 0) * (1 - mi(0, 0)) / ((1 + mi(0, 0)) * (1 - 2 * mi(0, 0)));
			C.resize(4, 4);
			C(0, 0) = C(1, 1) = C(3, 3) = k * 1;
			C(1, 0) = C(0, 1) = C(3, 0) = C(0, 3) = C(1, 3) = C(3, 1) = k * (mi(0, 0) / (1 - mi(0, 0)));
			C(2, 2) = k * ((1 - 2 * mi(0, 0)) / (2 * (1 - mi(0, 0))));
			return;
		}
		}
	}

	case MATERIAL_MODEL::LINEAR_ELASTIC_ORTHOTROPIC:
	{
		return;
	}

	case MATERIAL_MODEL::LINEAR_ELASTIC_ANISOTROPIC:
	{
		return;
	}

	default:
		ESINFO(ERROR) << "Linear elasticity 2D not supports set material model";
	}
}

static void processElement(DenseMatrix &Ke, std::vector<double> &fe, const espreso::Mesh &mesh, const Element* element, const LinearElasticity2DConfiguration &configuration)
{
	DenseMatrix Ce(4, 4), XY(1, 2), coordinates, J, invJ, dND, B, epsilon, rhsT;
	DenseMatrix
			matDENS(element->nodes(), 1), matE(element->nodes(), 2), matMI(element->nodes(), 1),
			matTE(element->nodes(), 2), matT(element->nodes(), 1), matThickness(element->nodes(), 1),
			matInitT(element->nodes(), 1), inertia(element->nodes(), 2);
	DenseMatrix gpDENS(1, 1), gpE(1, 2), gpMI(1, 1), gpTE(1, 2), gpT(1, 1), gpInitT(1, 1), gpThickness(1, 1), gpInertia(1, 2);
	double detJ;

	const Material* material = mesh.materials()[element->param(Element::MATERIAL)];
	const std::vector<DenseMatrix> &dN = element->dN();
	const std::vector<DenseMatrix> &N = element->N();
	const std::vector<double> &weighFactor = element->weighFactor();

	coordinates.resize(element->nodes(), 2);

	inertia = 0;
	for (size_t i = 0; i < element->nodes(); i++) {
		switch (material->getModel(PHYSICS::LINEAR_ELASTICITY_2D)) {
		case MATERIAL_MODEL::LINEAR_ELASTIC_ISOTROPIC:
			matE(i, 0) = matE(i, 1) = material->get(MATERIAL_PARAMETER::YOUNG_MODULUS_X)->evaluate(element->node(i));
			matMI(i, 0) = material->get(MATERIAL_PARAMETER::POISSON_RATIO_XY)->evaluate(element->node(i));
			matTE(i, 0) = matTE(i, 1) = material->get(MATERIAL_PARAMETER::THERMAL_EXPANSION_X)->evaluate(element->node(i));
			break;
		case MATERIAL_MODEL::LINEAR_ELASTIC_ANISOTROPIC:
			ESINFO(ERROR) << "Implement ANISOTROPIC MATERIAL";
			break;
		case MATERIAL_MODEL::LINEAR_ELASTIC_ORTHOTROPIC:
			matE(i, 0) = material->get(MATERIAL_PARAMETER::YOUNG_MODULUS_X)->evaluate(element->node(i));
			matE(i, 1) = material->get(MATERIAL_PARAMETER::YOUNG_MODULUS_Y)->evaluate(element->node(i));

			matMI(i, 0) = material->get(MATERIAL_PARAMETER::POISSON_RATIO_XY)->evaluate(element->node(i));

			matTE(i, 0) = material->get(MATERIAL_PARAMETER::THERMAL_EXPANSION_X)->evaluate(element->node(i));
			matTE(i, 1) = material->get(MATERIAL_PARAMETER::THERMAL_EXPANSION_Y)->evaluate(element->node(i));
		default:
			ESINFO(ERROR) << "Linear elasticity 2D not supports set material model";
		}

		matInitT(i, 0) = element->getProperty(Property::INITIAL_TEMPERATURE, i, 0, 0);
		matT(i, 0) = element->getProperty(Property::TEMPERATURE, i, 0, matInitT(i, 0));
		inertia(i, 0) = element->sumProperty(Property::ACCELERATION_X, i, 0, 0);
		inertia(i, 1) = element->sumProperty(Property::ACCELERATION_Y, i, 0, 0);
		matThickness(i, 0) = element->getProperty(Property::THICKNESS, i, 0, 1);
		coordinates(i, 0) = mesh.coordinates()[element->node(i)].x;
		coordinates(i, 1) = mesh.coordinates()[element->node(i)].y;
	}

	eslocal Ksize = 2 * element->nodes();

	Ke.resize(Ksize, Ksize);
	Ke = 0;
	fe.resize(Ksize);
	std::fill(fe.begin(), fe.end(), 0);
	rhsT.resize(Ksize, 1);
	rhsT = 0;

	for (size_t gp = 0; gp < element->gaussePoints(); gp++) {
		J.multiply(dN[gp], coordinates);
		detJ = determinant2x2(J);
		inverse(J, invJ, detJ);
		dND.multiply(invJ, dN[gp]);
		XY.multiply(N[gp], coordinates);
		gpDENS.multiply(N[gp], matDENS);
		gpE.multiply(N[gp], matE);
		gpMI.multiply(N[gp], matMI);
		gpTE.multiply(N[gp], matTE);
		gpT.multiply(N[gp], matT);
		gpInitT.multiply(N[gp], matInitT);
		gpInertia.multiply(N[gp], inertia);

		gpThickness.multiply(N[gp], matThickness, 1, 0);

		fillC(Ce, configuration.element_behaviour, material->getModel(PHYSICS::LINEAR_ELASTICITY_2D), gpDENS, gpE, gpMI, gpT);
		B.resize(Ce.rows(), Ksize);
		epsilon.resize(Ce.rows(), 1);
		switch (configuration.element_behaviour) {

		case LinearElasticity2DConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRESS:
		case LinearElasticity2DConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRAIN:
			gpThickness(0, 0) = 1;
		case LinearElasticity2DConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRESS_WITH_THICKNESS:
			epsilon(0, 0) = (gpT(0, 0) - gpInitT(0, 0)) * gpTE(0, 0);
			epsilon(1, 0) = (gpT(0, 0) - gpInitT(0, 0)) * gpTE(0, 1);
			epsilon(2, 0) = 0;
			distribute(B, dND);
			Ke.multiply(B, Ce * B, detJ * weighFactor[gp] * gpThickness(0, 0), 1, true);
			rhsT.multiply(B, Ce * epsilon, detJ * weighFactor[gp] * gpThickness(0, 0), 0, true, false);
			for (eslocal i = 0; i < Ksize; i++) {
				fe[i] += gpDENS(0, 0) * detJ * weighFactor[gp] * gpThickness(0, 0) * N[gp](0, i % element->nodes()) * gpInertia(0, i / element->nodes());
				fe[i] += gpDENS(0, 0) * detJ * weighFactor[gp] * gpThickness(0, 0) * N[gp](0, i % element->nodes()) * XY(0, i / element->nodes()) * pow(configuration.angular_velocity_z, 2);
				fe[i] += rhsT(i, 0);
			}
			break;

		case LinearElasticity2DConfiguration::ELEMENT_BEHAVIOUR::AXISYMMETRIC:
			epsilon(0, 0) = (gpT(0, 0) - gpInitT(0, 0)) * gpTE(0, 0);
			epsilon(1, 0) = (gpT(0, 0) - gpInitT(0, 0)) * gpTE(0, 1);
			epsilon(2, 0) = epsilon(3, 0) = 0;
			distribute(B, dND, N[gp], XY);
			Ke.multiply(B, Ce * B, detJ * weighFactor[gp] * 2 * M_PI * XY(0, 0), 1, true);
			rhsT.multiply(B, Ce * epsilon, detJ * weighFactor[gp] * 2 * M_PI * XY(0, 0), 0, true, false);
			for (eslocal i = 0; i < Ksize; i++) {
				fe[i] += gpDENS(0, 0) * detJ * weighFactor[gp] * 2 * M_PI * XY(0, 0) * N[gp](0, i % element->nodes()) * gpInertia(0, i / element->nodes());
				fe[i] += rhsT(i, 0);
			}
			for (eslocal i = 0; i < Ksize / 2; i++) {
				fe[i] += gpDENS(0, 0) * detJ * weighFactor[gp] * 2 * M_PI * XY(0, 0) * N[gp](0, i % element->nodes()) * XY(0, 0) * pow(configuration.angular_velocity_y, 2);
				fe[Ksize / 2 + i] += gpDENS(0, 0) * detJ * weighFactor[gp] * 2 * M_PI * XY(0, 0) * N[gp](0, i % element->nodes()) * XY(0, 1) * pow(configuration.angular_velocity_x, 2);
			}
			break;
		}
	}
}

static void processEdge(std::vector<double> &fe, const espreso::Mesh &mesh, const Element* edge, const LinearElasticity2DConfiguration &configuration)
{
	DenseMatrix coordinates(edge->nodes(), 2), dND(1, 2), P(edge->nodes(), 1), normal(2, 1), matThickness(edge->nodes(), 1), XY(1, 2);
	DenseMatrix gpP(1, 1), gpQ(1, 2), gpThickness(1, 1);

	const std::vector<DenseMatrix> &dN = edge->dN();
	const std::vector<DenseMatrix> &N = edge->N();
	const std::vector<double> &weighFactor = edge->weighFactor();

	for (size_t n = 0; n < edge->nodes(); n++) {
		coordinates(n, 0) = mesh.coordinates()[edge->node(n)].x;
		coordinates(n, 1) = mesh.coordinates()[edge->node(n)].y;
		P(n, 0) = edge->getProperty(Property::PRESSURE, n, 0, 0);
		matThickness(n, 0) = edge->getProperty(Property::THICKNESS, n, 0, 1);
	}

	eslocal Ksize = 2 * edge->nodes();
	fe.resize(Ksize);
	std::fill(fe.begin(), fe.end(), 0);

	for (size_t gp = 0; gp < edge->gaussePoints(); gp++) {
		dND.multiply(dN[gp], coordinates);
		double J = dND.norm();
		Point n(-dND(0, 1), dND(0, 0), 0);
		edge->rotateOutside(edge->parentElements()[0], mesh.coordinates(), n);
		normal(0, 0) = n.x / J;
		normal(1, 0) = n.y / J;
		gpP.multiply(N[gp], P);
		gpQ.multiply(normal, gpP);
		gpThickness.multiply(N[gp], matThickness);

		switch (configuration.element_behaviour) {

		case LinearElasticity2DConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRESS:
		case LinearElasticity2DConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRAIN:
			gpThickness(0, 0) = 1;
		case LinearElasticity2DConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRESS_WITH_THICKNESS:
			for (eslocal i = 0; i < Ksize; i++) {
				fe[i] += gpThickness(0, 0) * J * weighFactor[gp] * N[gp](0, i % edge->nodes()) * gpQ(0, i / edge->nodes());
			}
			break;

		case LinearElasticity2DConfiguration::ELEMENT_BEHAVIOUR::AXISYMMETRIC:
			XY.multiply(N[gp], coordinates);
			for (eslocal i = 0; i < Ksize; i++) {
				fe[i] += gpThickness(0, 0) * J * weighFactor[gp] * 2 * M_PI * XY(0, 0) * N[gp](0, i % edge->nodes()) * gpQ(0, i / edge->nodes());
			}
			break;
		}
	}
}

static void analyticsKernels(SparseMatrix &R1, const Coordinates &coordinates, size_t subdomain)
{
	size_t nodes = coordinates.localSize(subdomain);
	R1.rows = 2 * nodes;
	R1.cols = 3;
	R1.nnz = R1.rows * R1.cols;
	R1.type = 'G';

	R1.dense_values.reserve(R1.nnz);

	for (size_t c = 0; c < 2; c++) {
		std::vector<double> kernel = { 0, 0 };
		kernel[c] = 1;
		for (size_t i = 0; i < nodes; i++) {
			R1.dense_values.insert(R1.dense_values.end(), kernel.begin(), kernel.end());
		}
	}

	for (size_t i = 0; i < coordinates.localSize(subdomain); i++) {
		const Point &p = coordinates.get(i, subdomain);
		R1.dense_values.push_back(-p.y);
		R1.dense_values.push_back( p.x);
	}
}

static void analyticsRegMat(SparseMatrix &K, SparseMatrix &RegMat, const std::vector<Element*> &fixPoints, const Coordinates &coordinates, size_t subdomain)
{
	ESTEST(MANDATORY) << "Too few FIX POINTS: " << fixPoints.size() << (fixPoints.size() > 3 ? TEST_PASSED : TEST_FAILED);

	SparseMatrix Nt; // CSR matice s DOFY
	Nt.rows = 3;
	Nt.cols = K.cols;
	Nt.nnz  = 4 * fixPoints.size();
	Nt.type = 'G';

	std::vector<eslocal> &ROWS = Nt.CSR_I_row_indices;
	std::vector<eslocal> &COLS = Nt.CSR_J_col_indices;
	std::vector<double>  &VALS = Nt.CSR_V_values;

	ROWS.reserve(Nt.rows + 1);
	COLS.reserve(Nt.nnz);
	VALS.reserve(Nt.nnz);

	ROWS.push_back(1);
	ROWS.push_back(ROWS.back() + fixPoints.size());
	ROWS.push_back(ROWS.back() + fixPoints.size());
	ROWS.push_back(ROWS.back() + 2 * fixPoints.size());

	for (size_t c = 0; c < 2; c++) {
		std::vector<double> kernel = { 0, 0 };

		kernel[c] = 1;
		for (size_t i = 0; i < fixPoints.size(); i++) {
			COLS.push_back(fixPoints[i]->DOFIndex(subdomain, 0) + 1);
			COLS.push_back(fixPoints[i]->DOFIndex(subdomain, 1) + 1);
			VALS.insert(VALS.end(), kernel.begin(), kernel.end());
		}
	}

	for (size_t i = 0; i < fixPoints.size(); i++) {
		const Point &p = coordinates[fixPoints[i]->node(0)];
		COLS.push_back(fixPoints[i]->DOFIndex(subdomain, 0) + 1);
		COLS.push_back(fixPoints[i]->DOFIndex(subdomain, 1) + 1);
		VALS.push_back(-p.y);
		VALS.push_back( p.x);
	}

	SparseMatrix N;
	Nt.MatTranspose( N );

	RegMat.MatMat(Nt, 'N', N);
	RegMat.MatTranspose();
	RegMat.RemoveLower();

	SparseSolverCPU NtN;
	NtN.ImportMatrix(RegMat);
	RegMat.Clear();

	NtN.Factorization("Create RegMat");
	NtN.SolveMat_Sparse(Nt);
	NtN.Clear();

	RegMat.MatMat(N, 'N', Nt);
	RegMat.MatScale(K.getDiagonalMaximum());
}

static void algebraicKernelsAndRegularization(SparseMatrix &K, SparseMatrix &RegMat, SparseMatrix &R, size_t subdomain)
{
	double norm;
	eslocal defect;

	K.get_kernel_from_K(K, RegMat, R, norm, defect, subdomain);
}

void Elasticity2D::assembleStiffnessMatrix(const Element* e, DenseMatrix &Ke, std::vector<double> &fe, std::vector<eslocal> &dofs) const
{
	processElement(Ke, fe, _mesh, e, _configuration);
	dofs.resize(e->nodes() * pointDOFs.size());
	for (size_t dof = 0, i = 0; dof < pointDOFs.size(); dof++) {
		for (size_t n = 0; n < e->nodes(); n++, i++) {
			dofs[i] = e->node(n) * pointDOFs.size() + dof;
		}
	}


	std::vector<Property> forces = { Property::FORCE_X, Property::FORCE_Y };
	for (size_t n = 0; n < e->nodes(); n++) {
		for (size_t dof = 0; dof < pointDOFs.size(); dof++) {
			fe[n * pointDOFs.size() + dof] = e->sumProperty(forces[dof], n, 0, 0) / _mesh.nodes()[e->node(n)]->domains().size();
		}
	}
}

void Elasticity2D::makeStiffnessMatricesRegular()
{
	for (size_t subdomain = 0; subdomain < K.size(); subdomain++) {
		switch (_solverConfiguration.regularization) {
		case REGULARIZATION::FIX_POINTS:
			analyticsKernels(R1[subdomain], _mesh.coordinates(), subdomain);
			analyticsRegMat(K[subdomain], RegMat[subdomain], _mesh.fixPoints(subdomain), _mesh.coordinates(), subdomain);
			K[subdomain].RemoveLower();
			RegMat[subdomain].RemoveLower();
			K[subdomain].MatAddInPlace(RegMat[subdomain], 'N', 1);
			RegMat[subdomain].ConvertToCOO(1);
			break;
		case REGULARIZATION::NULL_PIVOTS:
			K[subdomain].RemoveLower();
			algebraicKernelsAndRegularization(K[subdomain], RegMat[subdomain], R1[subdomain], subdomain);
			break;
		}
	}
}

void Elasticity2D::composeSubdomain(size_t subdomain)
{
	if (!_configuration.post_process) {
		return;
	}

	SparseVVPMatrix<eslocal> _K;
	DenseMatrix Ke;
	std::vector<double> fe;

	_K.resize(matrixSize[subdomain], matrixSize[subdomain]);
	f[subdomain].resize(matrixSize[subdomain]);

	const std::vector<eslocal> &partition = _mesh.getPartition();
	const std::vector<Element*> &elements = _mesh.elements();
	const std::vector<Element*> &nodes = _mesh.nodes();

	for (eslocal e = partition[subdomain]; e < partition[subdomain + 1]; e++) {

		processElement(Ke, fe, _mesh, elements[e], _configuration);

		for (size_t nx = 0; nx < elements[e]->nodes(); nx++) {
			for (size_t dx = 0; dx < pointDOFs.size(); dx++) {
				size_t row = nodes[elements[e]->node(nx)]->DOFIndex(subdomain, dx);
				for (size_t ny = 0; ny < elements[e]->nodes(); ny++) {
					for (size_t dy = 0; dy < pointDOFs.size(); dy++) {
						size_t column = nodes[elements[e]->node(ny)]->DOFIndex(subdomain, dy);
						_K(row, column) = Ke(dx * elements[e]->nodes() + nx, dy * elements[e]->nodes() + ny);
					}
				}
				f[subdomain][row] += fe[dx * elements[e]->nodes() + nx];
			}
		}
	}

	for (size_t i = 0; i < _mesh.edges().size(); i++) {
		if (_mesh.edges()[i]->inDomain(subdomain) && _mesh.edges()[i]->regions().size()) {
			processEdge(fe, _mesh, _mesh.edges()[i], _configuration);

			for (size_t nx = 0; nx < _mesh.edges()[i]->nodes(); nx++) {
				for (size_t dx = 0; dx < pointDOFs.size(); dx++) {
					size_t row = nodes[_mesh.edges()[i]->node(nx)]->DOFIndex(subdomain, dx);
					f[subdomain][row] += fe[dx * _mesh.edges()[i]->nodes() + nx];
				}
			}

		}
	}

	// TODO: make it direct
	SparseCSRMatrix<eslocal> csrK = _K;
	K[subdomain] = csrK;
}

}




