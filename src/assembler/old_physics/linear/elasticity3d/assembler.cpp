
#include "../../../old_physics/linear/elasticity3d/assembler.h"

#include "../../../../basis/matrices/denseMatrix.h"
#include "../../../../basis/matrices/sparseVVPMatrix.h"
#include "../../../../basis/matrices/sparseCSRMatrix.h"
#include "../../../../configuration/physics/linearelasticity3d.h"
#include "../../../../solver/generic/SparseMatrix.h"
#include "../../../../solver/specific/sparsesolvers.h"

#include "../../../../mesh/elements/element.h"
#include "../../../../mesh/settings/evaluator.h"
#include "../../../../mesh/structures/mesh.h"
#include "../../../../mesh/structures/material.h"
#include "../../../../mesh/structures/elementtypes.h"

#include "../../../../mesh/elements/plane/square4.h"
#include "../../../../mesh/elements/plane/square8.h"
#include "../../../../mesh/elements/plane/triangle3.h"
#include "../../../../mesh/elements/plane/triangle6.h"

#include "../../../../mesh/elements/volume/hexahedron20.h"
#include "../../../../mesh/elements/volume/hexahedron8.h"
#include "../../../../mesh/elements/volume/prisma15.h"
#include "../../../../mesh/elements/volume/prisma6.h"
#include "../../../../mesh/elements/volume/pyramid13.h"
#include "../../../../mesh/elements/volume/pyramid5.h"
#include "../../../../mesh/elements/volume/tetrahedron10.h"
#include "../../../../mesh/elements/volume/tetrahedron4.h"

#include "../../../constraints/equalityconstraints.h"
#include "../../../constraints/inequalityconstraints.h"

#include "../../../../output/resultstore.h"


namespace espreso {

std::vector<Property> LinearElasticity3D::elementDOFs;
std::vector<Property> LinearElasticity3D::faceDOFs;
std::vector<Property> LinearElasticity3D::edgeDOFs;
std::vector<Property> LinearElasticity3D::pointDOFs = { Property::DISPLACEMENT_X, Property::DISPLACEMENT_Y, Property::DISPLACEMENT_Z };
std::vector<Property> LinearElasticity3D::midPointDOFs = { Property::DISPLACEMENT_X, Property::DISPLACEMENT_Y, Property::DISPLACEMENT_Z };

LinearElasticity3D::LinearElasticity3D(Mesh &mesh, Constraints &constraints, const LinearElasticity3DConfiguration &configuration)
: LinearPhysics(
		mesh, constraints, configuration.espreso,
		MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE,
		elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs),
  _configuration(configuration) {};

void LinearElasticity3D::prepareMeshStructures()
{
	Hexahedron8::setDOFs(elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs);
	Hexahedron20::setDOFs(elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs);
	Tetrahedron4::setDOFs(elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs);
	Tetrahedron10::setDOFs(elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs);
	Prisma6::setDOFs(elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs);
	Prisma15::setDOFs(elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs);
	Pyramid5::setDOFs(elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs);
	Pyramid13::setDOFs(elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs);

	Square4::setDOFs(elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs);
	Square8::setDOFs(elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs);
	Triangle3::setDOFs(elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs);
	Triangle6::setDOFs(elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs);

	std::vector<size_t> DOFsOffsets;
	matrixSize = _mesh.assignUniformDOFsIndicesToNodes(matrixSize, pointDOFs, DOFsOffsets);
	_mesh.computeNodesDOFsCounters(pointDOFs);

	if (_solverConfiguration.regularization == REGULARIZATION::FIX_POINTS) {
		_mesh.computeFixPoints(8);
	}

	if (_solverConfiguration.method == ESPRESO_METHOD::HYBRID_FETI) {
		switch (_solverConfiguration.B0_type) {
		case B0_TYPE::CORNERS:
			_mesh.computeVolumeCorners(1, true, true, false);
			break;
		case B0_TYPE::KERNELS:
			_mesh.computeFacesSharedByDomains();
			break;
		case B0_TYPE::COMBINED:
			_mesh.computeFacesSharedByDomains();
			if (!_mesh.corners().size()) {
				_mesh.computeEdgesOnBordersOfFacesSharedByDomains();
				_mesh.computeCornersOnEdges(1, true, true);
			}
			break;
		default:
			ESINFO(GLOBAL_ERROR) << "Not implemented B0";
		}
	}

	_constraints.initMatrices(matrixSize);

	_mesh.removeDuplicateRegions();
}

void LinearElasticity3D::saveMeshProperties(output::ResultStore &store)
{
//	store.storeProperty("displacement", { Property::DISPLACEMENT_X, Property::DISPLACEMENT_Y, Property::DISPLACEMENT_Z }, ElementType::NODES);
//	store.storeProperty("forces", { Property::FORCE_X, Property::FORCE_Y, Property::FORCE_Z }, ElementType::NODES);
//	store.storeProperty("obstacle", { Property::OBSTACLE }, ElementType::NODES);
//	store.storeProperty("normal_direction", { Property::NORMAL_DIRECTION }, ElementType::NODES);
	if (_solverConfiguration.regularization == REGULARIZATION::FIX_POINTS) {
//		store::VTK::fixPoints(store.configuration(), _mesh, "fixPoints");
	}
	if (_solverConfiguration.method == ESPRESO_METHOD::HYBRID_FETI) {
		switch (_solverConfiguration.B0_type) {
		case B0_TYPE::CORNERS:
		case B0_TYPE::COMBINED:
//			store::VTK::mesh(store.configuration(), _mesh, "faces", store::ElementType::FACES);
//			store::VTK::mesh(store.configuration(), _mesh, "edges", store::ElementType::EDGES);
//			store::VTK::corners(store.configuration(), _mesh, "corners");
			break;
		case B0_TYPE::KERNELS:
//			store::VTK::mesh(store.configuration(), _mesh, "faces", store::ElementType::FACES);
			break;
		default:
			ESINFO(GLOBAL_ERROR) << "Not implemented saving properties of B0";
		}
	}
}

void LinearElasticity3D::saveMeshResults(output::ResultStore &store, const std::vector<std::vector<double> > &results)
{
	store.storeValues("displacement", 3, results, ElementType::NODES);
}

void LinearElasticity3D::assembleB1()
{
	EqualityConstraints::insertDirichletToB1(_constraints, _mesh.nodes(), pointDOFs);
	EqualityConstraints::insertElementGluingToB1(_constraints, _mesh.nodes(), pointDOFs, K);
	EqualityConstraints::insertMortarGluingToB1(_constraints, _mesh.faces(), pointDOFs);

	for (size_t i = 0; i < _mesh.evaluators().size(); i++) {
		if (_mesh.evaluators()[i]->property() == Property::OBSTACLE) {
			InequalityConstraints::insertLowerBoundToB1(_constraints, pointDOFs, { Property::OBSTACLE });
		}
	}
}

void LinearElasticity3D::assembleB0()
{
	if (_solverConfiguration.method == ESPRESO_METHOD::HYBRID_FETI) {
		switch (_solverConfiguration.B0_type) {
		case B0_TYPE::CORNERS:
			EqualityConstraints::insertDomainGluingToB0(_constraints, _mesh.corners(), pointDOFs);
			break;
		case B0_TYPE::KERNELS:
			EqualityConstraints::insertKernelsToB0(_constraints, _mesh.faces(), pointDOFs, R1);
			break;
		case B0_TYPE::COMBINED:
			EqualityConstraints::insertKernelsToB0(_constraints, _mesh.faces(), pointDOFs, R1);
			EqualityConstraints::insertDomainGluingToB0(_constraints, _mesh.corners(), pointDOFs);
			break;
		default:
			ESINFO(GLOBAL_ERROR) << "Not implemented construction of B0";
		}
	}
}

static double determinant3x3(DenseMatrix &m)
{
	const double *values = m.values();
	return fabs(
		values[0] * values[4] * values[8] +
		values[1] * values[5] * values[6] +
		values[2] * values[3] * values[7] -
		values[2] * values[4] * values[6] -
		values[1] * values[3] * values[8] -
		values[0] * values[5] * values[7]
   );
}

static void inverse(const DenseMatrix &m, DenseMatrix &inv, double det)
{
	const double *values = m.values();
	inv.resize(m.rows(), m.columns());
	double *invj = inv.values();
	double detJx = 1 / det;
	invj[0] = detJx * (values[8] * values[4] - values[7] * values[5]);
	invj[1] = detJx * (-values[8] * values[1] + values[7] * values[2]);
	invj[2] = detJx * (values[5] * values[1] - values[4] * values[2]);
	invj[3] = detJx * (-values[8] * values[3] + values[6] * values[5]);
	invj[4] = detJx * (values[8] * values[0] - values[6] * values[2]);
	invj[5] = detJx * (-values[5] * values[0] + values[3] * values[2]);
	invj[6] = detJx * (values[7] * values[3] - values[6] * values[4]);
	invj[7] = detJx * (-values[7] * values[0] + values[6] * values[1]);
	invj[8] = detJx * (values[4] * values[0] - values[3] * values[1]);
}

// B =
// dX   0   0
//  0  dY   0
//  0   0  dZ
// dY  dX   0
//  0  dZ  dY
// dZ   0  dX
static void distribute(DenseMatrix &B, DenseMatrix &dND)
{
	eslocal columns = dND.rows() * dND.columns();
	const double *dNDx = dND.values();
	const double *dNDy = dND.values() + dND.columns();
	const double *dNDz = dND.values() + 2 * dND.columns();

	double *v = B.values();

	memcpy(&v[0], dNDx,                               sizeof(double) * dND.columns());
	memcpy(&v[3 * columns + dND.columns()],     dNDx, sizeof(double) * dND.columns());
	memcpy(&v[5 * columns + 2 * dND.columns()], dNDx, sizeof(double) * dND.columns());

	memcpy(&v[1 * columns + dND.columns()],     dNDy, sizeof(double) * dND.columns());
	memcpy(&v[3 * columns],                     dNDy, sizeof(double) * dND.columns());
	memcpy(&v[4 * columns + 2 * dND.columns()], dNDy, sizeof(double) * dND.columns());

	memcpy(&v[2 * columns + 2 * dND.columns()], dNDz, sizeof(double) * dND.columns());
	memcpy(&v[4 * columns + dND.columns()],     dNDz, sizeof(double) * dND.columns());
	memcpy(&v[5 * columns],                     dNDz, sizeof(double) * dND.columns());
}

static void processElement(DenseMatrix &Ke, std::vector<double> &fe, const espreso::Mesh &mesh, const Element* element)
{
	DenseMatrix Ce(6, 6), coordinates, J, invJ, dND, B;
	std::vector<double> inertia(3, 0);
	double detJ;

	const Material* material = mesh.materials()[element->param(Element::MATERIAL)];
	const std::vector<DenseMatrix> &dN = element->dN();
	const std::vector<DenseMatrix> &N = element->N();
	const std::vector<double> &weighFactor = element->weighFactor();

	// TODO: set the omega from example
	Point omega(50, 50, 0);

	double ex = material->get(MATERIAL_PARAMETER::YOUNG_MODULUS_X)->evaluate(0);
	double mi = material->get(MATERIAL_PARAMETER::POISSON_RATIO_XY)->evaluate(0);

	double E = ex / ((1 + mi) * (1 - 2 * mi));

	Ce(0, 1) = Ce(0, 2) = Ce(1, 0) = Ce(1, 2) = Ce(2, 0) = Ce(2, 1) = E * mi;
	Ce(0, 0) = Ce(1, 1) = Ce(2, 2) = E * (1.0 - mi);
	Ce(3, 3) = Ce(4, 4) = Ce(5, 5) = E * (0.5 - mi);

	inertia[0] = inertia[1] = 0; // inertia[2] = 0;
	inertia[2] = 9.8066 * material->get(MATERIAL_PARAMETER::DENSITY)->evaluate(0);

	coordinates.resize(element->nodes(), 3);

	Point mid;
	for (size_t i = 0; i < element->nodes(); i++) {
		coordinates(i, 0) = mesh.coordinates()[element->node(i)].x;
		coordinates(i, 1) = mesh.coordinates()[element->node(i)].y;
		coordinates(i, 2) = mesh.coordinates()[element->node(i)].z;
		mid += mesh.coordinates()[element->node(i)];
	}
	mid /= element->nodes();

	eslocal Ksize = 3 * element->nodes();

	Ke.resize(Ksize, Ksize);
	Ke = 0;
	fe.resize(Ksize);
	fill(fe.begin(), fe.end(), 0);

	for (size_t gp = 0; gp < element->gaussePoints(); gp++) {
		J.multiply(dN[gp], coordinates);
		detJ = determinant3x3(J);
		inverse(J, invJ, detJ);

		dND.multiply(invJ, dN[gp]);
		B.resize(Ce.rows(), Ksize);
		distribute(B, dND);
		Ke.multiply(B, Ce * B, detJ * weighFactor[gp], 1, true);

		for (eslocal i = 0; i < Ksize; i++) {
			// TODO: set rotation from example
			fe[i] += detJ * weighFactor[gp] * N[gp](0, i % element->nodes()) * inertia[i / element->nodes()];
			//fe[i] += detJ * weighFactor[gp] * N[gp](0, i % e->size()) * 7850 * rotation[i / e->size()];
		}
	}
}

static void analyticsKernels(SparseMatrix &R1, const Coordinates &coordinates, size_t subdomain)
{
	size_t nodes = coordinates.localSize(subdomain);
	R1.rows = 3 * nodes;
	R1.cols = 6;
	R1.nnz = R1.rows * R1.cols;
	R1.type = 'G';

	R1.dense_values.reserve(R1.nnz);

	for (size_t c = 0; c < 3; c++) {
		std::vector<double> kernel = { 0, 0, 0 };
		kernel[c] = 1;
		for (size_t i = 0; i < nodes; i++) {
			R1.dense_values.insert(R1.dense_values.end(), kernel.begin(), kernel.end());
		}
	}

	for (size_t i = 0; i < coordinates.localSize(subdomain); i++) {
		const Point &p = coordinates.get(i, subdomain);
		R1.dense_values.push_back(-p.y);
		R1.dense_values.push_back( p.x);
		R1.dense_values.push_back(   0);
	}

	for (size_t i = 0; i < coordinates.localSize(subdomain); i++) {
		const Point &p = coordinates.get(i, subdomain);
		R1.dense_values.push_back(-p.z);
		R1.dense_values.push_back(   0);
		R1.dense_values.push_back( p.x);
	}

	for (size_t i = 0; i < coordinates.localSize(subdomain); i++) {
		const Point &p = coordinates.get(i, subdomain);
		R1.dense_values.push_back(   0);
		R1.dense_values.push_back(-p.z);
		R1.dense_values.push_back( p.y);
	}
}

static void analyticsRegMat(SparseMatrix &K, SparseMatrix &RegMat, const std::vector<Element*> &fixPoints, const Coordinates &coordinates, size_t subdomain)
{
	ESTEST(MANDATORY) << "Too few FIX POINTS: " << fixPoints.size() << (fixPoints.size() > 3 ? TEST_PASSED : TEST_FAILED);

	SparseMatrix Nt; // CSR matice s DOFY
	Nt.rows = 6;
	Nt.cols = K.cols;
	Nt.nnz  = 9 * fixPoints.size();
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
	ROWS.push_back(ROWS.back() + fixPoints.size());
	ROWS.push_back(ROWS.back() + 2 * fixPoints.size());
	ROWS.push_back(ROWS.back() + 2 * fixPoints.size());
	ROWS.push_back(ROWS.back() + 2 * fixPoints.size());

	for (size_t c = 0; c < 3; c++) {
		for (size_t i = 0; i < fixPoints.size(); i++) {
			COLS.push_back(fixPoints[i]->DOFIndex(subdomain, c) + 1);
		}
	}
	VALS.insert(VALS.end(), 3 * fixPoints.size(), 1);

	for (size_t i = 0; i < fixPoints.size(); i++) {
		const Point &p = coordinates[fixPoints[i]->node(0)];
		COLS.push_back(fixPoints[i]->DOFIndex(subdomain, 0) + 1);
		COLS.push_back(fixPoints[i]->DOFIndex(subdomain, 1) + 1);
		VALS.push_back(-p.y);
		VALS.push_back( p.x);
	}

	for (size_t i = 0; i < fixPoints.size(); i++) {
		const Point &p = coordinates[fixPoints[i]->node(0)];
		COLS.push_back(fixPoints[i]->DOFIndex(subdomain, 0) + 1);
		COLS.push_back(fixPoints[i]->DOFIndex(subdomain, 2) + 1);
		VALS.push_back(-p.z);
		VALS.push_back( p.x);
	}

	for (size_t i = 0; i < fixPoints.size(); i++) {
		const Point &p = coordinates[fixPoints[i]->node(0)];
		COLS.push_back(fixPoints[i]->DOFIndex(subdomain, 1) + 1);
		COLS.push_back(fixPoints[i]->DOFIndex(subdomain, 2) + 1);
		VALS.push_back(-p.z);
		VALS.push_back( p.y);
	}

	SparseMatrix N;
	Nt.MatTranspose( N );
	RegMat.MatMat(Nt, 'N', N);
	RegMat.MatTranspose();
	RegMat.RemoveLower();
	RegMat.mtype = MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE;

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

void LinearElasticity3D::assembleStiffnessMatrix(const Element* e, DenseMatrix &Ke, std::vector<double> &fe, std::vector<eslocal> &dofs) const
{
	processElement(Ke, fe, _mesh, e);
	dofs.resize(e->nodes() * pointDOFs.size());
	for (size_t dof = 0, i = 0; dof < pointDOFs.size(); dof++) {
		for (size_t n = 0; n < e->nodes(); n++, i++) {
			dofs[i] = e->node(n) * pointDOFs.size() + dof;
		}
	}


	std::vector<Property> forces = { Property::FORCE_X, Property::FORCE_Y, Property::FORCE_Z };
	for (size_t n = 0; n < e->nodes(); n++) {
		for (size_t dof = 0; dof < pointDOFs.size(); dof++) {
			fe[n * pointDOFs.size() + dof] = e->sumProperty(forces[dof],n, 0, 0);
		}
	}
}

void LinearElasticity3D::makeStiffnessMatricesRegular()
{
	#pragma omp parallel for
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

void LinearElasticity3D::composeSubdomain(size_t subdomain)
{
	SparseVVPMatrix<eslocal> _K;
	DenseMatrix Ke;
	std::vector<double> fe;

	_K.resize(matrixSize[subdomain], matrixSize[subdomain]);
	f[subdomain].resize(matrixSize[subdomain]);

	const std::vector<eslocal> &partition = _mesh.getPartition();
	const std::vector<Element*> &elements = _mesh.elements();
	const std::vector<Element*> &nodes = _mesh.nodes();

	for (eslocal e = partition[subdomain]; e < partition[subdomain + 1]; e++) {

		processElement(Ke, fe, _mesh, elements[e]);

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

	std::vector<Property> forces = { Property::FORCE_X, Property::FORCE_Y, Property::FORCE_Z };
	for (size_t n = 0; n < _mesh.coordinates().localSize(subdomain); n++) {
		Element *node = _mesh.nodes()[_mesh.coordinates().clusterIndex(n, subdomain)];
		for (size_t dof = 0; dof < pointDOFs.size(); dof++) {
			f[subdomain][node->DOFIndex(subdomain, dof)] += node->sumProperty(forces[dof], 0, 0, 0) / node->numberOfGlobalDomainsWithDOF(dof);
		}
	}

	// TODO: make it direct
	SparseCSRMatrix<eslocal> csrK = _K;
	K[subdomain] = csrK;
}

}



