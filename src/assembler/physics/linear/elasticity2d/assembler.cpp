
#include "assembler.h"

using namespace espreso;

std::vector<Property> LinearElasticity2D::elementDOFs;
std::vector<Property> LinearElasticity2D::faceDOFs;
std::vector<Property> LinearElasticity2D::edgeDOFs;
std::vector<Property> LinearElasticity2D::pointDOFs = { Property::DISPLACEMENT_X, Property::DISPLACEMENT_Y };
std::vector<Property> LinearElasticity2D::midPointDOFs = { Property::DISPLACEMENT_X, Property::DISPLACEMENT_Y };

LinearElasticity2D::ELEMENT_BEHAVIOUR LinearElasticity2D::elementBehaviour = ELEMENT_BEHAVIOUR::PLANE_STRESS;
Point LinearElasticity2D::angularVelocity(0, 0, 0);

void LinearElasticity2D::prepareMeshStructures()
{
	Square4::setDOFs(elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs);
	Square8::setDOFs(elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs);
	Triangle3::setDOFs(elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs);
	Triangle6::setDOFs(elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs);

	matrixSize = _mesh.assignUniformDOFsIndicesToNodes(matrixSize, pointDOFs);
	_mesh.computeNodesDOFsCounters(pointDOFs);

	if (config::solver::REGULARIZATION == config::solver::REGULARIZATIONalternative::FIX_POINTS) {
		_mesh.computeFixPoints(config::mesh::FIX_POINTS / 2);
	}

	if (config::solver::FETI_METHOD == config::solver::FETI_METHODalternative::HYBRID_FETI) {
		switch (config::solver::B0_TYPE) {
		case config::solver::B0_TYPEalternative::CORNERS:
			_mesh.computePlaneCorners(config::mesh::CORNERS, config::mesh::VERTEX_CORNERS, config::mesh::EDGE_CORNERS);
			break;
		case config::solver::B0_TYPEalternative::KERNELS:
			_mesh.computeEdgesSharedByDomains();
			break;
		}
	}
}

void LinearElasticity2D::saveMeshProperties(output::Store &store)
{

}

void LinearElasticity2D::saveMeshResults(output::Store &store, const std::vector<std::vector<double> > &results)
{
	store.storeValues("displacement", 2, results, output::Store::ElementType::NODES);
}

void LinearElasticity2D::assembleGluingMatrices()
{
	_constraints.initMatrices(matrixSize);

	_constraints.insertDirichletToB1(_mesh.nodes(), pointDOFs);
	_constraints.insertElementGluingToB1(_mesh.nodes(), pointDOFs);

	if (config::solver::FETI_METHOD == config::solver::FETI_METHODalternative::HYBRID_FETI) {
		switch (config::solver::B0_TYPE) {
		case config::solver::B0_TYPEalternative::CORNERS:
			_constraints.insertDomainGluingToB0(_mesh.corners(), pointDOFs);
			break;
		case config::solver::B0_TYPEalternative::KERNELS:
			_constraints.insertKernelsToB0(_mesh.edges(), pointDOFs, R1);
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

static void fillC(DenseMatrix &C, LinearElasticity2D::ELEMENT_BEHAVIOUR behaviour, Material::MODEL model, DenseMatrix &dens, DenseMatrix &E, DenseMatrix &mi, DenseMatrix &T)
{
	switch (model) {

	case Material::MODEL::LINEAR_ELASTIC_ISOTROPIC:
	{

		switch (behaviour) {

		case LinearElasticity2D::ELEMENT_BEHAVIOUR::PLANE_STRAIN:
		{
			double k = E(0, 0) * (1 - mi(0, 0)) / ((1 + mi(0, 0)) * (1 - 2 * mi(0, 0)));
			C.resize(3, 3);
			C(0, 0) = C(1, 1) = k * 1;
			C(1, 0) = C(0, 1) = k * (mi(0, 0) / (1 - mi(0, 0)));
			C(2, 2) = k * ((1 - 2 * mi(0, 0)) / (2 * (1 - mi(0, 0))));
			return;
		}

		case LinearElasticity2D::ELEMENT_BEHAVIOUR::PLANE_STRESS:
		case LinearElasticity2D::ELEMENT_BEHAVIOUR::PLANE_STRESS_WITH_THICKNESS:
		{
			double k = E(0, 0) / (1 - mi(0, 0) * mi(0, 0));
			C.resize(3, 3);
			C(0, 0) = C(1, 1) = k * 1;
			C(1, 0) = C(0, 1) = k * mi(0, 0);
			C(2, 2) = k * ((1 -  mi(0, 0)) / 2);
			return;
		}

		case LinearElasticity2D::ELEMENT_BEHAVIOUR::AXISYMMETRIC:
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

	case Material::MODEL::LINEAR_ELASTIC_ORTHOTROPIC:
	{
		return;
	}

	case Material::MODEL::LINEAR_ELASTIC_ANISOTROPIC:
	{
		return;
	}
	}
}

static void processElement(DenseMatrix &Ke, std::vector<double> &fe, const espreso::Mesh &mesh, size_t subdomain, const Element* element)
{
	DenseMatrix Ce(4, 4), XY(1, 2), coordinates, J, invJ, dND, B, epsilon, rhsT;
	DenseMatrix
			matDENS(element->nodes(), 1), matE(element->nodes(), 2), matMI(element->nodes(), 1),
			matTE(element->nodes(), 2), matT(element->nodes(), 1), matThickness(element->nodes(), 1),
			matInitT(element->nodes(), 1);
	DenseMatrix gpDENS(1, 1), gpE(1, 2), gpMI(1, 1), gpTE(1, 2), gpT(1, 1), gpInitT(1, 1), gpThickness(1, 1);
	double detJ, inertia[2];

	const Material &material = mesh.materials()[element->param(Element::MATERIAL)];
	const std::vector<DenseMatrix> &dN = element->dN();
	const std::vector<DenseMatrix> &N = element->N();
	const std::vector<double> &weighFactor = element->weighFactor();

	coordinates.resize(element->nodes(), 2);

	for (size_t i = 0; i < element->nodes(); i++) {
		matDENS(i, 0) = material.density(element->node(i));
		matE(i, 0) = material.youngModulusX(element->node(i));
		matE(i, 1) = material.youngModulusY(element->node(i));
		matMI(i, 0) = material.poissonRatioXY(element->node(i));
		matTE(i, 0) = material.termalExpansionX(element->node(i));
		matTE(i, 1) = material.termalExpansionY(element->node(i));
		matInitT(i, 0) = element->settings(Property::INITIAL_TEMPERATURE).back()->evaluate(element->node(i));
		if (mesh.nodes()[element->node(i)]->settings().isSet(Property::TEMPERATURE)) {
			matT(i, 0) =  mesh.nodes()[element->node(i)]->settings(Property::TEMPERATURE).back()->evaluate(element->node(i));
		} else {
			matT(i, 0) = matInitT(i, 0);
		}

		inertia[0] = element->settings(Property::ACCELERATION_X).back()->evaluate(element->node(i));
		inertia[1] = element->settings(Property::ACCELERATION_Y).back()->evaluate(element->node(i));

		matThickness(i, 0) = mesh.nodes()[element->node(i)]->settings(Property::THICKNESS).back()->evaluate(element->node(i));

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

	for (eslocal gp = 0; gp < element->gaussePoints(); gp++) {
		J.multiply(dN[gp], coordinates);
		detJ = determinant2x2(J);
		inverse(J, invJ, detJ);
		dND.multiply(invJ, dN[gp]);
		XY.multiply(N[gp], coordinates);
		gpDENS.multiply(N[gp], matDENS, 1, 0);
		gpE.multiply(N[gp], matE, 1, 0);
		gpMI.multiply(N[gp], matMI, 1, 0);
		gpTE.multiply(N[gp], matTE, 1, 0);
		gpT.multiply(N[gp], matT, 1, 0);
		gpInitT.multiply(N[gp], matInitT, 1, 0);

		gpThickness.multiply(N[gp], matThickness, 1, 0);

		fillC(Ce, LinearElasticity2D::elementBehaviour, material.model(), gpDENS, gpE, gpMI, gpT);
		B.resize(Ce.rows(), Ksize);
		epsilon.resize(Ce.rows(), 1);
		switch (LinearElasticity2D::elementBehaviour) {

		case LinearElasticity2D::ELEMENT_BEHAVIOUR::PLANE_STRESS:
		case LinearElasticity2D::ELEMENT_BEHAVIOUR::PLANE_STRAIN:
			gpThickness(0, 0) = 1;
		case LinearElasticity2D::ELEMENT_BEHAVIOUR::PLANE_STRESS_WITH_THICKNESS:
			epsilon(0, 0) = (gpT(0, 0) - gpInitT(0, 0)) * gpTE(0, 0);
			epsilon(1, 0) = (gpT(0, 0) - gpInitT(0, 0)) * gpTE(0, 1);
			epsilon(2, 0) = 0;
			distribute(B, dND);
			Ke.multiply(B, Ce * B, detJ * weighFactor[gp] * gpThickness(0, 0), 1, true);
			rhsT.multiply(B, Ce * epsilon, detJ * weighFactor[gp] * gpThickness(0, 0), 0, true, false);
			for (eslocal i = 0; i < Ksize; i++) {
				fe[i] += gpDENS(0, 0) * detJ * weighFactor[gp] * gpThickness(0, 0) * N[gp](0, i % element->nodes()) * inertia[i / element->nodes()];
				fe[i] += gpDENS(0, 0) * detJ * weighFactor[gp] * gpThickness(0, 0) * N[gp](0, i % element->nodes()) * XY(0, i / element->nodes()) * pow(LinearElasticity2D::angularVelocity.z, 2);
				fe[i] += rhsT(i, 0);
			}
			break;

		case LinearElasticity2D::ELEMENT_BEHAVIOUR::AXISYMMETRIC:
			epsilon(0, 0) = (gpT(0, 0) - gpInitT(0, 0)) * gpTE(0, 0);
			epsilon(1, 0) = (gpT(0, 0) - gpInitT(0, 0)) * gpTE(0, 1);
			epsilon(2, 0) = epsilon(3, 0) = 0;
			distribute(B, dND, N[gp], XY);
			Ke.multiply(B, Ce * B, detJ * weighFactor[gp] * 2 * M_PI * XY(0, 0), 1, true);
			rhsT.multiply(B, Ce * epsilon, detJ * weighFactor[gp] * 2 * M_PI * XY(0, 0), 0, true, false);
			for (eslocal i = 0; i < Ksize; i++) {
				fe[i] += gpDENS(0, 0) * detJ * weighFactor[gp] * 2 * M_PI * XY(0, 0) * N[gp](0, i % element->nodes()) * inertia[i / element->nodes()];
				fe[i] += rhsT(i, 0);
			}
			for (eslocal i = 0; i < Ksize / 2; i++) {
				fe[i] += gpDENS(0, 0) * detJ * weighFactor[gp] * 2 * M_PI * XY(0, 0) * N[gp](0, i % element->nodes()) * XY(0, 0) * pow(LinearElasticity2D::angularVelocity.y, 2);
				fe[Ksize / 2 + i] += gpDENS(0, 0) * detJ * weighFactor[gp] * 2 * M_PI * XY(0, 0) * N[gp](0, i % element->nodes()) * XY(0, 1) * pow(LinearElasticity2D::angularVelocity.x, 2);
			}
			break;
		}
	}
}

static void processEdge(std::vector<double> &fe, const espreso::Mesh &mesh, size_t subdomain, const Element* edge)
{
	DenseMatrix coordinates(edge->nodes(), 2), dND(1, 2), P(edge->nodes(), 1), normal(2, 1), matThickness(edge->nodes(), 1), XY(1, 2);
	DenseMatrix gpP(1, 1), gpQ(1, 2), gpThickness(1, 1);

	const std::vector<DenseMatrix> &dN = edge->dN();
	const std::vector<DenseMatrix> &N = edge->N();
	const std::vector<double> &weighFactor = edge->weighFactor();

	for (size_t n = 0; n < edge->nodes(); n++) {
		coordinates(n, 0) = mesh.coordinates()[edge->node(n)].x;
		coordinates(n, 1) = mesh.coordinates()[edge->node(n)].y;
		P(n, 0) = edge->settings(Property::PRESSURE).back()->evaluate(edge->node(n));
		matThickness(n, 0) = mesh.nodes()[edge->node(n)]->settings(Property::THICKNESS).back()->evaluate(edge->node(n));
	}

	eslocal Ksize = 2 * edge->nodes();
	fe.resize(Ksize);
	std::fill(fe.begin(), fe.end(), 0);

	for (eslocal gp = 0; gp < edge->gaussePoints(); gp++) {
		dND.multiply(dN[gp], coordinates);
		double J = dND.norm();
		normal(0, 0) = -dND(0, 1) / J;
		normal(1, 0) =  dND(0, 0) / J;
		gpP.multiply(N[gp], P);
		gpQ.multiply(normal, gpP);
		gpThickness.multiply(N[gp], matThickness);

		switch (LinearElasticity2D::elementBehaviour) {

		case LinearElasticity2D::ELEMENT_BEHAVIOUR::PLANE_STRESS:
		case LinearElasticity2D::ELEMENT_BEHAVIOUR::PLANE_STRAIN:
			gpThickness(0, 0) = 1;
		case LinearElasticity2D::ELEMENT_BEHAVIOUR::PLANE_STRESS_WITH_THICKNESS:
			for (eslocal i = 0; i < Ksize; i++) {
				fe[i] += gpThickness(0, 0) * J * weighFactor[gp] * N[gp](0, i % edge->nodes()) * gpQ(0, i / edge->nodes());
			}
			break;

		case LinearElasticity2D::ELEMENT_BEHAVIOUR::AXISYMMETRIC:
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
			COLS.push_back(fixPoints[i]->DOFIndex(subdomain, 0) + IJVMatrixIndexing);
			COLS.push_back(fixPoints[i]->DOFIndex(subdomain, 1) + IJVMatrixIndexing);
			VALS.insert(VALS.end(), kernel.begin(), kernel.end());
		}
	}

	for (size_t i = 0; i < fixPoints.size(); i++) {
		const Point &p = coordinates[fixPoints[i]->node(0)];
		COLS.push_back(fixPoints[i]->DOFIndex(subdomain, 0) + IJVMatrixIndexing);
		COLS.push_back(fixPoints[i]->DOFIndex(subdomain, 1) + IJVMatrixIndexing);
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

void LinearElasticity2D::assembleStiffnessMatrix(const Element* e, DenseMatrix &Ke, std::vector<double> &fe, std::vector<eslocal> &dofs) const
{
	ESINFO(GLOBAL_ERROR) << "Implement assembleStiffnessMatrix";
}

void LinearElasticity2D::makeStiffnessMatricesRegular()
{
	for (size_t subdomain = 0; subdomain < K.size(); subdomain++) {
		switch (config::solver::REGULARIZATION) {
		case config::solver::REGULARIZATIONalternative::FIX_POINTS:
			if (config::assembler::DOFS_ORDER == config::assembler::DOFS_ORDERalternative::GROUP_DOFS) {
				ESINFO(GLOBAL_ERROR) << "Implement regularization for GROUP_DOFS alternative";
			}

			analyticsKernels(R1[subdomain], _mesh.coordinates(), subdomain);
			analyticsRegMat(K[subdomain], RegMat[subdomain], _mesh.fixPoints(subdomain), _mesh.coordinates(), subdomain);
			K[subdomain].RemoveLower();
			RegMat[subdomain].RemoveLower();
			K[subdomain].MatAddInPlace(RegMat[subdomain], 'N', 1);
			RegMat[subdomain].ConvertToCOO(1);
			break;
		case config::solver::REGULARIZATIONalternative::NULL_PIVOTS:
			K[subdomain].RemoveLower();
			algebraicKernelsAndRegularization(K[subdomain], RegMat[subdomain], R1[subdomain], subdomain);
			break;
		}
	}
}

void LinearElasticity2D::composeSubdomain(size_t subdomain)
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

		processElement(Ke, fe, _mesh, subdomain, elements[e]);

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
		if (_mesh.edges()[i]->inDomain(subdomain) && _mesh.edges()[i]->settings().size()) {
			processEdge(fe, _mesh, subdomain, _mesh.edges()[i]);

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




