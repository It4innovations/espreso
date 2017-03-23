
#include "assembler.h"

namespace espreso {

std::vector<Property> AdvectionDiffusion3D::elementDOFs;
std::vector<Property> AdvectionDiffusion3D::faceDOFs;
std::vector<Property> AdvectionDiffusion3D::edgeDOFs;
std::vector<Property> AdvectionDiffusion3D::pointDOFs = { Property::TEMPERATURE };
std::vector<Property> AdvectionDiffusion3D::midPointDOFs = { Property::TEMPERATURE };

void AdvectionDiffusion3D::prepareMeshStructures()
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

	matrixSize = _mesh.assignUniformDOFsIndicesToNodes(matrixSize, pointDOFs);
	_mesh.computeNodesDOFsCounters(pointDOFs);

	if (config::solver::FETI_METHOD == config::solver::FETI_METHODalternative::HYBRID_FETI) {
		switch (config::solver::B0_TYPE) {
		case config::solver::B0_TYPEalternative::CORNERS:
			_mesh.computeVolumeCorners(config::mesh::CORNERS, config::mesh::VERTEX_CORNERS, config::mesh::EDGE_CORNERS, config::mesh::FACE_CORNERS);
			break;
		case config::solver::B0_TYPEalternative::KERNELS:
			_mesh.computeFacesSharedByDomains();
			break;
		default:
			break;
		}
	}

	_constraints.initMatrices(matrixSize);
}

void AdvectionDiffusion3D::assembleB1()
{
	EqualityConstraints::insertDirichletToB1(_constraints, _mesh.nodes(), pointDOFs);
	EqualityConstraints::insertElementGluingToB1(_constraints, _mesh.nodes(), pointDOFs, K);
}

void AdvectionDiffusion3D::assembleB0()
{
	if (config::solver::FETI_METHOD == config::solver::FETI_METHODalternative::HYBRID_FETI) {
		switch (config::solver::B0_TYPE) {
		case config::solver::B0_TYPEalternative::CORNERS:
			EqualityConstraints::insertDomainGluingToB0(_constraints, _mesh.corners(), pointDOFs);
			break;
		case config::solver::B0_TYPEalternative::KERNELS:
			std::for_each(R1.begin(), R1.end(), [] (SparseMatrix &m) { m.ConvertCSRToDense(0); });
			EqualityConstraints::insertKernelsToB0(_constraints, _mesh.faces(), pointDOFs, R1);
			break;
		default:
			break;
		}
	}
}

void AdvectionDiffusion3D::saveMeshProperties(output::Store &store)
{
	//store.storeProperty("translationMotion", { Property::TRANSLATION_MOTION_X, Property::TRANSLATION_MOTION_Y, Property::TRANSLATION_MOTION_Z }, output::Store::ElementType::ELEMENTS);
	store.storeProperty("headSource", { Property::HEAT_SOURCE }, output::Store::ElementType::ELEMENTS);
	store.storeProperty("temperature", { Property::TEMPERATURE }, output::Store::ElementType::NODES);
}

void AdvectionDiffusion3D::saveMeshResults(output::Store &store, const std::vector<std::vector<double> > &results)
{
	store.storeValues("temperature", 1, results, output::Store::ElementType::NODES);
	store.finalize();
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

static void processElement(DenseMatrix &Ke, std::vector<double> &fe, const espreso::Mesh &mesh, const Element* element)
{
	DenseMatrix Ce(3, 3), coordinates, J, invJ, dND;
	double detJ, inertia;

	const std::vector<DenseMatrix> &dN = element->dN();
	const std::vector<DenseMatrix> &N = element->N();
	const std::vector<double> &weighFactor = element->weighFactor();

	Ce(0, 0) = Ce(1, 1) = Ce(2, 2) = 1;
	inertia = 0;

	coordinates.resize(element->nodes(), 3);
	for (size_t i = 0; i < element->nodes(); i++) {
		const Point &p = mesh.coordinates()[element->node(i)];
		coordinates(i, 0) = p.x;
		coordinates(i, 1) = p.y;
		coordinates(i, 2) = p.z;
	}

	eslocal Ksize = element->nodes();

	Ke.resize(Ksize, Ksize);
	Ke = 0;
	fe.resize(Ksize);
	fill(fe.begin(), fe.end(), 0);

	for (size_t gp = 0; gp < element->gaussePoints(); gp++) {
		J.multiply(dN[gp], coordinates);
		detJ = determinant3x3(J);
		inverse(J, invJ, detJ);

		dND.multiply(invJ, dN[gp]);
		Ke.multiply(dND, Ce * dND, detJ * weighFactor[gp], 1, true);

		for (eslocal i = 0; i < Ksize; i++) {
			fe[i] += detJ * weighFactor[gp] * N[gp](0, i) * inertia;
		}
	}
}

static void algebraicKernelsAndRegularization(SparseMatrix &K, SparseMatrix &R1, SparseMatrix &R2, SparseMatrix &RegMat, size_t subdomain)
{
	double norm;
	eslocal defect;

	K.get_kernel_from_K(K, RegMat, R1, norm, defect, subdomain);
}

void AdvectionDiffusion3D::assembleStiffnessMatrix(const Element* e, DenseMatrix &Ke, std::vector<double> &fe, std::vector<eslocal> &dofs) const
{
	processElement(Ke, fe, _mesh, e);
	dofs.resize(e->nodes());
	for (size_t n = 0; n < e->nodes(); n++) {
		dofs[n] = e->node(n);
	}
}

static void analyticsKernels(SparseMatrix &R1, size_t nodes)
{
	R1.rows = nodes;
	R1.cols = 1;
	R1.nnz = R1.rows * R1.cols;
	R1.type = 'G';

	R1.dense_values.resize(R1.nnz, 1 / sqrt(nodes));
}

static void analyticsRegMat(SparseMatrix &K, SparseMatrix &RegMat)
{
	RegMat.rows = K.rows;
	RegMat.cols = K.cols;
	RegMat.nnz  = 1;
	RegMat.type = K.type;

	RegMat.I_row_indices.push_back(1);
	RegMat.J_col_indices.push_back(1);
	RegMat.V_values.push_back(K.getDiagonalMaximum());
	RegMat.ConvertToCSR(1);
}

void AdvectionDiffusion3D::makeStiffnessMatricesRegular()
{
	#pragma omp parallel for
	for (size_t subdomain = 0; subdomain < K.size(); subdomain++) {
		switch (config::solver::REGULARIZATION) {
		case config::solver::REGULARIZATIONalternative::FIX_POINTS:
			analyticsKernels(R1[subdomain], _mesh.coordinates().localSize(subdomain));
			analyticsRegMat(K[subdomain], RegMat[subdomain]);
			K[subdomain].RemoveLower();
			RegMat[subdomain].RemoveLower();
			K[subdomain].MatAddInPlace(RegMat[subdomain], 'N', 1);
			RegMat[subdomain].ConvertToCOO(1);
			break;
		case config::solver::REGULARIZATIONalternative::NULL_PIVOTS:
            std::cout << "HOOVNOO" <<std::endl;
			K[subdomain].RemoveLower();
			algebraicKernelsAndRegularization(K[subdomain], R1[subdomain], R2[subdomain], RegMat[subdomain], subdomain);
			break;
		default:
			break;
		}
	}
}

void AdvectionDiffusion3D::composeSubdomain(size_t subdomain)
{
	SparseVVPMatrix<eslocal> _K;
	DenseMatrix Ke;
	std::vector<double> fe;

	_K.resize(matrixSize[subdomain], matrixSize[subdomain]);
	f[subdomain].resize(matrixSize[subdomain]);

	const std::vector<eslocal> &partition = _mesh.getPartition();

	for (eslocal i = partition[subdomain]; i < partition[subdomain + 1]; i++) {

		const Element* e = _mesh.elements()[i];
		processElement(Ke, fe, _mesh, e);

		for (size_t i = 0; i < e->nodes(); i++) {
			eslocal row = _mesh.nodes()[e->node(i)]->DOFIndex(subdomain, 0);
			for (size_t j = 0; j < e->nodes(); j++) {
				eslocal column = _mesh.nodes()[e->node(j)]->DOFIndex(subdomain, 0);
				_K(row, column) = Ke(i, j);
			}
			f[subdomain][row] += fe[i];
		}
	}

	// TODO: make it direct
	SparseCSRMatrix<eslocal> csrK = _K;
	K[subdomain] = csrK;
}

}

