
#include "assembler.h"

namespace espreso {

std::vector<Property> Temperature::elementDOFs;
std::vector<Property> Temperature::faceDOFs;
std::vector<Property> Temperature::edgeDOFs;
std::vector<Property> Temperature::pointDOFs = { Property::TEMPERATURE };
std::vector<Property> Temperature::midPointDOFs = { Property::TEMPERATURE };

void Temperature::prepareMeshStructures()
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
}

void Temperature::saveMeshProperties(store::Store &store)
{

}

void Temperature::saveMeshResults(store::Store &store, const std::vector<std::vector<double> > &results)
{

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

static void processElement(DenseMatrix &Ke, std::vector<double> &fe, const espreso::Mesh &mesh, size_t subdomain, const Element* element)
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
		coordinates(i, 0) = mesh.coordinates().get(element->node(i), subdomain).x;
		coordinates(i, 1) = mesh.coordinates().get(element->node(i), subdomain).y;
		coordinates(i, 2) = mesh.coordinates().get(element->node(i), subdomain).z;
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

static void algebraicKernelsAndRegularization(SparseMatrix &K, SparseMatrix &RegMat, SparseMatrix &R, size_t subdomain)
{
	double norm;
	eslocal defect;

	K.get_kernel_from_K(K, RegMat, R, norm, defect, subdomain);
}

void Temperature::assembleStiffnessMatrix(const Element* e, DenseMatrix &Ke, std::vector<double> &fe, std::vector<eslocal> &dofs) const
{
	ESINFO(GLOBAL_ERROR) << "Implement assembleStiffnessMatrix";
}

void Temperature::makeStiffnessMatricesRegular()
{

}

void Temperature::composeSubdomain(size_t subdomain)
{
	eslocal subdomainSize = _mesh.coordinates().localSize(subdomain);
	SparseVVPMatrix<eslocal> _K;
	DenseMatrix Ke;
	std::vector<double> fe;

	_K.resize(subdomainSize, subdomainSize);
	f[subdomain].resize(subdomainSize);

	const std::vector<eslocal> &partition = _mesh.getPartition();

	for (eslocal i = partition[subdomain]; i < partition[subdomain + 1]; i++) {

		const Element* e = _mesh.elements()[i];
		processElement(Ke, fe, _mesh, subdomain, e);

		for (size_t i = 0; i < e->nodes(); i++) {
			eslocal row = e->node(i);
			for (size_t j = 0; j < e->nodes(); j++) {
				eslocal column = e->node(j);
				_K(row, column) = Ke(i, j);
			}
			f[subdomain][row] += fe[i];
		}
	}

	// TODO: make it direct
	SparseCSRMatrix<eslocal> csrK = _K;
	K[subdomain] = csrK;

	switch (_configuration.regularization) {
	case REGULARIZATION::FIX_POINTS:
		analyticsKernels(R1[subdomain], _mesh.coordinates().localSize(subdomain));
		analyticsRegMat(K[subdomain], RegMat[subdomain]);
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


