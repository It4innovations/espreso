
#include "assembler.h"

using namespace espreso;

std::vector<Property> UniformSymmetric3DOFs::elementDOFs;
std::vector<Property> UniformSymmetric3DOFs::faceDOFs;
std::vector<Property> UniformSymmetric3DOFs::edgeDOFs;
std::vector<Property> UniformSymmetric3DOFs::pointDOFs = { Property::DISPLACEMENT_X, Property::DISPLACEMENT_Y, Property::DISPLACEMENT_Z };
std::vector<Property> UniformSymmetric3DOFs::midPointDOFs = { Property::DISPLACEMENT_X, Property::DISPLACEMENT_Y, Property::DISPLACEMENT_Z };

void UniformSymmetric3DOFs::prepareMeshStructures()
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

void UniformSymmetric3DOFs::assembleStiffnessMatrix(const Element* e, DenseMatrix &Ke, std::vector<double> &fe)
{
	ESINFO(GLOBAL_ERROR) << "Implement assembleStiffnessMatrix";
}

static void algebraicKernelsAndRegularization(SparseMatrix &K, SparseMatrix &RegMat, SparseMatrix &R, size_t subdomain)
{
	double norm;
	eslocal defect;

	K.get_kernel_from_K(K, RegMat, R, norm, defect, subdomain);
}

void UniformSymmetric3DOFs::makeStiffnessMatricesRegular()
{
	for (size_t subdomain = 0; subdomain < K.size(); subdomain++) {
		algebraicKernelsAndRegularization(K[subdomain], RegMat[subdomain], R1[subdomain], subdomain);
	}
}

void UniformSymmetric3DOFs::composeSubdomain(size_t subdomain)
{
	SparseVVPMatrix<eslocal> _K;
	eslocal nK = _mesh.coordinates().localSize(subdomain) * pointDOFs.size();
	_K.resize(nK, nK);

	const std::vector<eslocal> &parts = _mesh.getPartition();
	const std::vector<Element*> &elements = _mesh.elements();

	size_t row, column;

	for (size_t e = parts[subdomain]; e < parts[subdomain + 1]; e++) {
		for (size_t i = 0; i < elements[e]->nodes() * pointDOFs.size(); i++) {
			row = pointDOFs.size() * elements[e]->node(i / pointDOFs.size()) + (i % pointDOFs.size());
			for (size_t j = 0; j < elements[e]->nodes() * pointDOFs.size(); j++) {
				column = pointDOFs.size() * elements[e]->node(j / pointDOFs.size()) + (j % pointDOFs.size());
				_K(row, column) = _apimesh.eMatrix(e)[i * elements[e]->nodes() * pointDOFs.size() + j];
			}
		}
	}

	SparseCSRMatrix<eslocal> csrK = _K;
	K[subdomain] = csrK;

	for (size_t p = 0; p < _mesh.parts(); p++) {
		const std::vector<eslocal> &l2c = _mesh.coordinates().localToCluster(p);
		f[p].resize(_mesh.coordinates().localSize(p) * pointDOFs.size(), 0);
		for (size_t i = 0; i < l2c.size() * pointDOFs.size(); i++) {
			f[p][i] = rhs[pointDOFs.size() * l2c[i / pointDOFs.size()] + i % pointDOFs.size()] / _mesh.nodes()[l2c[i / pointDOFs.size()]]->domains().size();
		}
	}
}
