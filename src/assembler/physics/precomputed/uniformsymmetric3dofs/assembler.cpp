
#include "assembler.h"

using namespace espreso;

std::vector<Property> UniformSymmetric3DOFs::elementDOFs;
std::vector<Property> UniformSymmetric3DOFs::faceDOFs;
std::vector<Property> UniformSymmetric3DOFs::edgeDOFs;
std::vector<Property> UniformSymmetric3DOFs::pointDOFs = { Property::DISPLACEMENT_X, Property::DISPLACEMENT_Y, Property::DISPLACEMENT_Z };
std::vector<Property> UniformSymmetric3DOFs::midPointDOFs = { Property::DISPLACEMENT_X, Property::DISPLACEMENT_Y, Property::DISPLACEMENT_Z };

void UniformSymmetric3DOFs::init()
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

void UniformSymmetric3DOFs::composeSubdomain(size_t subdomain)
{
	SparseVVPMatrix<eslocal> _K;
	eslocal nK = _mesh.coordinates().localSize(subdomain) * DOFs.size();
	_K.resize(nK, nK);

	const std::vector<eslocal> &parts = _mesh.getPartition();
	const std::vector<Element*> &elements = _mesh.elements();
	const std::vector<std::vector<double> > &matrices = _apimesh.getMatrices();

	size_t row, column;

	for (size_t e = parts[subdomain]; e < parts[subdomain + 1]; e++) {
		for (size_t i = 0; i < elements[e]->nodes() * DOFs.size(); i++) {
			row = DOFs.size() * elements[e]->node(i / DOFs.size()) + (i % DOFs.size());
			for (size_t j = 0; j < elements[e]->nodes() * DOFs.size(); j++) {
				column = DOFs.size() * elements[e]->node(j / DOFs.size()) + (j % DOFs.size());
				_K(row, column) = matrices[e][i * elements[e]->nodes() * DOFs.size() + j];
			}
		}
	}

	SparseCSRMatrix<eslocal> csrK = _K;
	K[subdomain] = csrK;

	for (size_t p = 0; p < _mesh.parts(); p++) {
		const std::vector<eslocal> &l2c = _mesh.coordinates().localToCluster(p);
		f[p].resize(_mesh.coordinates().localSize(p) * DOFs.size(), 0);
		for (size_t i = 0; i < l2c.size() * DOFs.size(); i++) {
			f[p][i] = rhs[DOFs.size() * l2c[i / DOFs.size()] + i % DOFs.size()] / _mesh.nodes()[l2c[i / DOFs.size()]]->domains().size();
		}
	}
}
