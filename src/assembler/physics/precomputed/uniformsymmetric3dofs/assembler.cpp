
#include "assembler.h"

using namespace espreso;

void UniformSymmetric3DOFs::composeSubdomain(size_t subdomain)
{
	SparseVVPMatrix<eslocal> _K;
	eslocal nK = _mesh.coordinates().localSize(subdomain) * DOFs;
	_K.resize(nK, nK);

	const std::vector<eslocal> &parts = _mesh.getPartition();
	const std::vector<Element*> &elements = _mesh.getElements();
	const std::vector<std::vector<double> > &matrices = _apimesh.getMatrices();

	size_t row, column;

	for (size_t e = parts[subdomain]; e < parts[subdomain + 1]; e++) {
		for (size_t i = 0; i < elements[e]->size() * DOFs; i++) {
			row = DOFs * elements[e]->node(i / DOFs) + (i % DOFs);
			for (size_t j = 0; j < elements[e]->size() * DOFs; j++) {
				column = DOFs * elements[e]->node(j / DOFs) + (j % DOFs);
				_K(row, column) = matrices[e][i * elements[e]->size() * DOFs + j];
			}
		}
	}

	SparseCSRMatrix<eslocal> csrK = _K;
	K[subdomain] = csrK;
	K[subdomain].mtype = SparseMatrix::MatrixType::REAL_UNSYMMETRIC;

	for (size_t p = 0; p < _mesh.parts(); p++) {
		const std::vector<eslocal> &l2c = _mesh.coordinates().localToCluster(p);
		const Boundaries &boundaries = _mesh.subdomainBoundaries();
		f[p].resize(_mesh.coordinates().localSize(p) * DOFs, 0);
		for (size_t i = 0; i < l2c.size() * DOFs; i++) {
			f[p][i] = rhs[DOFs * l2c[i / DOFs] + i % DOFs] / boundaries[l2c[i / DOFs]].size();
		}
	}
}
