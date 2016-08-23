
#include "assembler.h"
//#include "esbem.h"

using namespace espreso;

void LameSteklovPoincare::composeSubdomain(size_t subdomain)
{
	const std::vector<eslocal> &partition = _mesh.getPartition();
	const std::vector<Element*> &elements = _mesh.getElements();

	DenseMatrix _K;
	eslocal nK = _mesh.coordinates().localSize(subdomain) * Point::size();
	eslocal eSize = partition[subdomain + 1] - partition[subdomain];
	_K.resize(nK, nK);
	std::vector<double> nodes(nK);
	std::vector<eslocal> elems(3 * eSize);

	for (size_t i = 0; i < _mesh.coordinates().localSize(subdomain); i++) {
		nodes[i * Point::size() + 0] = _mesh.coordinates().get(i, subdomain).x;
		nodes[i * Point::size() + 1] = _mesh.coordinates().get(i, subdomain).y;
		nodes[i * Point::size() + 2] = _mesh.coordinates().get(i, subdomain).z;
	}
	for (size_t i = partition[subdomain], index = 0; i < partition[subdomain + 1]; i++, index++) {
		for (size_t j = 0; j < elements[i]->size(); j++) {
			elems[3 * index + j] = elements[i]->node(j);
		}
	}
	ESINFO(GLOBAL_ERROR) << "missing BEM library";
/*
	bem4i::getLameSteklovPoincare<eslocal, double>(
			_K.values(),
			_mesh.coordinates().localSize(subdomain),
			&nodes[0],
			eSize,
			&elems[0],
			0.3,			// nu
			2.1e5,			// E
			3,				// order near
			4,				// order far
			false			// verbose
			);
*/
	DenseMatrix tmp = _K;
	eslocal n = _K.rows();
	for (eslocal i = 0; i < n / 3; i++) {
		for (eslocal j = 0; j < n; j++) {
			tmp(3 * i + 0, j) = _K(0 * (n / 3) + i, j);
			tmp(3 * i + 1, j) = _K(1 * (n / 3) + i, j);
			tmp(3 * i + 2, j) = _K(2 * (n / 3) + i, j);
		}
	}

	for (eslocal i = 0; i < n / 3; i++) {
		for (eslocal j = 0; j < n; j++) {
			_K(j, 3 * i + 0) = tmp(j, 0 * (n / 3) + i);
			_K(j, 3 * i + 1) = tmp(j, 1 * (n / 3) + i);
			_K(j, 3 * i + 2) = tmp(j, 2 * (n / 3) + i);
		}
	}

	// TODO: make it direct
	SparseCSRMatrix<eslocal> csrK = _K;
	K[subdomain] = csrK;

	f[subdomain].clear();
	f[subdomain].resize(_K.rows(), 0);
}
