
#include "linear.h"

namespace assembler {

template <>
size_t Linear<BEM>::subdomains()
{
	return this->_surface.parts();
}

template <>
void Linear<BEM>::KMf(size_t part, bool dynamics)
{
	// TODO: temperature??
	const std::vector<eslocal> &partition = _surface.getPartition();
	const std::vector<mesh::Element*> &elements = _surface.getElements();

	DenseMatrix _K;
	eslocal nK = _surface.coordinates().localSize(part) * mesh::Point::size();
	eslocal eSize = partition[part + 1] - partition[part];
	_K.resize(nK, nK);
	std::vector<double> nodes(nK);
	std::vector<eslocal> elems(3 * eSize);

	for (size_t i = 0; i < _surface.coordinates().localSize(part); i++) {
		&nodes[i * mesh::Point::size()] << _surface.coordinates().get(i, part);
	}
	for (size_t i = partition[part], index = 0; i < partition[part + 1]; i++, index++) {
		for (size_t j = 0; j < elements[i]->size(); j++) {
			elems[3 * index + j] = elements[i]->node(j);
		}
	}

	bem4i::getLameSteklovPoincare(
			_K.values(),
			_surface.coordinates().localSize(part),
			&nodes[0],
			eSize,
			&elems[0],
			0.3,			// nu
			2.1e5,			// E
			3,				// order near
			4,				// order far
			false			// verbose
			);

	DenseMatrix tmp;
	size_t n = _K.rows();
    for (int i = 0; i < n / 3; i++) {
        for (int j = 0; j < n; j++) {
            tmp(3 * i + 0, j) = _K(0 * (n / 3) + i, j);
            tmp(3 * i + 1, j) = _K(1 * (n / 3) + i, j);
            tmp(3 * i + 2, j) = _K(2 * (n / 3) + i, j);
        }
    }

    for (int i = 0; i < n / 3; i++) {
        for (int j = 0; j < n; j++) {
            _K(j, 3 * i + 0) = tmp(j, 0 * (n / 3) + i);
            _K(j, 3 * i + 1) = tmp(j, 1 * (n / 3) + i);
            _K(j, 3 * i + 2) = tmp(j, 2 * (n / 3) + i);
        }
    }

	// TODO: make it direct
	SparseCSRMatrix<eslocal> csrK = _K;
	this->_K[part] = csrK;

	// TODO: fill M, f
	_f[part].clear();
	_f[part].resize(_K.rows(), 0);
}

template <>
void Linear<BEM>::saveResult()
{
	_surface.store(mesh::VTK_FULL, "mesh", _prim_solution, 0.95, 0.9);
}

template <>
void Linear<BEM>::initSolver()
{
	_lin_solver.init(
		_surface,
		_K,
		_globalB,
		_localB,
		_lambda_map_sub_B1,
		_lambda_map_sub_B0,
		_lambda_map_sub_clst,
		_B1_duplicity,
		_f,
		_vec_c,
		_surface.getFixPoints(),
		_neighClusters
	);
}

}
