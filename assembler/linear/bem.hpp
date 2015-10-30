
#include "linear.h"

namespace assembler {

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
	this->_K[part] = csrK;

	// TODO: fill M, f
	_f[part].clear();
	_f[part].resize(_K.rows(), 0);
}

template <>
void Linear<BEM>::RHS()
{
	const std::map<eslocal, double> &forces_x = this->_mesh.coordinates().property(mesh::FORCES_X).values();
	const std::map<eslocal, double> &forces_y = this->_mesh.coordinates().property(mesh::FORCES_Y).values();
	const std::map<eslocal, double> &forces_z = this->_mesh.coordinates().property(mesh::FORCES_Z).values();

	for (size_t p = 0; p < this->_mesh.parts(); p++) {
		const std::vector<eslocal> &l2g = this->_surface.coordinates().localToCluster(p);
		for (eslocal i = 0; i < l2g.size(); i++) {
			if (forces_x.find(l2g[i]) != forces_x.end()) {
				_f[p][3 * i + 0] = forces_x.at(l2g[i]);
			}
			if (forces_y.find(l2g[i]) != forces_y.end()) {
				_f[p][3 * i + 1] = forces_y.at(l2g[i]);
			}
			if (forces_z.find(l2g[i]) != forces_z.end()) {
				_f[p][3 * i + 2] = forces_z.at(l2g[i]);
			}
		}
	}
}

template <>
void Linear<BEM>::saveResult()
{
	_surface.store(mesh::VTK_FULL, "surface", _prim_solution, 0.95, 0.9);
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
