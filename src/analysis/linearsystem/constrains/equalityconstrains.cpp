
#include "equalityconstrains.h"

#include "basis/containers/allocators.h"
#include "math/physics/vector_distributed.h"

#include <algorithm>
#include <numeric>

namespace espreso {

template <typename T>
void EqualityConstrains<T>::set(step::Step &step, const Vector_Distributed<Vector_Sparse, T> &dirichlet)
{
	typename FETI<T>::EqualityConstraints &eq = feti.equalityConstraints;
	eq.domain.resize(feti.K.domains.size());

	size_t maxMultiplicity = 2;
	std::vector<std::vector<esint> > COLS(feti.K.domains.size()), D2C(feti.K.domains.size());
	std::vector<std::vector<T> > VALS(feti.K.domains.size());

	struct __lambda__ { int dof, size; };
	std::vector<__lambda__> permutation;

	esint dindex = 0, dof = 0;
	for (auto dmap = feti.decomposition->dmap->cbegin(); dmap != feti.decomposition->dmap->cend(); ++dmap, ++dof) {
		while (dindex < dirichlet.cluster.nnz && dirichlet.cluster.indices[dindex] < dof) { ++dindex; }
		if (dindex < dirichlet.cluster.nnz && dirichlet.cluster.indices[dindex] == dof) {
			for (auto di = dmap->begin(); di != dmap->end(); ++di) {
				if (feti.decomposition->ismy(di->domain)) {
					COLS[di->domain - feti.decomposition->dbegin].push_back(di->index);
					VALS[di->domain - feti.decomposition->dbegin].push_back(1);
				}
			}
		} else {
			maxMultiplicity = std::max(maxMultiplicity, dmap->size());
			size_t i = 0;
			while (i < dmap->size() && !feti.decomposition->ismy(dmap->at(i).domain)) { ++i; }
			for (i = std::max(i, 1UL); i < dmap->size(); ++i) {
				permutation.push_back(__lambda__{ dof, (int)i + 1 });
			}
		}
	}

	std::sort(permutation.begin(), permutation.end(), [&] (const __lambda__ &i, const __lambda__ &j) {
		auto imap = feti.decomposition->dmap->cbegin() + i.dof;
		auto jmap = feti.decomposition->dmap->cbegin() + j.dof;
		int size = std::min(i.size, j.size);
		for (int k = 0; k < size; ++k) {
			if (imap->at(k).domain != jmap->at(k).domain) { return imap->at(k).domain < jmap->at(k).domain; }
		}
		if (i.size != j.size) {
			return i.size < j.size;
		} else {
			return i.dof < j.dof;
		}
	});

	Matrix_Dense<double> lambdas;
	lambdas.resize(maxMultiplicity - 1, maxMultiplicity);
	math::set(lambdas.nrows * lambdas.ncols, lambdas.vals, 1, T{0});
	for (esint r = 0, nc = 1; r < lambdas.nrows; ++r, ++nc) {
		double scale = std::sqrt(1 + (double)nc / (nc * nc));
		for (esint c = 0; c < nc; ++c) {
			lambdas.vals[r * maxMultiplicity + c] = scale / (nc + 1);
		}
		lambdas.vals[r * maxMultiplicity + nc] = -scale * nc / (nc + 1);
	}

	esint lambda = 0;
	for (size_t d = 0; d < feti.K.domains.size(); ++d) {
		for (size_t i = 0; i < COLS[d].size(); ++i) {
			D2C[d].push_back(lambda++);
		}
		eq.cmap.push_back(COLS[d].size());
		eq.cmap.push_back(1);
		eq.cmap.push_back(feti.decomposition->dbegin + d);
	}
	eq.dirichlet = eq.nhalo = lambda;

	for (size_t i = 0, prev = eq.cmap.size(); i < permutation.size(); ++i, ++lambda) {
		auto dmap = feti.decomposition->dmap->cbegin() + permutation[i].dof;
		size_t cbegin = eq.cmap.size();
		eq.cmap.push_back(1);
		eq.cmap.push_back(permutation[i].size);
		if (dmap->at(0).domain < feti.decomposition->dbegin) {
			eq.nhalo = lambda + 1;
		}
		for (int c = 0, r = permutation[i].size - 2; c < permutation[i].size; ++c) {
			eq.cmap.push_back(dmap->at(c).domain);
			if (feti.decomposition->ismy(dmap->at(c).domain)) {
				D2C [dmap->at(c).domain - feti.decomposition->dbegin].push_back(lambda);
				COLS[dmap->at(c).domain - feti.decomposition->dbegin].push_back(dmap->at(c).index);
				VALS[dmap->at(c).domain - feti.decomposition->dbegin].push_back(lambdas.vals[r * maxMultiplicity + c]);
			}
		}
		if (eq.cmap[prev + 1] == eq.cmap[cbegin + 1]) { // domains sizes
			if (std::equal(eq.cmap.cbegin() + cbegin + 1, eq.cmap.cend(), eq.cmap.cbegin() + prev + 1)) {
				if (prev != cbegin) {
					++eq.cmap[prev];
					eq.cmap.resize(cbegin);
				}
			}
		}
	}
	eq.size = lambda;

	#pragma omp parallel for
	for (size_t d = 0; d < feti.K.domains.size(); ++d) {
		eq.domain[d].B1.resize(COLS[d].size(), feti.K.domains[d].nrows, COLS[d].size());
		std::iota(eq.domain[d].B1.rows, eq.domain[d].B1.rows + COLS[d].size() + 1, 0); // B1 is indexed from 0
		std::copy(COLS[d].begin(), COLS[d].end(), eq.domain[d].B1.cols);
		std::copy(VALS[d].begin(), VALS[d].end(), eq.domain[d].B1.vals);
		eq.domain[d].D2C = D2C[d];
	}

	eq.c.resize(eq.size);
//
//	for (size_t i = 0; i < eq.cmap.size(); ) {
//		printf("%dx:", eq.cmap[i]);
//		for (esint d = 0; d < eq.cmap[i + 1]; ++d) {
//			printf(" %d", eq.cmap[i + 2 + d]);
//		}
//		printf("\n");
//		i += eq.cmap[i + 1] + 2;
//	}
}

template <typename T>
void EqualityConstrains<T>::update(step::Step &step, const Vector_Distributed<Vector_Sparse, T> &dirichlet)
{
	typename FETI<T>::EqualityConstraints &eq = feti.equalityConstraints;

	math::set(eq.c, T{0});
	std::vector<esint> dindex(feti.K.domains.size());

	for (size_t i = 0, index = 0; i < eq.cmap.size(); ) {
		dindex[eq.cmap[i + 2] - feti.decomposition->dbegin] = index;
		index += eq.cmap[i];
		if ((esint)index == eq.dirichlet) {
			break;
		}
		i += eq.cmap[i + 1] + 2;
	}

	if (dirichlet.cluster.nnz) {
		auto map = feti.decomposition->dmap->cbegin();
		for (esint i = 0, prev = 0; i < dirichlet.cluster.nnz; prev = dirichlet.cluster.indices[i++]) {
			map += dirichlet.cluster.indices[i] - prev;
			for (auto di = map->begin(); di != map->end(); ++di) {
				if (feti.decomposition->ismy(di->domain)) {
					eq.c.vals[dindex[di->domain - feti.decomposition->dbegin]++] = dirichlet.cluster.vals[i];
				}
			}
		}
	}
}

template struct EqualityConstrains<double>;

}


