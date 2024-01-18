
#include "equalityconstrains.h"

#include "basis/containers/allocators.h"
#include "math/physics/vector_distributed.h"

#include <algorithm>
#include <numeric>

namespace espreso {

template struct EqualityConstrains<double>;

template <typename T>
void EqualityConstrains<T>::set(step::Step &step, const Vector_Distributed<Vector_Sparse, T> &dirichlet)
{
	feti.B1.resize(feti.K.size());
	feti.D2C.resize(feti.K.size());

	size_t maxMultiplicity = 2;
	std::vector<std::vector<esint> > COLS(feti.K.size()), D2C(feti.K.size());
	std::vector<std::vector<T> > VALS(feti.K.size());

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
	for (size_t d = 0; d < feti.K.size(); ++d) {
		for (size_t i = 0; i < COLS[d].size(); ++i) {
			D2C[d].push_back(lambda++);
		}
		feti.lambdas.cmap.push_back(COLS[d].size());
		feti.lambdas.cmap.push_back(1);
		feti.lambdas.cmap.push_back(feti.decomposition->dbegin + d);
	}
	feti.lambdas.dirichlet = feti.lambdas.nhalo = lambda;

	for (size_t i = 0, prev = feti.lambdas.cmap.size(); i < permutation.size(); ++i, ++lambda) {
		auto dmap = feti.decomposition->dmap->cbegin() + permutation[i].dof;
		size_t cbegin = feti.lambdas.cmap.size();
		feti.lambdas.cmap.push_back(1);
		feti.lambdas.cmap.push_back(permutation[i].size);
		if (dmap->at(0).domain < feti.decomposition->dbegin) {
			feti.lambdas.nhalo = lambda + 1;
		}
		for (int c = 0, r = permutation[i].size - 2; c < permutation[i].size; ++c) {
			feti.lambdas.cmap.push_back(dmap->at(c).domain);
			if (feti.decomposition->ismy(dmap->at(c).domain)) {
				D2C [dmap->at(c).domain - feti.decomposition->dbegin].push_back(lambda);
				COLS[dmap->at(c).domain - feti.decomposition->dbegin].push_back(dmap->at(c).index);
				VALS[dmap->at(c).domain - feti.decomposition->dbegin].push_back(lambdas.vals[r * maxMultiplicity + c]);
			}
		}
		if (feti.lambdas.cmap[prev + 1] == feti.lambdas.cmap[cbegin + 1]) { // domains sizes
			if (std::equal(feti.lambdas.cmap.cbegin() + cbegin + 1, feti.lambdas.cmap.cend(), feti.lambdas.cmap.cbegin() + prev + 1)) {
				if (prev != cbegin) {
					++feti.lambdas.cmap[prev];
					feti.lambdas.cmap.resize(cbegin);
				}
			}
		}
	}
	feti.lambdas.size = lambda;

	#pragma omp parallel for
	for (size_t d = 0; d < feti.K.size(); ++d) {
		feti.B1[d].resize(COLS[d].size(), feti.K[d].nrows, COLS[d].size());
		std::iota(feti.B1[d].rows, feti.B1[d].rows + COLS[d].size() + 1, 0); // B1 is indexed from 0
		std::copy(COLS[d].begin(), COLS[d].end(), feti.B1[d].cols);
		std::copy(VALS[d].begin(), VALS[d].end(), feti.B1[d].vals);
		feti.D2C[d] = D2C[d];
	}

	feti.c.resize(feti.lambdas.size);
//
//	for (size_t i = 0; i < feti.lambdas.cmap.size(); ) {
//		printf("%dx:", feti.lambdas.cmap[i]);
//		for (esint d = 0; d < feti.lambdas.cmap[i + 1]; ++d) {
//			printf(" %d", feti.lambdas.cmap[i + 2 + d]);
//		}
//		printf("\n");
//		i += feti.lambdas.cmap[i + 1] + 2;
//	}
}

template <typename T>
void EqualityConstrains<T>::update(step::Step &step, const Vector_Distributed<Vector_Sparse, T> &dirichlet)
{
	math::set(feti.c, T{0});
	std::vector<esint> dindex(feti.K.size());

	for (size_t i = 0, index = 0; i < feti.lambdas.cmap.size(); ) {
		dindex[feti.lambdas.cmap[i + 2] - feti.decomposition->dbegin] = index;
		index += feti.lambdas.cmap[i];
		if ((esint)index == feti.lambdas.dirichlet) {
			break;
		}
		i += feti.lambdas.cmap[i + 1] + 2;
	}

	if (dirichlet.cluster.nnz) {
		auto map = feti.decomposition->dmap->cbegin();
		for (esint i = 0, prev = 0; i < dirichlet.cluster.nnz; prev = dirichlet.cluster.indices[i++]) {
			map += dirichlet.cluster.indices[i] - prev;
			for (auto di = map->begin(); di != map->end(); ++di) {
				if (feti.decomposition->ismy(di->domain)) {
					feti.c.vals[dindex[di->domain - feti.decomposition->dbegin]++] = dirichlet.cluster.vals[i];
				}
			}
		}
	}
}

}


