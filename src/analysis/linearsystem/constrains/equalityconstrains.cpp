
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
	bool redundantLagrange = true;

	struct HasDirichlet {
		esint index = 0;
		const Vector_Distributed<Vector_Sparse, T> &dirichlet;

		HasDirichlet(const Vector_Distributed<Vector_Sparse, T> &dirichlet): dirichlet(dirichlet) {}

		bool operator()(esint dof) {
			while (index < dirichlet.cluster.nnz && dirichlet.cluster.indices[index] < dof) { ++index; }
			return (index < dirichlet.cluster.nnz && dirichlet.cluster.indices[index] == dof);
		}
	};

	eq.global = eq.nhalo = eq.paired = eq.nn = eq.local = 0;
	eq.domain.resize(feti.K.domains.size());

	std::vector<std::vector<esint> > COLS(feti.K.domains.size()), D2C(feti.K.domains.size());
	std::vector<std::vector<T> > VALS(feti.K.domains.size());
	std::vector<std::vector<esint, initless_allocator<esint> > > dpermutation(feti.K.domains.size());

	if (redundantLagrange && feti.K.decomposition->sharedDOFs.size()) {
		int local = feti.K.decomposition->neighbors.size(); // we need to sort them at the end
		HasDirichlet hasDirichlet(dirichlet);
		struct __map__ { int from, to, offset, neigh; };
		std::vector<__map__> lmap;

		auto map = feti.K.decomposition->dmap->cbegin();
		for (size_t n = 0, prev = 0; n < feti.K.decomposition->sharedDOFs.size(); prev = feti.K.decomposition->sharedDOFs[n++]) {
			map += feti.K.decomposition->sharedDOFs[n] - prev;
			if (hasDirichlet(feti.K.decomposition->sharedDOFs[n])) { continue; } // DOFs with Dirichlet are not glued

			for (auto di1 = map->begin(); di1 != map->end(); ++di1) {
				for (auto di2 = di1 + 1; di2 != map->end(); ++di2) {
					if (feti.K.decomposition->ismy(di1->domain) || feti.K.decomposition->ismy(di2->domain)) {

						lmap.push_back({ di1->domain, di2->domain, eq.paired, local });
						if (feti.K.decomposition->ismy(di1->domain)) {
							D2C[di1->domain - feti.K.decomposition->dbegin].push_back(eq.paired);
							COLS[di1->domain - feti.K.decomposition->dbegin].push_back(di1->index);
							VALS[di1->domain - feti.K.decomposition->dbegin].push_back(1);
						} else {
							lmap.back().neigh = feti.K.decomposition->noffset(di1->domain);
						}
						if (feti.K.decomposition->ismy(di2->domain)) {
							D2C[di2->domain - feti.K.decomposition->dbegin].push_back(eq.paired);
							COLS[di2->domain - feti.K.decomposition->dbegin].push_back(di2->index);
							VALS[di2->domain - feti.K.decomposition->dbegin].push_back(-1);
						} else {
							lmap.back().neigh = feti.K.decomposition->noffset(di2->domain);
						}
						++eq.paired;
					}
				}
			}
		}

		std::vector<esint, initless_allocator<esint> > permutation(lmap.size()), backpermutation(lmap.size());
		std::iota(permutation.begin(), permutation.end(), 0);
		std::sort(permutation.begin(), permutation.end(), [&] (esint i, esint j) {
			if (lmap[i].neigh != lmap[j].neigh) {
				return lmap[i].neigh < lmap[j].neigh;
			}
			if (lmap[i].from == lmap[j].from) {
				if (lmap[i].to == lmap[j].to) {
					return lmap[i].offset < lmap[j].offset;
				}
				return lmap[i].to < lmap[j].to;
			}
			return lmap[i].from < lmap[j].from;
		});

		eq.lmap.push_back({ lmap[*permutation.begin()].from, lmap[*permutation.begin()].to, 0, lmap[*permutation.begin()].neigh == local ? LMAP::LOCAL : lmap[*permutation.begin()].neigh });
		eq.nhalo = feti.K.decomposition->dbegin <= eq.lmap.back().from ? 0L : permutation.size();
		for (auto p = permutation.begin(); p != permutation.end(); ++p) {
			esint offset = backpermutation[*p] = p - permutation.begin();
			if (eq.lmap.back().from != lmap[*p].from || eq.lmap.back().to != lmap[*p].to) {
				eq.lmap.push_back({ lmap[*p].from, lmap[*p].to, offset, lmap[*p].neigh == local ? LMAP::LOCAL : lmap[*p].neigh });
				if (feti.K.decomposition->dbegin <= eq.lmap.back().from) {
					eq.nhalo = std::min(eq.nhalo, offset);
				}
			}
		}

		#pragma omp parallel for
		for (size_t d = 0; d < feti.K.domains.size(); ++d) {
			dpermutation[d].resize(COLS[d].size());
			std::iota(dpermutation[d].begin(), dpermutation[d].end(), 0);
			std::sort(dpermutation[d].begin(), dpermutation[d].end(), [&] (esint i, esint j) {
				return backpermutation[D2C[d][i]] < backpermutation[D2C[d][j]];
			});
			eq.domain[d].D2C.reserve(dpermutation[d].size());
			for (auto p = dpermutation[d].begin(); p != dpermutation[d].end(); ++p) {
				eq.domain[d].D2C.push_back(backpermutation[D2C[d][*p]]);
			}
		}
	}

	if (dirichlet.cluster.nnz) {
		auto map = feti.K.decomposition->dmap->cbegin();
		for (esint i = 0, prev = 0; i < dirichlet.cluster.nnz; prev = dirichlet.cluster.indices[i++]) {
			map += dirichlet.cluster.indices[i] - prev;
			for (auto di = map->begin(); di != map->end(); ++di) {
				if (feti.K.decomposition->ismy(di->domain)) {
					COLS[di->domain - feti.K.decomposition->dbegin].push_back(di->index);
					VALS[di->domain - feti.K.decomposition->dbegin].push_back(1);
					++eq.local;
				}
			}
		}
	}
	eq.c.resize(eq.paired + eq.local);
	std::vector<esint> doffset(feti.K.domains.size());
	for (esint d = 0, offset = eq.paired; d < (esint)feti.K.domains.size(); ++d) {
		if (COLS[d].size() != dpermutation[d].size()) {
			eq.lmap.push_back({ feti.K.decomposition->dbegin + d, feti.K.decomposition->dbegin + d, offset, LMAP::DIRICHLET });
			doffset[d] = offset;
			offset += COLS[d].size() - dpermutation[d].size();
		}
	}
	#pragma omp parallel for
	for (size_t d = 0; d < feti.K.domains.size(); ++d) {
		for (size_t i = dpermutation[d].size(), j = 0; i < COLS[d].size(); ++i) {
			eq.domain[d].D2C.push_back(doffset[d] + j++);
		}
	}

	#pragma omp parallel for
	for (size_t d = 0; d < feti.K.domains.size(); ++d) {
		eq.domain[d].B1.resize(COLS[d].size(), feti.K.domains[d].nrows, COLS[d].size());
		eq.domain[d].duplication.resize(COLS[d].size());
		math::set(eq.domain[d].duplication, T{.5});

		std::iota(eq.domain[d].B1.rows, eq.domain[d].B1.rows + COLS[d].size() + 1, 0); // B1 is indexed from 0
		for (size_t i = 0; i < dpermutation[d].size(); ++i) {
			eq.domain[d].B1.cols[i] = COLS[d][dpermutation[d][i]];
			eq.domain[d].B1.vals[i] = VALS[d][dpermutation[d][i]];
		}
		for (size_t i = dpermutation[d].size(); i < COLS[d].size(); ++i) {
			eq.domain[d].B1.cols[i] = COLS[d][i];
			eq.domain[d].B1.vals[i] = 1;
		}
	}

	eq.ordered.reserve(eq.lmap.size());
	for (size_t i = 0; i < eq.lmap.size(); ++i) {
		eq.ordered.push_back(i);
	}
	std::sort(eq.ordered.begin(), eq.ordered.end(), [&] (esint i, esint j) {
		if (eq.lmap[i].from == eq.lmap[j].from) {
			return eq.lmap[i].to < eq.lmap[j].to;
		}
		return eq.lmap[i].from < eq.lmap[j].from;
	});
}

template <typename T>
void EqualityConstrains<T>::update(step::Step &step, const Vector_Distributed<Vector_Sparse, T> &dirichlet)
{
	typename FETI<T>::EqualityConstraints &eq = feti.equalityConstraints;
	bool redundantLagrange = true;

	// TODO: store Dirichlet directly to 'c'
	math::set(eq.c, T{0});
	if (dirichlet.cluster.nnz) {
		std::vector<std::vector<T> > C(feti.K.domains.size());
		auto map = feti.K.decomposition->dmap->cbegin();
		for (esint i = 0, prev = 0; i < dirichlet.cluster.nnz; prev = dirichlet.cluster.indices[i++]) {
			map += dirichlet.cluster.indices[i] - prev;
			for (auto di = map->begin(); di != map->end(); ++di) {
				if (feti.K.decomposition->ismy(di->domain)) {
					C[di->domain - feti.K.decomposition->dbegin].push_back(dirichlet.cluster.vals[i]);
					++eq.local;
				}
			}
		}
		std::vector<esint> doffset(feti.K.domains.size());
		for (esint d = 0, offset = eq.paired; d < (esint)feti.K.domains.size(); ++d) {
			std::copy(C[d].begin(), C[d].end(), eq.c.vals + offset);
			offset += C[d].size();
		}
	}
}

template struct EqualityConstrains<double>;

}


