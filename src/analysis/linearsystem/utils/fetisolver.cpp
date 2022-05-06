
#include "analysis/linearsystem/fetisystem.h"
#include "basis/containers/allocators.h"
#include "basis/utilities/utils.h"
#include "esinfo/meshinfo.h"
#include "math/physics/matrix_feti.h"
#include "math/primitives/matrix_dense.h"
#include "mesh/store/domainstore.h"

#include <algorithm>
#include <numeric>

namespace espreso {

template <typename T>
void _composeEqualityConstraints(const Matrix_FETI<Matrix_CSR, T> &K, const Vector_Distributed<Vector_Sparse, T> &dirichlet, typename FETI<T>::EqualityConstraints &eq, bool redundantLagrange)
{
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
	eq.domain.resize(K.domains.size());

	std::vector<std::vector<esint> > COLS(K.domains.size()), D2C(K.domains.size());
	std::vector<std::vector<T> > VALS(K.domains.size());
	std::vector<std::vector<esint, initless_allocator<esint> > > dpermutation(K.domains.size());

	if (redundantLagrange && K.decomposition->sharedDOFs.size()) {
		int local = K.decomposition->neighbors.size(); // we need to sort them at the end
		HasDirichlet hasDirichlet(dirichlet);
		struct __map__ { int from, to, offset, neigh; };
		std::vector<__map__> lmap;

		auto map = K.decomposition->dmap->cbegin();
		for (size_t n = 0, prev = 0; n < K.decomposition->sharedDOFs.size(); prev = K.decomposition->sharedDOFs[n++]) {
			map += K.decomposition->sharedDOFs[n] - prev;
			if (hasDirichlet(K.decomposition->sharedDOFs[n])) { continue; } // DOFs with Dirichlet are not glued

			for (auto di1 = map->begin(); di1 != map->end(); ++di1) {
				for (auto di2 = di1 + 1; di2 != map->end(); ++di2) {
					if (K.decomposition->ismy(di1->domain) || K.decomposition->ismy(di2->domain)) {

						lmap.push_back({ di1->domain, di2->domain, eq.paired, local });
						if (K.decomposition->ismy(di1->domain)) {
							D2C[di1->domain - K.decomposition->dbegin].push_back(eq.paired);
							COLS[di1->domain - K.decomposition->dbegin].push_back(di1->index);
							VALS[di1->domain - K.decomposition->dbegin].push_back(1);
						} else {
							lmap.back().neigh = K.decomposition->noffset(di1->domain);
						}
						if (K.decomposition->ismy(di2->domain)) {
							D2C[di2->domain - K.decomposition->dbegin].push_back(eq.paired);
							COLS[di2->domain - K.decomposition->dbegin].push_back(di2->index);
							VALS[di2->domain - K.decomposition->dbegin].push_back(-1);
						} else {
							lmap.back().neigh = K.decomposition->noffset(di2->domain);
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
		eq.nhalo = K.decomposition->dbegin <= eq.lmap.back().from ? 0L : permutation.size();
		for (auto p = permutation.begin(); p != permutation.end(); ++p) {
			esint offset = backpermutation[*p] = p - permutation.begin();
			if (eq.lmap.back().from != lmap[*p].from || eq.lmap.back().to != lmap[*p].to) {
				eq.lmap.push_back({ lmap[*p].from, lmap[*p].to, offset, lmap[*p].neigh == local ? LMAP::LOCAL : lmap[*p].neigh });
				if (K.decomposition->dbegin <= eq.lmap.back().from) {
					eq.nhalo = std::min(eq.nhalo, offset);
				}
			}
		}

		#pragma omp parallel for
		for (size_t d = 0; d < K.domains.size(); ++d) {
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
		auto map = K.decomposition->dmap->cbegin();
		for (esint i = 0, prev = 0; i < dirichlet.cluster.nnz; prev = dirichlet.cluster.indices[i++]) {
			map += dirichlet.cluster.indices[i] - prev;
			for (auto di = map->begin(); di != map->end(); ++di) {
				if (K.decomposition->ismy(di->domain)) {
					COLS[di->domain - K.decomposition->dbegin].push_back(di->index);
					VALS[di->domain - K.decomposition->dbegin].push_back(1);
					++eq.local;
				}
			}
		}
	}
	eq.c.resize(eq.paired + eq.local);
	std::vector<esint> doffset(K.domains.size());
	for (esint d = 0, offset = eq.paired; d < (esint)K.domains.size(); ++d) {
		if (COLS[d].size() != dpermutation[d].size()) {
			eq.lmap.push_back({ K.decomposition->dbegin + d, K.decomposition->dbegin + d, offset, LMAP::DIRICHLET });
			doffset[d] = offset;
			offset += COLS[d].size() - dpermutation[d].size();
		}
	}
	#pragma omp parallel for
	for (size_t d = 0; d < K.domains.size(); ++d) {
		for (size_t i = dpermutation[d].size(), j = 0; i < COLS[d].size(); ++i) {
			eq.domain[d].D2C.push_back(doffset[d] + j++);
		}
	}

	#pragma omp parallel for
	for (size_t d = 0; d < K.domains.size(); ++d) {
		eq.domain[d].B1.resize(COLS[d].size(), K.domains[d].nrows, COLS[d].size());
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
void _evaluateEqualityConstraints(const Matrix_FETI<Matrix_CSR, T> &K, const Vector_Distributed<Vector_Sparse, T> &dirichlet, typename FETI<T>::EqualityConstraints &eq, bool redundantLagrange)
{
	// TODO: store Dirichlet directly to 'c'
	math::set(eq.c, T{0});
	if (dirichlet.cluster.nnz) {
		std::vector<std::vector<T> > C(K.domains.size());
		auto map = K.decomposition->dmap->cbegin();
		for (esint i = 0, prev = 0; i < dirichlet.cluster.nnz; prev = dirichlet.cluster.indices[i++]) {
			map += dirichlet.cluster.indices[i] - prev;
			for (auto di = map->begin(); di != map->end(); ++di) {
				if (K.decomposition->ismy(di->domain)) {
					C[di->domain - K.decomposition->dbegin].push_back(dirichlet.cluster.vals[i]);
					++eq.local;
				}
			}
		}
		std::vector<esint> doffset(K.domains.size());
		for (esint d = 0, offset = eq.paired; d < (esint)K.domains.size(); ++d) {
			std::copy(C[d].begin(), C[d].end(), eq.c.vals + offset);
			offset += C[d].size();
		}
	}
}

void composeHeatTransferKernel(const Matrix_CSR<double> &K, Matrix_Dense<double> &R, Matrix_CSR<double> &RegMat)
{
	R.resize(K.nrows, 1);
	R.type = Matrix_Type::REAL_NONSYMMETRIC;
	R.shape = Matrix_Shape::FULL;

	RegMat.resize(K.nrows, K.ncols, 1);
	RegMat.type = K.type;
	RegMat.shape = K.shape;

	RegMat.rows[0] = RegMat.cols[0] = Indexing::CSR;
	std::fill(RegMat.rows + 1, RegMat.rows + RegMat.nrows + 1, RegMat.rows[0] + 1);
}

void evaluateHeatTransferKernel(const Matrix_CSR<double> &K, Matrix_Dense<double> &R, Matrix_CSR<double> &RegMat)
{
	RegMat.vals[0] = math::getDiagonalMax(K);
	math::set(R, 1.0 / std::sqrt(K.nrows));

//	if (solver.configuration.method != FETIConfiguration::METHOD::HYBRID_FETI) { // N1 orthogonal for whole cluster
//		esint crows = 0;
//		for (esint dd = 0; dd < info::mesh->domains->size; dd++) {
//			if (info::mesh->domains->cluster[d] == info::mesh->domains->cluster[dd]) {
//				crows += solver.A.domains[dd].nrows;
//			}
//		}
//		math::fill(solver.N1.domains[d], 1.0 / std::sqrt(crows));
//	} else {
//		math::fill(solver.N1.domains[d], 1.0 / std::sqrt(solver.A.domains[d].nrows));
//	}
}

template <> void composeEqualityConstraints(const Matrix_FETI<Matrix_CSR, double> &K, const Vector_Distributed<Vector_Sparse, double> &dirichlet, FETI<double>::EqualityConstraints &equalityConstraints, bool redundantLagrange) { _composeEqualityConstraints(K, dirichlet, equalityConstraints, redundantLagrange); }
template <> void composeEqualityConstraints(const Matrix_FETI<Matrix_CSR, std::complex<double> > &K, const Vector_Distributed<Vector_Sparse, std::complex<double> > &dirichlet, FETI<std::complex<double> >::EqualityConstraints &equalityConstraints, bool redundantLagrange) { _composeEqualityConstraints(K, dirichlet, equalityConstraints, redundantLagrange); }

template <> void evaluateEqualityConstraints(const Matrix_FETI<Matrix_CSR, double> &K, const Vector_Distributed<Vector_Sparse, double> &dirichlet, FETI<double>::EqualityConstraints &equalityConstraints, bool redundantLagrange) { _evaluateEqualityConstraints(K, dirichlet, equalityConstraints, redundantLagrange); }
template <> void evaluateEqualityConstraints(const Matrix_FETI<Matrix_CSR, std::complex<double> > &K, const Vector_Distributed<Vector_Sparse, std::complex<double> > &dirichlet, FETI<std::complex<double> >::EqualityConstraints &equalityConstraints, bool redundantLagrange) { _evaluateEqualityConstraints(K, dirichlet, equalityConstraints, redundantLagrange); }

}
