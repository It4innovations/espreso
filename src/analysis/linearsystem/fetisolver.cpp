
#include "fetisolver.h"

#include "basis/containers/allocators.h"
#include "analysis/physics/heat.steadystate.linear.h"
#include "analysis/physics/heat.steadystate.nonlinear.h"

#include <numeric>

namespace espreso {

template struct FETILinearSystemSolver<double, HeatSteadyStateLinear>;
template struct FETILinearSystemSolver<double, HeatSteadyStateNonLinear>;

template <typename T>
static void _setEqualityConstraints(const Matrix_FETI<Matrix_CSR, T> &K, const Vector_Distributed<Vector_Sparse, T> &dirichlet, typename FETI<T>::EqualityConstraints &eq, bool redundantLagrange)
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
static void _updateEqualityConstraints(const Matrix_FETI<Matrix_CSR, T> &K, const Vector_Distributed<Vector_Sparse, T> &dirichlet, typename FETI<T>::EqualityConstraints &eq, bool redundantLagrange)
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

template <typename T, class Physics>
void FETILinearSystemSolver<T, Physics>::setEqualityConstraints(step::Step &step)
{
	_setEqualityConstraints(feti.K, dirichlet, feti.equalityConstraints, true);
}

template <typename T, class Physics>
void FETILinearSystemSolver<T, Physics>::updateEqualityConstraints(step::Step &step)
{
	_updateEqualityConstraints(feti.K, dirichlet, feti.equalityConstraints, true);
	if (info::ecf->output.print_matrices) {
		eslog::storedata(" STORE: system/{B1, B1c, B1Duplication, D2C, LMAP}\n");
		math::store(feti.equalityConstraints.c, utils::filename(utils::debugDirectory(step) + "/system", "B1c").c_str());
		for (size_t d = 0; d < feti.equalityConstraints.domain.size(); ++d) {
			math::store(feti.equalityConstraints.domain[d].B1, utils::filename(utils::debugDirectory(step) + "/system", "B1" + std::to_string(d)).c_str());
			math::store(feti.equalityConstraints.domain[d].duplication, utils::filename(utils::debugDirectory(step) + "/system", "B1Duplication" + std::to_string(d)).c_str());
			math::store(feti.equalityConstraints.domain[d].D2C, utils::filename(utils::debugDirectory(step) + "/system", "D2C" + std::to_string(d)).c_str());
			math::store(feti.equalityConstraints.lmap, utils::filename(utils::debugDirectory(step) + "/system", "LMAP").c_str());
		}
	}
}

template <typename T, class Physics> struct AnalyticKernel {
	static void set(const Matrix_CSR<double> &K, Matrix_Dense<double> &R, Matrix_CSR<double> &RegMat)
	{

	}

	static void update(const Matrix_CSR<double> &K, Matrix_Dense<double> &R, Matrix_CSR<double> &RegMat)
	{

	}
};

template <typename T>
struct AnalyticKernel<T, HeatSteadyStateLinear>
{
	static void set(const Matrix_CSR<double> &K, Matrix_Dense<double> &R, Matrix_CSR<double> &RegMat)
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

	static void update(const Matrix_CSR<double> &K, Matrix_Dense<double> &R, Matrix_CSR<double> &RegMat)
	{
		RegMat.vals[0] = math::getDiagonalMax(K);
		math::set(R, 1.0 / std::sqrt(K.nrows));
	}
};

template <typename T>
struct AnalyticKernel<T, HeatSteadyStateNonLinear>
{
	static void set(const Matrix_CSR<double> &K, Matrix_Dense<double> &R, Matrix_CSR<double> &RegMat)
	{
		AnalyticKernel<T, HeatSteadyStateLinear>::set(K, R, RegMat);
	}

	static void update(const Matrix_CSR<double> &K, Matrix_Dense<double> &R, Matrix_CSR<double> &RegMat)
	{
		AnalyticKernel<T, HeatSteadyStateLinear>::update(K, R, RegMat);
	}
};

template <typename T, class Physics>
void FETILinearSystemSolver<T, Physics>::setKernel(step::Step &step)
{
	eslog::info(" = REGULARIZATION                                                                   ANALYTIC = \n");
	feti.regularization.R1.domains.resize(feti.K.domains.size());
	feti.regularization.R2.domains.resize(feti.K.domains.size());
	feti.regularization.RegMat.domains.resize(feti.K.domains.size());

	feti.regularization.RegMat.type = feti.K.type;
	feti.regularization.RegMat.shape = feti.K.shape;
	#pragma omp parallel for
	for (size_t d = 0; d < feti.K.domains.size(); ++d) {
		// TODO: some domains can be without kernel
		AnalyticKernel<T, Physics>::set(feti.K.domains[d], feti.regularization.R1.domains[d], feti.regularization.RegMat.domains[d]);
//		switch (feti.configuration.regularization) {
//		case FETIConfiguration::REGULARIZATION::ANALYTIC: AnalyticKernel<T, Physics>::set(K.domains[d], regularization.R1.domains[d], regularization.RegMat.domains[d]); break;
//		case FETIConfiguration::REGULARIZATION::ALGEBRAIC: break;
//		}
	}
}

template <typename T, class Physics>
void FETILinearSystemSolver<T, Physics>::updateKernel(step::Step &step)
{
	#pragma omp parallel for
	for (size_t d = 0; d < feti.K.domains.size(); ++d) {
		AnalyticKernel<T, Physics>::update(feti.K.domains[d], feti.regularization.R1.domains[d], feti.regularization.RegMat.domains[d]);
//		switch (feti.configuration.regularization) {
//		case FETIConfiguration::REGULARIZATION::ANALYTIC: AnalyticKernel<T, Physics>::update(K.domains[d], regularization.R1.domains[d], regularization.RegMat.domains[d]); break;
//		case FETIConfiguration::REGULARIZATION::ALGEBRAIC: break;
//		}
	}
	if (info::ecf->output.print_matrices) {
		eslog::storedata(" STORE: system/{R, RegMat}\n");
		math::store(feti.regularization.R1, utils::filename(utils::debugDirectory(step) + "/system", "R").c_str());
		math::store(feti.regularization.RegMat, utils::filename(utils::debugDirectory(step) + "/system", "RegMat").c_str());
	}
}

}

