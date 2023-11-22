

#include "tfetisymmetric.h"
#include "basis/containers/serializededata.h"
#include "basis/utilities/sysutils.h"
#include "basis/utilities/utils.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/envinfo.h"
#include "esinfo/eslog.hpp"

#include <vector>

namespace espreso {

template<typename T>
OrthogonalTFETISymmetric<T>::OrthogonalTFETISymmetric(FETI<T> &feti)
: Projector<T>(feti)
{
	e.resize();
	Gx.resize();
	iGGtGx.resize(feti.sinfo.R1size);

	_setG();
	if (feti.equalityConstraints.global) {
		_setDenseGGt();
	} else {
		_setSparseGGt();
	}
}

template<typename T>
OrthogonalTFETISymmetric<T>::~OrthogonalTFETISymmetric()
{

}

template<typename T>
void OrthogonalTFETISymmetric<T>::info()
{
	esint nnz = 2 * (GGt.nnz - GGt.nrows) + GGt.nrows;

	eslog::info(" = ORTHOGONAL PROJECTOR PROPERTIES                                                           = \n");
	eslog::info(" =   GGT ROWS                                                                      %9d = \n", GGtSolver.rows);
	eslog::info(" =   GGT NNZ                                                                       %9d = \n", nnz);
	eslog::info(" =   GGT FACTORS NNZ                                                               %9d = \n", GGtSolver.nnzL);
	if (feti.configuration.exhaustive_info) {
		// PPt = eye
	}
	eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
}

template<typename T>
void OrthogonalTFETISymmetric<T>::update(const step::Step &step)
{
	const typename FETI<T>::Regularization &R = feti.regularization;

	#pragma omp parallel for
	for (size_t d = 0; d < R.R1.domains.size(); ++d) {
		Vector_Dense<T> _e;
		_e.size = R.R1.domains[d].ncols;
		_e.vals = e.vals + Roffset[d] + Vector_Kernel<T>::offset;
		Matrix_Dense<T> _Rt;
		_Rt.nrows = R.R1.domains[d].ncols;
		_Rt.ncols = R.R1.domains[d].nrows;
		_Rt.vals = R.R1.domains[d].vals;
		math::apply(_e, T{1}, _Rt, T{0}, feti.f.domains[d]);
	}
	e.synchronize();
	eslog::checkpointln("FETI: COMPUTE DUAL RHS [e]");

	_updateG();
	if (feti.equalityConstraints.global) {
		_updateDenseGGt();
	} else {
		_updateSparseGGt();
	}

	_print(step);
}

template<typename T>
void OrthogonalTFETISymmetric<T>::apply(const Vector_Dual<T> &x, Vector_Dual<T> &y)
{
	x.copyToWithoutHalo(y);
	_applyG(x, Gx);
	_applyInvGGt(Gx, iGGtGx);
	_applyGt(iGGtGx, T{-1}, y);
}

template<typename T>
void OrthogonalTFETISymmetric<T>::applyGtInvGGt(const Vector_Kernel<T> &x, Vector_Dual<T> &y)
{
	math::set(y, T{0});
	_applyInvGGt(x, iGGtGx);
	_applyGt(iGGtGx, T{-1}, y);
}

template<typename T>
void OrthogonalTFETISymmetric<T>::applyRInvGGtG(const Vector_Dual<T> &x, Vector_FETI<Vector_Dense, T> &y)
{
	_applyG(x, Gx);
	_applyInvGGt(Gx, iGGtGx);
	_applyR(iGGtGx, y);
}

template<typename T>
void OrthogonalTFETISymmetric<T>::_applyG(const Vector_Dual<T> &in, Vector_Kernel<T> &out)
{
	#pragma omp parallel for
	for (int t = 0; t < info::env::threads; ++t) {
		for (size_t r = Vector_Kernel<T>::distribution[t]; r < Vector_Kernel<T>::distribution[t + 1]; ++r) {
			out.vals[r + Vector_Kernel<T>::offset] = T{0};
			for (esint c = G.rows[r]; c < G.rows[r + 1]; ++c) {
				out.vals[r + Vector_Kernel<T>::offset] += G.vals[c] * in.vals[G.cols[c]];
			}
		}
	}
	out.synchronize();
}

template<typename T>
void OrthogonalTFETISymmetric<T>::_applyInvGGt(const Vector_Kernel<T> &in, Vector_Dense<T> &out)
{
	#pragma omp parallel for
	for (int t = 0; t < info::env::threads; ++t) {
		Matrix_Dense<T> a;
		Vector_Dense<T> y;
		a.ncols = invGGt.ncols;
		a.nrows = y.size = Vector_Kernel<T>::distribution[t + 1] - Vector_Kernel<T>::distribution[t];

		a.vals = invGGt.vals + invGGt.ncols * Vector_Kernel<T>::distribution[t];
		y.vals = out.vals + Vector_Kernel<T>::distribution[t];

		math::apply(y, T{1}, a, T{0}, in);
	}
}

// TODO: threaded implementation: utilize properties of B of gluing two domains
template<typename T>
void OrthogonalTFETISymmetric<T>::_applyGt(const Vector_Dense<T> &in, const T &alpha, Vector_Dual<T> &out)
{
	for (esint r = 0; r < G.nrows; ++r) {
		for (esint c = G.rows[r]; c < G.rows[r + 1]; ++c) {
			out.vals[G.cols[c]] += alpha * G.vals[c] * in.vals[r];
		}
	}
	out.synchronize();
}

template<typename T>
void OrthogonalTFETISymmetric<T>::_applyR(const Vector_Dense<T> &in, Vector_FETI<Vector_Dense, T> &out)
{
	#pragma omp parallel for
	for (size_t d = 0; d < out.domains.size(); ++d) {
		Vector_Dense<T> y;
		y.size = feti.regularization.R1.domains[d].ncols;
		y.vals = in.vals + Roffset[d];

		math::apply(out.domains[d], T{1}, feti.regularization.R1.domains[d], T{0}, y);
	}
}

struct LMAPSize {
	LMAPSize(std::vector<LMAP>::const_iterator end, esint lambdas): end(end), lambdas(lambdas) { }

	int operator()(std::vector<LMAP>::const_iterator it) {
		std::vector<LMAP>::const_iterator next = it + 1;
		if (next != end) {
			return next->offset - it->offset;
		}
		return lambdas - it->offset;
	}

	std::vector<LMAP>::const_iterator end;
	esint lambdas;
};

template<typename T>
void OrthogonalTFETISymmetric<T>::_setG()
{
	// G is stored with 0-based in indexing
	const Matrix_FETI<Matrix_CSR, T> &K = feti.K;
	const typename FETI<T>::Regularization &R = feti.regularization;
	const typename FETI<T>::EqualityConstraints &L = feti.equalityConstraints;

	esint nrows = 0;
	Roffset.resize(K.domains.size());
	for (size_t d = 0; d < K.domains.size(); ++d) {
		Roffset[d] = nrows;
		nrows += R.R1.domains[d].ncols;
	}
	G.resize(nrows, feti.sinfo.lambdasLocal, 0);
	math::set(G.nrows + 1, G.rows, 1, 0);

	LMAPSize lsize(L.lmap.cend(), feti.sinfo.lambdasLocal);
	for (auto map = L.lmap.cbegin(); map != L.lmap.cend(); ++map) {
		auto add = [&] (esint domain) {
			domain -= K.decomposition->dbegin;
			esint kernels = R.R1.domains[domain].ncols;
			esint cols = lsize(map);
			for (esint k = 0; k < kernels; ++k) {
				G.rows[Roffset[domain] + k] += cols;
			}
		};

		if (map->neigh == LMAP::DIRICHLET) {
			add(map->from);
		} else {
			if (K.decomposition->ismy(map->from)) {
				add(map->from);
			}
			if (K.decomposition->ismy(map->to)) {
				add(map->to);
			}
		}
	}
	utils::sizesToOffsets(G.rows, G.rows + G.nrows + 1);
	G.resize(nrows, feti.sinfo.lambdasLocal, G.rows[G.nrows]);

	std::vector<esint> rpointer(G.nrows);
	Goffset.reserve(L.lmap.size());
	for (auto map = L.lmap.cbegin(); map != L.lmap.cend(); ++map) {
		auto add = [&] (esint domain) {
			domain -= K.decomposition->dbegin;
			esint kernels = R.R1.domains[domain].ncols;
			esint cols = lsize(map);
			for (esint k = 0; k < kernels; ++k) {
				for (esint c = 0, r = Roffset[domain] + k; c < cols; ++c) {
					G.cols[G.rows[r] + rpointer[r]++] = map->offset + c;
				}
			}
		};

		Goffset.push_back(std::make_pair(-1, -1));
		if (map->neigh == LMAP::DIRICHLET) {
			add(map->from);
			Goffset.back().first = rpointer[Roffset[map->from - K.decomposition->dbegin]];
		} else {
			if (K.decomposition->ismy(map->from)) {
				Goffset.back().first = rpointer[Roffset[map->from - K.decomposition->dbegin]];
				add(map->from);
			}
			if (K.decomposition->ismy(map->to)) {
				Goffset.back().second = rpointer[Roffset[map->to - K.decomposition->dbegin]];
				add(map->to);
			}
		}
	}
	eslog::checkpointln("FETI: SET G");
}

template<typename T>
void OrthogonalTFETISymmetric<T>::_setSparseGGt()
{
	const Matrix_FETI<Matrix_CSR, T> &K = feti.K;
	const typename FETI<T>::Regularization &R = feti.regularization;
	const typename FETI<T>::EqualityConstraints &L = feti.equalityConstraints;
	const int IDX = Indexing::CSR;

	nKernels.resize(K.decomposition->neighbors.size());
	sBuffer.resize(K.decomposition->neighbors.size());
	rBuffer.resize(K.decomposition->neighbors.size());

	LMAPSize lsize(L.lmap.cend(), feti.sinfo.lambdasLocal);
	std::vector<std::vector<std::pair<esint, esint> > > kernels(K.decomposition->neighbors.size());
	std::vector<esint> nsize(K.decomposition->neighbors.size());
	for (auto p = L.ordered.begin(); p != L.ordered.end(); ++p) {
		if (L.lmap[*p].from < K.decomposition->dbegin) { // halo
			kernels[L.lmap[*p].neigh].push_back(std::make_pair(R.R1.domains[L.lmap[*p].to - K.decomposition->dbegin].ncols, feti.sinfo.R1offset + Roffset[L.lmap[*p].to - K.decomposition->dbegin]));
			nsize[L.lmap[*p].neigh] += lsize(L.lmap.cbegin() + *p) * R.R1.domains[L.lmap[*p].to - K.decomposition->dbegin].ncols;
		}
		if (K.decomposition->dend <= L.lmap[*p].to) { // holder
			nKernels[L.lmap[*p].neigh].push_back(std::make_pair(-1, L.lmap[*p].from));
		}
	}

	if (!Communication::receiveUpperKnownSize(kernels, nKernels, K.decomposition->neighbors)) {
		eslog::error("cannot exchange kernels sizes.\n");
	}

	GGtOffset = G.nrows;
	std::vector<esint> noffset(K.decomposition->neighbors.size());
	for (auto p = L.ordered.begin(); p != L.ordered.end(); ++p) {
		if (K.decomposition->dend <= L.lmap[*p].to) { // holder
			GGtOffset += R.R1.domains[L.lmap[*p].from - K.decomposition->dbegin].ncols * nKernels[L.lmap[*p].neigh][noffset[L.lmap[*p].neigh]++].first;
		}
		if (L.lmap[*p].neigh == LMAP::LOCAL) {
			GGtOffset += R.R1.domains[L.lmap[*p].from - K.decomposition->dbegin].ncols * R.R1.domains[L.lmap[*p].to - K.decomposition->dbegin].ncols;
		}
	}
	GGtSize = GGtOffset;
	Communication::exscan(&GGtOffset, nullptr, 1, MPITools::getType(GGtOffset).mpitype, MPI_SUM);
	if (info::mpi::rank == 0) {
		GGtOffset = 0;
	}
	esint nnz = GGtOffset + GGtSize;
	Communication::broadcast(&nnz, 1, MPITools::getType(nnz).mpitype, info::mpi::size - 1);
	GGt.resize(feti.sinfo.R1totalSize, feti.sinfo.R1totalSize, nnz);

	esint GGtOffsetIterator = GGtOffset;
	GGt.rows[0] = IDX;
	std::fill(noffset.begin(), noffset.end(), 0);
	for (size_t d = 0, i = 0, r = feti.sinfo.R1offset + 1; d < K.domains.size(); ++d) {
		while (i < L.ordered.size() && L.lmap[L.ordered[i]].from - K.decomposition->dbegin < (esint)d) { ++i; }
		for (esint k = 0; k < R.R1.domains[d].ncols; ++k, ++r) {
			size_t j = i;
			GGt.cols[GGtOffsetIterator++] = feti.sinfo.R1offset + Roffset[d] + IDX; // diagonal is always filled
			while (j < L.ordered.size() && L.lmap[L.ordered[j]].from - K.decomposition->dbegin == (esint)d) {
				auto next = L.lmap.cbegin() + L.ordered[j];
				if (next->neigh == LMAP::LOCAL) {
					for (esint c = 0; c < R.R1.domains[next->to - K.decomposition->dbegin].ncols; ++c) {
						GGt.cols[GGtOffsetIterator++] = feti.sinfo.R1offset + Roffset[next->to - K.decomposition->dbegin] + c + IDX;
					}
				} else if (next->neigh != LMAP::DIRICHLET) {
					for (esint c = 0; c < nKernels[next->neigh][noffset[next->neigh]].first; ++c) {
						GGt.cols[GGtOffsetIterator++] = nKernels[next->neigh][noffset[next->neigh]].second + c + IDX;
					}
					nsize[next->neigh] += lsize(next) * nKernels[next->neigh][noffset[next->neigh]].first;
					++noffset[next->neigh];
				}
				++j;
			}
			GGt.rows[r] = GGtOffsetIterator + IDX;
		}
	}

	for (size_t n = 0; n < K.decomposition->neighbors.size(); ++n) {
		if (K.decomposition->neighbors[n] < info::mpi::rank) {
			sBuffer[n].resize(nsize[n]);
		} else {
			rBuffer[n].resize(nsize[n]);
		}
	}

	if (!Communication::allGatherInplace(GGt.rows, feti.sinfo.R1offset + 1, G.nrows)) {
		eslog::error("cannot gather GGt rows.\n");
	}
	if (!Communication::allGatherInplace(GGt.cols, GGtOffset, GGtSize)) {
		eslog::error("cannot gather GGt cols.\n");
	}

	eslog::checkpointln("FETI: GATHER GGT INDICES");

	GGt.shape = Matrix_Shape::UPPER;
	GGt.type = Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE;
	GGtSolver.commit(GGt);
	GGtSolver.symbolicFactorization();

	invGGt.resize(G.nrows, GGt.ncols);
	eslog::checkpointln("FETI: GGT SYMBOLIC FACTORIZATION");
}

template<typename T>
void OrthogonalTFETISymmetric<T>::_setDenseGGt()
{

}

template<typename T>
void OrthogonalTFETISymmetric<T>::_updateG()
{
	// G is stored with 0-based in indexing
	const Matrix_FETI<Matrix_CSR, T> &K = feti.K;
	const typename FETI<T>::Regularization &R = feti.regularization;
	const typename FETI<T>::EqualityConstraints &L = feti.equalityConstraints;

	#pragma omp parallel for
	for (size_t d = 0; d < K.domains.size(); ++d) {
		for (esint k = 0, r = Roffset[d]; k < R.R1.domains[d].ncols; ++k, ++r) {
			for (esint c = 0; c < L.domain[d].B1.nrows; ++c) {
				G.vals[G.rows[r] + c] = 0;
				for (esint i = L.domain[d].B1.rows[c]; i < L.domain[d].B1.rows[c + 1]; ++i) {
					G.vals[G.rows[r] + c] -= R.R1.domains[d].vals[R.R1.domains[d].nrows * k + L.domain[d].B1.cols[i]] * L.domain[d].B1.vals[i];
				}
			}
		}
	}
	eslog::checkpointln("FETI: UPDATE G");
}

template<typename T>
void OrthogonalTFETISymmetric<T>::_updateSparseGGt()
{
	const Matrix_FETI<Matrix_CSR, T> &K = feti.K;
	const typename FETI<T>::Regularization &R = feti.regularization;
	const typename FETI<T>::EqualityConstraints &L = feti.equalityConstraints;

	for (size_t n = 0; n < sBuffer.size(); ++n) {
		sBuffer[n].clear();
	}

	LMAPSize lsize(L.lmap.cend(), feti.sinfo.lambdasLocal);
	for (auto p = L.ordered.begin(); p != L.ordered.end(); ++p) {
		if (L.lmap[*p].from < K.decomposition->dbegin) { // halo
			esint r = Roffset[L.lmap[*p].to - K.decomposition->dbegin];
			esint size = lsize(L.lmap.cbegin() + *p);
			for (esint k = 0; k < R.R1.domains[L.lmap[*p].to - K.decomposition->dbegin].ncols; ++k) {
				esint rbegin = G.rows[r + k] + Goffset[*p].second;
				sBuffer[L.lmap[*p].neigh].insert(sBuffer[L.lmap[*p].neigh].end(), G.vals + rbegin, G.vals + rbegin + size);
			}
		}
	}

	if (!Communication::receiveUpperKnownSize(sBuffer, rBuffer, K.decomposition->neighbors)) {
		eslog::error("cannot exchange G values.\n");
	}

	esint GGtOffsetIterator = GGtOffset;
	std::vector<esint> noffset(K.decomposition->neighbors.size()), npointer(K.decomposition->neighbors.size());
	for (size_t d = 0, i = 0, r = 0; d < K.domains.size(); ++d) {
		while (i < L.ordered.size() && L.lmap[L.ordered[i]].from - K.decomposition->dbegin < (esint)d) { ++i; }
		for (esint k = 0; k < R.R1.domains[d].ncols; ++k, ++r) {
			auto mult = [] (const T *a, const T *b, esint size) {
				T sum = 0;
				for (esint i = 0; i < size; ++i) {
					sum += a[i] * b[i];
				}
				return sum;
			};

			size_t j = i;
			GGt.vals[GGtOffsetIterator++] = mult(G.vals + G.rows[r], G.vals + G.rows[r], G.rows[r + 1] - G.rows[r]);
			while (j < L.ordered.size() && L.lmap[L.ordered[j]].from - K.decomposition->dbegin == (esint)d) {
				esint fbegin = G.rows[r + k] + Goffset[L.ordered[j]].first;
				auto next = L.lmap.cbegin() + L.ordered[j];
				if (next->neigh == LMAP::LOCAL) {
					esint to = next->to - K.decomposition->dbegin, size = lsize(next);
					for (esint c = 0; c < R.R1.domains[to].ncols; ++c) {
						esint tbegin = G.rows[Roffset[to] + c] + Goffset[L.ordered[j]].second;
						GGt.vals[GGtOffsetIterator++] = mult(G.vals + fbegin, G.vals + tbegin, size);
					}
				} else if (next->neigh != LMAP::DIRICHLET) {
					esint size = lsize(next);
					for (esint c = 0; c < nKernels[next->neigh][noffset[next->neigh]].first; ++c) {
						GGt.vals[GGtOffsetIterator++] = mult(G.vals + fbegin, rBuffer[next->neigh].data() + npointer[next->neigh], size);
						npointer[next->neigh] += size;
					}
				}
				++j;
			}
		}
	}

	if (!Communication::allGatherInplace(GGt.vals, GGtOffset, GGtSize)) {
		eslog::error("cannot gather GGt vals.\n");
	}
	eslog::checkpointln("FETI: GATHER GGT VALUES");
	GGtSolver.numericalFactorization();
	eslog::checkpointln("FETI: GGT NUMERICAL FACTORIZATION");

	Matrix_Dense<T> eye;
	eye.resize(G.nrows, feti.sinfo.R1totalSize);
	math::set(eye, T{});
	for (esint r = 0; r < G.nrows; ++r) {
		eye.vals[r * feti.sinfo.R1totalSize + feti.sinfo.R1offset + r] = T{1};
	}
	GGtSolver.solve(eye, invGGt);
	eslog::checkpointln("FETI: COMPUTE GGT INVERSE");
}

template<typename T>
void OrthogonalTFETISymmetric<T>::_updateDenseGGt()
{

}

template<typename T>
void OrthogonalTFETISymmetric<T>::_print(const step::Step &step)
{
	if (info::ecf->output.print_matrices) {
		eslog::storedata(" STORE: feti/projector/{G, e, GGt, invGGt}\n");
		math::store(G, utils::filename(utils::debugDirectory(step) + "/feti/projector", "G").c_str());
		math::store(e, utils::filename(utils::debugDirectory(step) + "/feti/projector", "e").c_str());
		math::store(GGt, utils::filename(utils::debugDirectory(step) + "/feti/projector", "GGt").c_str());
		math::store(invGGt, utils::filename(utils::debugDirectory(step) + "/feti/projector", "invGGt").c_str());
	}
}

template struct OrthogonalTFETISymmetric<double>;
template struct OrthogonalTFETISymmetric<std::complex<double> >;

}

