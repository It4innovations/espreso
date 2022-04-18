

#include "tfetisymmetric.h"
#include "basis/containers/serializededata.h"
#include "basis/utilities/sysutils.h"
#include "basis/utilities/utils.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/envinfo.h"
#include "esinfo/eslog.hpp"

#include <vector>

namespace espreso {

template <typename T> static void _setG(OrthogonalTFETISymmetric<T> *projector);
template <typename T> static void _setSparseGGt(OrthogonalTFETISymmetric<T> *projector);
template <typename T> static void _setDenseGGt(OrthogonalTFETISymmetric<T> *projector);
template <typename T> static void _updateG(OrthogonalTFETISymmetric<T> *projector);
template <typename T> static void _updateSparseGGt(OrthogonalTFETISymmetric<T> *projector);
template <typename T> static void _updateDenseGGt(OrthogonalTFETISymmetric<T> *projector);

template <typename T> static void _applyG(OrthogonalTFETISymmetric<T> *projector, const Vector_Dual<T> &in, Vector_Kernel<T> &out);
template <typename T> static void _applyInvGGt(OrthogonalTFETISymmetric<T> *projector, const Vector_Kernel<T> &in, Vector_Kernel<T> &out);
template <typename T> static void _applyGt(OrthogonalTFETISymmetric<T> *projector, const Vector_Kernel<T> &in, Vector_Dual<T> &out);

template <typename T> static void _print(OrthogonalTFETISymmetric<T> *projector);

template <typename T>
static void _info(OrthogonalTFETISymmetric<T> *projector)
{
	esint nnz = 2 * (projector->GGt.nnz - projector->GGt.nrows) + projector->GGt.nrows;

	eslog::info(" = ORTHOGONAL PROJECTOR PROPERTIES                                                           = \n");
	eslog::info(" =   GGt SIZE                                                                      %9d = \n", projector->GGt.nrows);
	eslog::info(" =   GGt FILL-IN [\%]                                                                %8.4f = \n", 100.0 * nnz / (projector->GGt.nrows * projector->GGt.nrows));
	if (projector->feti->configuration.exhaustive_info) {
		// PPt = eye
	}
	eslog::info(" = ----------------------------------------------------------------------------------------- = \n");

}

template <typename T>
static void _set(OrthogonalTFETISymmetric<T> *projector)
{
	const typename AX_FETI<T>::EqualityConstraints *L = projector->feti->equalityConstraints;

	projector->e.resize();
	projector->Gx.resize();
	projector->iGGtGx.resize();

	projector->mapHalo = 0;
	while (L->lmap[projector->mapHalo].from < projector->feti->sinfo.R1offset) { ++projector->mapHalo; }
	projector->mapNN = projector->mapHalo;
	while (projector->mapNN < L->lmap.size() && L->lmap[projector->mapNN].neigh != LMAP::LOCAL && L->lmap[projector->mapNN].neigh != LMAP::DIRICHLET) { ++projector->mapNN; }
	projector->mapLocal = projector->mapNN;
	while (projector->mapLocal < L->lmap.size() && L->lmap[projector->mapLocal].neigh != LMAP::DIRICHLET) { ++projector->mapLocal; }

	_setG(projector);
	if (projector->feti->equalityConstraints->global) {
		_setDenseGGt(projector);
	} else {
		_setSparseGGt(projector);
	}
}

template <typename T>
static void _free(OrthogonalTFETISymmetric<T> *projector)
{
	math::free(projector->GGt);
}

template <typename T>
static void _update(OrthogonalTFETISymmetric<T> *projector)
{
	const typename AX_FETI<T>::Regularization *R = projector->feti->regularization;

	#pragma omp parallel for
	for (size_t d = 0; d < R->R1.domains.size(); ++d) {
		Vector_Dense<T> _e;
		_e.size = R->R1.domains[d].ncols;
		_e.vals = projector->e.vals + projector->Roffset[d] + Vector_Kernel<T>::offset;
		Matrix_Dense<T> _Rt;
		_Rt.nrows = R->R1.domains[d].ncols;
		_Rt.ncols = R->R1.domains[d].nrows;
		_Rt.vals = R->R1.domains[d].vals;
		math::apply(_e, T{1}, _Rt, T{0}, projector->feti->f->domains[d]);
	}
	projector->e.synchronize();
	eslog::checkpointln("FETI: COMPUTE DUAL RHS [e]");

	_updateG(projector);
	if (projector->feti->equalityConstraints->global) {
		_updateDenseGGt(projector);
	} else {
		_updateSparseGGt(projector);
	}

	_print(projector);
}

template <typename T>
static void _apply(OrthogonalTFETISymmetric<T> *projector, const Vector_Dual<T> &x, Vector_Dual<T> &y)
{
	x.copyTo(y);
	_applyG(projector, x, projector->Gx);
	_applyInvGGt(projector, projector->Gx, projector->iGGtGx);
	_applyGt(projector, projector->iGGtGx, y);
}

template <typename T>
static void _applyGtInvGGt(OrthogonalTFETISymmetric<T> *projector, const Vector_Kernel<T> &x, Vector_Dual<T> &y)
{
	math::set(y, T{0});
	_applyInvGGt(projector, x, projector->iGGtGx);
	_applyGt(projector, projector->iGGtGx, y);
}

template <typename T>
static void _applyInvGGtG(OrthogonalTFETISymmetric<T> *projector, const Vector_Dual<T> &x, Vector_Kernel<T> &y)
{
	_applyG(projector, x, projector->Gx);
	_applyInvGGt(projector, projector->Gx, y);
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

struct LMAPOrder {
	LMAPOrder(const std::vector<LMAP> &map, size_t halo, size_t nn): map(map), halo(halo), nn(nn) { }

	std::vector<LMAP>::const_iterator next(esint domain) {
		std::vector<LMAP>::const_iterator hpair = map.cend(), nnpair = map.cend();
		if (
				halo < map.size() &&
				(map.cbegin() + halo)->neigh != LMAP::LOCAL &&
				(map.cbegin() + halo)->neigh != LMAP::DIRICHLET &&
				(map.cbegin() + halo)->from == domain) {
			hpair = map.cbegin() + halo;
		}
		if (
				nn < map.size() &&
				(map.cbegin() + nn)->neigh == LMAP::LOCAL &&
				(map.cbegin() + nn)->from == domain) {
			nnpair = map.cbegin() + nn;
		}
		if (hpair == map.cend() || hpair->from != domain) {
			if (nn < map.size() && (map.cbegin() + nn)->from == domain) {
				++nn;
				return nnpair;
			}
			return map.cend();
		}
		if (nnpair == map.cend() || nnpair->from != domain) {
			if (halo < map.size() && (map.cbegin() + halo)->from == domain) {
				++halo;
				return hpair;
			}
			return map.cend();
		}
		if (hpair->to < nnpair->to) {
			return map.cbegin() + halo++;
		} else {
			return map.cbegin() + nn++;
		}
	}

	const std::vector<LMAP> &map;
	size_t halo, nn;
};

template <typename T>
static void _setG(OrthogonalTFETISymmetric<T> *projector)
{
	// G is stored with 0-based in indexing
	const Matrix_FETI<Matrix_CSR, T> *K = projector->feti->K;
	const typename AX_FETI<T>::Regularization *R = projector->feti->regularization;
	const typename AX_FETI<T>::EqualityConstraints *L = projector->feti->equalityConstraints;

	esint nrows = 0;
	projector->Roffset.resize(K->domains.size());
	for (size_t d = 0; d < K->domains.size(); ++d) {
		projector->Roffset[d] = nrows;
		nrows += R->R1.domains[d].ncols;
	}
	projector->G.resize(nrows, projector->feti->sinfo.lambdasLocal, 0);
	math::set(projector->G.nrows + 1, projector->G.rows, 1, 0);

	LMAPSize lsize(L->lmap.cend(), projector->feti->sinfo.lambdasLocal);
	for (auto map = L->lmap.cbegin(); map != L->lmap.cend(); ++map) {
		auto add = [&] (esint domain) {
			domain -= K->decomposition->dbegin;
			esint kernels = R->R1.domains[domain].ncols;
			esint cols = lsize(map);
			for (esint k = 0; k < kernels; ++k) {
				projector->G.rows[projector->Roffset[domain] + k] += cols;
			}
		};

		if (map->neigh == LMAP::DIRICHLET) {
			add(map->from);
		} else {
			if (K->decomposition->ismy(map->from)) {
				add(map->from);
			}
			if (K->decomposition->ismy(map->to)) {
				add(map->to);
			}
		}
	}
	utils::sizesToOffsets(projector->G.rows, projector->G.rows + projector->G.nrows + 1);
	projector->G.resize(nrows, projector->feti->sinfo.lambdasLocal, projector->G.rows[projector->G.nrows]);

	std::vector<esint> rpointer(projector->G.nrows);
	for (auto map = L->lmap.cbegin(); map != L->lmap.cend(); ++map) {
		auto add = [&] (esint domain) {
			domain -= K->decomposition->dbegin;
			esint kernels = R->R1.domains[domain].ncols;
			esint cols = lsize(map);
			for (esint k = 0; k < kernels; ++k) {
				for (esint c = 0, r = projector->Roffset[domain] + k; c < cols; ++c) {
					projector->G.cols[projector->G.rows[r] + rpointer[r]++] = map->offset + c;
				}
			}
		};

		if (map->neigh == LMAP::DIRICHLET) {
			add(map->from);
		} else {
			if (K->decomposition->ismy(map->from)) {
				add(map->from);
			}
			if (K->decomposition->ismy(map->to)) {
				add(map->to);
			}
		}
	}
	eslog::checkpointln("FETI: SET G");
}

template <typename T>
static void _setSparseGGt(OrthogonalTFETISymmetric<T> *projector)
{
	const Matrix_FETI<Matrix_CSR, T> *K = projector->feti->K;
	const typename AX_FETI<T>::Regularization *R = projector->feti->regularization;
	const typename AX_FETI<T>::EqualityConstraints *L = projector->feti->equalityConstraints;
	const int IDX = Indexing::CSR;

	projector->nKernels.resize(K->decomposition->neighbors.size());
	projector->sBuffer.resize(K->decomposition->neighbors.size());
	projector->rBuffer.resize(K->decomposition->neighbors.size());

	LMAPSize lsize(L->lmap.cend(), projector->feti->sinfo.lambdasLocal);
	std::vector<std::vector<std::pair<esint, esint> > > kernels(K->decomposition->neighbors.size());
	std::vector<esint> nsize(K->decomposition->neighbors.size());
	for (auto map = L->lmap.cbegin(); map != L->lmap.cbegin() + projector->mapHalo; ++map) {
		kernels[map->neigh].push_back(std::make_pair(R->R1.domains[map->to - K->decomposition->dbegin].ncols, projector->feti->sinfo.R1offset + projector->Roffset[map->to - K->decomposition->dbegin]));
		nsize[map->neigh] += lsize(map) * R->R1.domains[map->to - K->decomposition->dbegin].ncols;
	}
	for (auto map = L->lmap.cbegin() + projector->mapHalo; map != L->lmap.cbegin() + projector->mapNN; ++map) {
		projector->nKernels[map->neigh].push_back(std::make_pair(-1, -1)); // dummy
	}

	if (!Communication::receiveUpperKnownSize(kernels, projector->nKernels, K->decomposition->neighbors)) {
		eslog::error("cannot exchange kernels sizes.\n");
	}

	projector->GGtOffset = projector->G.nrows;
	std::vector<esint> noffset(K->decomposition->neighbors.size());
	for (auto map = L->lmap.cbegin() + projector->mapHalo; map != L->lmap.cbegin() + projector->mapLocal; ++map) {
		if (map->neigh != LMAP::LOCAL) {
			projector->GGtOffset += R->R1.domains[map->from - K->decomposition->dbegin].ncols * projector->nKernels[map->neigh][noffset[map->neigh]++].first;
		} else {
			projector->GGtOffset += R->R1.domains[map->from - K->decomposition->dbegin].ncols * R->R1.domains[map->to - K->decomposition->dbegin].ncols;
		}
	}

	projector->GGtSize = projector->GGtOffset;
	Communication::exscan(&projector->GGtOffset, nullptr, 1, MPITools::getType(projector->GGtOffset).mpitype, MPI_SUM);
	if (info::mpi::rank == 0) {
		projector->GGtOffset = 0;
	}
	esint nnz = projector->GGtOffset + projector->GGtSize;
	Communication::broadcast(&nnz, 1, MPITools::getType(nnz).mpitype, info::mpi::size - 1);
	projector->GGt.resize(projector->feti->sinfo.R1totalSize, projector->feti->sinfo.R1totalSize, nnz);

	esint GGtOffset = projector->GGtOffset;
	projector->GGt.rows[0] = IDX;
	LMAPOrder order(L->lmap, projector->mapHalo, projector->mapNN);
	std::fill(noffset.begin(), noffset.end(), 0);
	for (size_t d = 0, r = projector->feti->sinfo.R1offset + 1; d < K->domains.size(); ++d) {
		for (esint k = 0; k < R->R1.domains[d].ncols; ++k, ++r) {
			projector->GGt.cols[GGtOffset++] = projector->feti->sinfo.R1offset + projector->Roffset[d] + IDX; // diagonal is always filled
			std::vector<LMAP>::const_iterator next;
			while ((next = order.next(K->decomposition->dbegin + d)) != L->lmap.cend()) {
				if (next->neigh == LMAP::LOCAL) {
					for (esint c = 0; c < R->R1.domains[next->to - K->decomposition->dbegin].ncols; ++c) {
						projector->GGt.cols[GGtOffset++] = projector->feti->sinfo.R1offset + projector->Roffset[next->to - K->decomposition->dbegin] + c + IDX;
					}
				} else {
					for (esint c = 0; c < projector->nKernels[next->neigh][noffset[next->neigh]].first; ++c) {
						projector->GGt.cols[GGtOffset++] = projector->nKernels[next->neigh][noffset[next->neigh]].second + c + IDX;
					}
					nsize[next->neigh] += lsize(next) * projector->nKernels[next->neigh][noffset[next->neigh]].first;
					++noffset[next->neigh];
				}
			}
			projector->GGt.rows[r] = GGtOffset + IDX;
		}
	}

	for (size_t n = 0; n < K->decomposition->neighbors.size(); ++n) {
		if (K->decomposition->neighbors[n] < info::mpi::rank) {
			projector->sBuffer[n].resize(nsize[n]);
		} else {
			projector->rBuffer[n].resize(nsize[n]);
		}
	}

	if (!Communication::allGatherInplace(projector->GGt.rows, projector->feti->sinfo.R1offset + 1, projector->G.nrows)) {
		eslog::error("cannot gather GGt rows.\n");
	}
	if (!Communication::allGatherInplace(projector->GGt.cols, projector->GGtOffset, projector->GGtSize)) {
		eslog::error("cannot gather GGt cols.\n");
	}

	eslog::checkpointln("FETI: GATHER GGT INDICES");

	projector->GGt.shape = Matrix_Shape::UPPER;
	projector->GGt.type = Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE;
	math::commit(projector->GGt);
	math::symbolicFactorization(projector->GGt);

	projector->invGGt.resize(projector->G.nrows, projector->GGt.ncols);
	eslog::checkpointln("FETI: GGT SYMBOLIC FACTORIZATION");
}

template <typename T>
static void _setDenseGGt(OrthogonalTFETISymmetric<T> *projector)
{

}

template <typename T>
static void _updateG(OrthogonalTFETISymmetric<T> *projector)
{
	// G is stored with 0-based in indexing
	const Matrix_FETI<Matrix_CSR, T> *K = projector->feti->K;
	const typename AX_FETI<T>::Regularization *R = projector->feti->regularization;
	const typename AX_FETI<T>::EqualityConstraints *L = projector->feti->equalityConstraints;

	#pragma omp parallel for
	for (size_t d = 0; d < K->domains.size(); ++d) {
		for (esint k = 0, r = projector->Roffset[d]; k < R->R1.domains[d].ncols; ++k, ++r) {
			for (esint c = 0; c < L->domain[d].B1.nrows; ++c) {
				projector->G.vals[projector->G.rows[r] + c] = 0;
				for (esint i = L->domain[d].B1.rows[c]; i < L->domain[d].B1.rows[c + 1]; ++i) {
					projector->G.vals[projector->G.rows[r] + c] -= R->R1.domains[d].vals[R->R1.domains[d].nrows * k + L->domain[d].B1.cols[i]] * L->domain[d].B1.vals[i];
				}
			}
		}
	}
	eslog::checkpointln("FETI: UPDATE G");
}

template <typename T>
static void _updateSparseGGt(OrthogonalTFETISymmetric<T> *projector)
{
	const Matrix_FETI<Matrix_CSR, T> *K = projector->feti->K;
	const typename AX_FETI<T>::Regularization *R = projector->feti->regularization;
	const typename AX_FETI<T>::EqualityConstraints *L = projector->feti->equalityConstraints;

	for (size_t n = 0; n < projector->sBuffer.size(); ++n) {
		projector->sBuffer[n].clear();
	}

	LMAPSize lsize(L->lmap.cend(), projector->feti->sinfo.lambdasLocal);
	std::vector<esint> rpointer(projector->G.rows, projector->G.rows + projector->G.nrows);
	for (auto map = L->lmap.cbegin(); map != L->lmap.cbegin() + projector->mapHalo; ++map) {
		esint r = projector->Roffset[map->to - K->decomposition->dbegin];
		while (projector->G.cols[rpointer[r]] < map->offset) { ++rpointer[r]; }
		esint size = lsize(map);
		for (esint k = 1; k < R->R1.domains[map->to - K->decomposition->dbegin].ncols; ++k) {
			rpointer[r + k] = rpointer[r];
		}
		for (esint k = 0; k < R->R1.domains[map->to - K->decomposition->dbegin].ncols; ++k) {
			projector->sBuffer[map->neigh].insert(projector->sBuffer[map->neigh].end(), projector->G.vals + rpointer[r + k], projector->G.vals + rpointer[r + k] + size);
			rpointer[r + k] += size;
		}
	}

	if (!Communication::receiveUpperKnownSize(projector->sBuffer, projector->rBuffer, K->decomposition->neighbors)) {
		eslog::error("cannot exchange G values.\n");
	}

	esint GGtOffset = projector->GGtOffset;
	LMAPOrder order(L->lmap, projector->mapHalo, projector->mapNN);
	std::vector<esint> noffset(K->decomposition->neighbors.size()), npointer(K->decomposition->neighbors.size());;
	for (size_t d = 0, r = 0; d < K->domains.size(); ++d) {
		for (esint k = 0; k < R->R1.domains[d].ncols; ++k, ++r) {
			auto mult = [] (const T *a, const T *b, esint size) {
				T sum = 0;
				for (esint i = 0; i < size; ++i) {
					sum += a[i] * b[i];
				}
				return sum;
			};

			projector->GGt.vals[GGtOffset++] = mult(projector->G.vals + projector->G.rows[r], projector->G.vals + projector->G.rows[r], projector->G.rows[r + 1] - projector->G.rows[r]);
			std::vector<LMAP>::const_iterator next;
			while ((next = order.next(K->decomposition->dbegin + d)) != L->lmap.cend()) {
				if (next->neigh == LMAP::LOCAL) {
					esint to = next->to - K->decomposition->dbegin, size = lsize(next);
					for (esint c = 0; c < R->R1.domains[to].ncols; ++c) {
						projector->GGt.vals[GGtOffset++] = mult(projector->G.vals + rpointer[r + c], projector->G.vals + rpointer[to + c], size);
						rpointer[r + c] += size;
						rpointer[to + c] += size;
					}
				} else {
					esint size = lsize(next);
					for (esint c = 0; c < projector->nKernels[next->neigh][noffset[next->neigh]].first; ++c) {
						projector->GGt.vals[GGtOffset++] = mult(projector->G.vals + rpointer[r + c], projector->rBuffer[next->neigh].data() + npointer[next->neigh], size);
						npointer[next->neigh] += size;
					}
				}
			}
		}
	}

	if (!Communication::allGatherInplace(projector->GGt.vals, projector->GGtOffset, projector->GGtSize)) {
		eslog::error("cannot gather GGt vals.\n");
	}
	eslog::checkpointln("FETI: GATHER GGT VALUES");
	math::numericalFactorization(projector->GGt);
	eslog::checkpointln("FETI: GGT NUMERICAL FACTORIZATION");

	Matrix_Dense<T> eye;
	eye.resize(projector->G.nrows, projector->feti->sinfo.R1totalSize);
	math::set(eye, T{});
	for (esint r = 0; r < projector->G.nrows; ++r) {
		eye.vals[r * projector->feti->sinfo.R1totalSize + projector->feti->sinfo.R1offset + r] = T{1};
	}
	math::solve(projector->GGt, eye, projector->invGGt);
	eslog::checkpointln("FETI: COMPUTE GGT INVERSE");
}

template <typename T>
static void _updateDenseGGt(OrthogonalTFETISymmetric<T> *projector)
{

}

template <typename T> static void _applyG(OrthogonalTFETISymmetric<T> *projector, const Vector_Dual<T> &in, Vector_Kernel<T> &out)
{
	#pragma omp parallel for
	for (int t = 0; t < info::env::threads; ++t) {
		for (size_t r = Vector_Kernel<T>::distribution[t]; r < Vector_Kernel<T>::distribution[t + 1]; ++r) {
			out.vals[r] = T{0};
			for (esint c = projector->G.rows[r]; c < projector->G.rows[r + 1]; ++c) {
				out.vals[r] += projector->G.vals[c] * in.vals[projector->G.cols[c]];
			}
		}
	}
	out.synchronize();
}

template <typename T> static void _applyInvGGt(OrthogonalTFETISymmetric<T> *projector, const Vector_Kernel<T> &in, Vector_Kernel<T> &out)
{
	#pragma omp parallel for
	for (int t = 0; t < info::env::threads; ++t) {
		Matrix_Dense<T> a;
		Vector_Dense<T> y;
		a.ncols = projector->invGGt.ncols;
		a.nrows = y.size = Vector_Kernel<T>::distribution[t + 1] - Vector_Kernel<T>::distribution[t];

		a.vals = projector->invGGt.vals + projector->invGGt.ncols * Vector_Kernel<T>::distribution[t];
		y.vals = out.vals + Vector_Kernel<T>::distribution[t];

		math::apply(y, T{-1}, a, T{0}, in);
	}
}

// TODO: threaded implementation: utilize properties of B of gluing two domains
template <typename T> static void _applyGt(OrthogonalTFETISymmetric<T> *projector, const Vector_Kernel<T> &in, Vector_Dual<T> &out)
{
	for (esint r = 0; r < projector->G.nrows; ++r) {
		for (esint c = projector->G.rows[r]; c < projector->G.rows[r + 1]; ++c) {
			out.vals[projector->G.cols[c]] += projector->G.vals[c] * in.vals[r];
		}
	}
}

template <typename T>
static void _print(OrthogonalTFETISymmetric<T> *projector)
{
	if (info::ecf->output.print_matrices) {
		eslog::storedata(" STORE: feti/projector/{G, e, GGt, invGGt}\n");
		math::store(projector->G, utils::filename(utils::debugDirectory(*projector->feti->step) + "/feti/projector", "G").c_str());
		math::store(projector->e, utils::filename(utils::debugDirectory(*projector->feti->step) + "/feti/projector", "e").c_str());
		math::store(projector->GGt, utils::filename(utils::debugDirectory(*projector->feti->step) + "/feti/projector", "GGt").c_str());
		math::store(projector->invGGt, utils::filename(utils::debugDirectory(*projector->feti->step) + "/feti/projector", "invGGt").c_str());
	}
}

template <> OrthogonalTFETISymmetric<double>::OrthogonalTFETISymmetric(AX_FETI<double> *feti): Projector(feti) { _set<double>(this); }
template <> OrthogonalTFETISymmetric<std::complex<double> >::OrthogonalTFETISymmetric(AX_FETI<std::complex<double> > *feti): Projector(feti) { _set<std::complex<double> >(this); }

template <> OrthogonalTFETISymmetric<double>::~OrthogonalTFETISymmetric() { _free<double>(this); }
template <> OrthogonalTFETISymmetric<std::complex<double> >::~OrthogonalTFETISymmetric() { _free<std::complex<double> >(this); }

template <> void OrthogonalTFETISymmetric<double>::info() { _info<double>(this); }
template <> void OrthogonalTFETISymmetric<std::complex<double> >::info() { _info<std::complex<double> >(this); }

template <> void OrthogonalTFETISymmetric<double>::update() { _update<double>(this); }
template <> void OrthogonalTFETISymmetric<std::complex<double> >::update() { _update<std::complex<double> >(this); }

template <> void OrthogonalTFETISymmetric<double>::apply(const Vector_Dual<double> &x, Vector_Dual<double> &y) { _apply(this, x, y); }
template <> void OrthogonalTFETISymmetric<std::complex<double> >::apply(const Vector_Dual<std::complex<double> > &x, Vector_Dual<std::complex<double> > &y) { _apply(this, x, y); }

template <> void OrthogonalTFETISymmetric<double>::applyGtInvGGt(const Vector_Kernel<double> &x, Vector_Dual<double> &y) { _applyGtInvGGt(this, x, y); }
template <> void OrthogonalTFETISymmetric<std::complex<double> >::applyGtInvGGt(const Vector_Kernel<std::complex<double> > &x, Vector_Dual<std::complex<double> > &y) { _applyGtInvGGt(this, x, y); }

template <> void OrthogonalTFETISymmetric<double>::applyInvGGtG(const Vector_Dual<double> &x, Vector_Kernel<double> &y) { _applyInvGGtG(this, x, y); }
template <> void OrthogonalTFETISymmetric<std::complex<double> >::applyInvGGtG(const Vector_Dual<std::complex<double> > &x, Vector_Kernel<std::complex<double> > &y) { _applyInvGGtG(this, x, y); }

}

