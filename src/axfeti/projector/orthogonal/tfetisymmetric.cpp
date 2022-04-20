

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
template <typename T> static void _applyInvGGt(OrthogonalTFETISymmetric<T> *projector, const Vector_Kernel<T> &in, Vector_Dense<T> &out);
template <typename T> static void _applyGt(OrthogonalTFETISymmetric<T> *projector, const Vector_Dense<T> &in, Vector_Dual<T> &out);

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
	projector->e.resize();
	projector->Gx.resize();
	projector->iGGtGx.resize(projector->feti->sinfo.R1size);

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
	x.copyToWithoutHalo(y);
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

template <typename T> static void _applyG(OrthogonalTFETISymmetric<T> *projector, const Vector_Dual<T> &in, Vector_Kernel<T> &out)
{
	#pragma omp parallel for
	for (int t = 0; t < info::env::threads; ++t) {
		for (size_t r = Vector_Kernel<T>::distribution[t]; r < Vector_Kernel<T>::distribution[t + 1]; ++r) {
			out.vals[r + Vector_Kernel<T>::offset] = T{0};
			for (esint c = projector->G.rows[r]; c < projector->G.rows[r + 1]; ++c) {
				out.vals[r + Vector_Kernel<T>::offset] += projector->G.vals[c] * in.vals[projector->G.cols[c]];
			}
		}
	}
	out.synchronize();
}

template <typename T> static void _applyInvGGt(OrthogonalTFETISymmetric<T> *projector, const Vector_Kernel<T> &in, Vector_Dense<T> &out)
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
template <typename T> static void _applyGt(OrthogonalTFETISymmetric<T> *projector, const Vector_Dense<T> &in, Vector_Dual<T> &out)
{
	for (esint r = 0; r < projector->G.nrows; ++r) {
		for (esint c = projector->G.rows[r]; c < projector->G.rows[r + 1]; ++c) {
			out.vals[projector->G.cols[c]] += projector->G.vals[c] * in.vals[r];
		}
	}
	out.synchronize();
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
	projector->Goffset.reserve(L->lmap.size());
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

		projector->Goffset.push_back(std::make_pair(-1, -1));
		if (map->neigh == LMAP::DIRICHLET) {
			add(map->from);
			projector->Goffset.back().first = rpointer[projector->Roffset[map->from - K->decomposition->dbegin]];
		} else {
			if (K->decomposition->ismy(map->from)) {
				projector->Goffset.back().first = rpointer[projector->Roffset[map->from - K->decomposition->dbegin]];
				add(map->from);
			}
			if (K->decomposition->ismy(map->to)) {
				projector->Goffset.back().second = rpointer[projector->Roffset[map->to - K->decomposition->dbegin]];
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
	for (auto p = L->ordered.begin(); p != L->ordered.end(); ++p) {
		if (L->lmap[*p].from < K->decomposition->dbegin) { // halo
			kernels[L->lmap[*p].neigh].push_back(std::make_pair(R->R1.domains[L->lmap[*p].to - K->decomposition->dbegin].ncols, projector->feti->sinfo.R1offset + projector->Roffset[L->lmap[*p].to - K->decomposition->dbegin]));
			nsize[L->lmap[*p].neigh] += lsize(L->lmap.cbegin() + *p) * R->R1.domains[L->lmap[*p].to - K->decomposition->dbegin].ncols;
		}
		if (K->decomposition->dend <= L->lmap[*p].to) { // holder
			projector->nKernels[L->lmap[*p].neigh].push_back(std::make_pair(-1, L->lmap[*p].from));
		}
	}

	if (!Communication::receiveUpperKnownSize(kernels, projector->nKernels, K->decomposition->neighbors)) {
		eslog::error("cannot exchange kernels sizes.\n");
	}

	projector->GGtOffset = projector->G.nrows;
	std::vector<esint> noffset(K->decomposition->neighbors.size());
	for (auto p = L->ordered.begin(); p != L->ordered.end(); ++p) {
		if (K->decomposition->dend <= L->lmap[*p].to) { // holder
			projector->GGtOffset += R->R1.domains[L->lmap[*p].from - K->decomposition->dbegin].ncols * projector->nKernels[L->lmap[*p].neigh][noffset[L->lmap[*p].neigh]++].first;
		}
		if (L->lmap[*p].neigh == LMAP::LOCAL) {
			projector->GGtOffset += R->R1.domains[L->lmap[*p].from - K->decomposition->dbegin].ncols * R->R1.domains[L->lmap[*p].to - K->decomposition->dbegin].ncols;
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
	std::fill(noffset.begin(), noffset.end(), 0);
	for (size_t d = 0, i = 0, r = projector->feti->sinfo.R1offset + 1; d < K->domains.size(); ++d) {
		while (i < L->ordered.size() && L->lmap[L->ordered[i]].from - K->decomposition->dbegin < (esint)d) { ++i; }
		for (esint k = 0; k < R->R1.domains[d].ncols; ++k, ++r) {
			size_t j = i;
			projector->GGt.cols[GGtOffset++] = projector->feti->sinfo.R1offset + projector->Roffset[d] + IDX; // diagonal is always filled
			while (j < L->ordered.size() && L->lmap[L->ordered[j]].from - K->decomposition->dbegin == (esint)d) {
				auto next = L->lmap.cbegin() + L->ordered[j];
				if (next->neigh == LMAP::LOCAL) {
					for (esint c = 0; c < R->R1.domains[next->to - K->decomposition->dbegin].ncols; ++c) {
						projector->GGt.cols[GGtOffset++] = projector->feti->sinfo.R1offset + projector->Roffset[next->to - K->decomposition->dbegin] + c + IDX;
					}
				} else if (next->neigh != LMAP::DIRICHLET) {
					for (esint c = 0; c < projector->nKernels[next->neigh][noffset[next->neigh]].first; ++c) {
						projector->GGt.cols[GGtOffset++] = projector->nKernels[next->neigh][noffset[next->neigh]].second + c + IDX;
					}
					nsize[next->neigh] += lsize(next) * projector->nKernels[next->neigh][noffset[next->neigh]].first;
					++noffset[next->neigh];
				}
				++j;
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
	for (auto p = L->ordered.begin(); p != L->ordered.end(); ++p) {
		if (L->lmap[*p].from < K->decomposition->dbegin) { // halo
			esint r = projector->Roffset[L->lmap[*p].to - K->decomposition->dbegin];
			esint size = lsize(L->lmap.cbegin() + *p);
			for (esint k = 0; k < R->R1.domains[L->lmap[*p].to - K->decomposition->dbegin].ncols; ++k) {
				esint rbegin = projector->G.rows[r + k] + projector->Goffset[*p].second;
				projector->sBuffer[L->lmap[*p].neigh].insert(projector->sBuffer[L->lmap[*p].neigh].end(), projector->G.vals + rbegin, projector->G.vals + rbegin + size);
			}
		}
	}

	if (!Communication::receiveUpperKnownSize(projector->sBuffer, projector->rBuffer, K->decomposition->neighbors)) {
		eslog::error("cannot exchange G values.\n");
	}

	esint GGtOffset = projector->GGtOffset;
	std::vector<esint> noffset(K->decomposition->neighbors.size()), npointer(K->decomposition->neighbors.size());
	for (size_t d = 0, i = 0, r = 0; d < K->domains.size(); ++d) {
		while (i < L->ordered.size() && L->lmap[L->ordered[i]].from - K->decomposition->dbegin < (esint)d) { ++i; }
		for (esint k = 0; k < R->R1.domains[d].ncols; ++k, ++r) {
			auto mult = [] (const T *a, const T *b, esint size) {
				T sum = 0;
				for (esint i = 0; i < size; ++i) {
					sum += a[i] * b[i];
				}
				return sum;
			};

			size_t j = i;
			projector->GGt.vals[GGtOffset++] = mult(projector->G.vals + projector->G.rows[r], projector->G.vals + projector->G.rows[r], projector->G.rows[r + 1] - projector->G.rows[r]);
			while (j < L->ordered.size() && L->lmap[L->ordered[j]].from - K->decomposition->dbegin == (esint)d) {
				esint fbegin = projector->G.rows[r + k] + projector->Goffset[L->ordered[j]].first;
				auto next = L->lmap.cbegin() + L->ordered[j];
				if (next->neigh == LMAP::LOCAL) {
					esint to = next->to - K->decomposition->dbegin, size = lsize(next);
					for (esint c = 0; c < R->R1.domains[to].ncols; ++c) {
						esint tbegin = projector->G.rows[projector->Roffset[to] + c] + projector->Goffset[L->ordered[j]].second;
						projector->GGt.vals[GGtOffset++] = mult(projector->G.vals + fbegin, projector->G.vals + tbegin, size);
					}
				} else if (next->neigh != LMAP::DIRICHLET) {
					esint size = lsize(next);
					for (esint c = 0; c < projector->nKernels[next->neigh][noffset[next->neigh]].first; ++c) {
						projector->GGt.vals[GGtOffset++] = mult(projector->G.vals + fbegin, projector->rBuffer[next->neigh].data() + npointer[next->neigh], size);
						npointer[next->neigh] += size;
					}
				}
				++j;
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

