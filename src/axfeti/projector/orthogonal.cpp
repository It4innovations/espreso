
#include "orthogonal.h"

#include "basis/containers/serializededata.h"
#include "basis/utilities/sysutils.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.hpp"

#include <vector>

#include "basis/utilities/print.h"

namespace espreso {

template <typename T>
static void _info(Orthogonal<T> *projector)
{
	eslog::info(" = G1 SIZE                                                                          %8d = \n", projector->feti->sinfo.R1size);
	eslog::info(" = LAMBDAS                                                                          %8d = \n", projector->feti->equalityConstraints->B1.domains.front().nrows);
}

/*
 * G = Rt * Bt
 * R -> stored in row major (Dense) but transposed (hence column major)
 * B -> stored in row major (IJV) but transposed
 *
 * G -> (number of kernels) x (number of lambdas)
 *
 *
 * 2 domains each 3 DOFs (R, B are transposed)
 * | r1 r2 r3    0  0  0 |  | l1  l2  0   0 |      | x1 x2 y1 x3 |
 * |  0  0  0   r4 r5 r6 |  |  0   0  0  l4 |      | x4 x5 y1 x6 |
 *                          |  0   0  0   0 |
 *                          |               |
 *                          |  0 -l3  0   0 |
 *                          |  0   0  0 -l5 |
 *                          |  0   0  0   0 |
 *
 * x1 = r1 * l1    x2 =  r1 * l2    x3 =  r2 * l4
 * x4 =       0    x5 = -r4 * l3    x6 = -r5 * l5
 * y1 = always zero, since lambda glues different domains (sparse format)
 *
 * number of NONZEROS for each row = LAMBDAS in a particular domain
 *
 * G
 *     L1 L2 L3 L4 L5       R1 R2 R3
 * R1   x  x     x  x   L1   x
 * R2      x  x  x      L2   x  x  x
 * R2      x  x     x   L3      x  x
 *                      L4   x  x
 *                      L5   x     x
 *
 * R1 -> R2: L2 L4
 * R1 -> R3: L2 L5
 *
 * LAMBDAS size is irrelevant here
 *
 * rows for each domain are independent on other domains' rows
 *
 *
 */

template <typename T>
static void _setG(Orthogonal<T> *projector)
{
	const Matrix_FETI<Matrix_CSR, T> *K = projector->feti->K;
	const typename AX_FETI<T>::Regularization *R = projector->feti->regularization;
	const typename AX_FETI<T>::EqualityConstraints *L = projector->feti->equalityConstraints;
	const int IDX = _Matrix_CSR_Pattern::Indexing;

	// sIndices
	// index of the first R, number of R, number of column indices, columns indices
	//
	// rIndices
	// offset to G, number of values to mem-copy

	projector->sIndices.resize(K->decomposition->neighbors.size());
	projector->rIndices.resize(K->decomposition->neighbors.size());
	projector->sBuffer.resize(K->decomposition->neighbors.size());
	projector->rBuffer.resize(K->decomposition->neighbors.size());

	// go through all lambdas and collect indices that will be communicated
	esint Roffset = projector->feti->sinfo.R1offset;
	std::vector<size_t> sOffset(K->decomposition->neighbors.size());
	for (size_t d = 0; d < K->domains.size(); ++d) {
		for (size_t n = 0; n < projector->sIndices.size(); ++n) { // here we assume that each domain send something to all
			sOffset[n] = projector->sIndices[n].size();
			projector->sIndices[n].push_back(Roffset);
			projector->sIndices[n].push_back(R->R1.domains[d].ncols);
			projector->sIndices[n].push_back(0);
		}
		for (size_t lambda = 0; lambda < L->D2C[d].size(); ++lambda) {
			auto neighs = L->L2MPI->cbegin() + L->D2C[d][lambda];
			for (auto n = neighs->begin(); n != neighs->end(); ++n) {
				++projector->sIndices[*n][sOffset[*n] + 2];
				projector->sIndices[*n].push_back(L->C2G[L->D2C[d][lambda]] + IDX);
			}
		}
		for (size_t n = 0; n < projector->sIndices.size(); ++n) {
			if (projector->sIndices[n][sOffset[n] + 2] == 0) { // these are no columns to send
				projector->sIndices[n].resize(sOffset[n]);
			}
		}
		Roffset += R->R1.domains[d].ncols;
	}

	if (!Communication::exchangeUnknownSize(projector->sIndices, projector->rIndices, K->decomposition->neighbors)) {
		eslog::error("cannot exchange sIndices of G.\n");
	}

	// build local G
	size_t n = 0;
	esint rows = 0, cols = 0, nnrows = 0, nncols = 0;
	std::vector<size_t> rDomains(K->decomposition->neighbors.size());
	for ( ; n < K->decomposition->neighbors.size() && K->decomposition->neighbors[n] < info::mpi::rank; ++n) {
		for (size_t i = 0; i < projector->rIndices[n].size(); i += projector->rIndices[n][i + 2] + 3) {
			nnrows += projector->rIndices[n][i + 1];
			nncols += projector->rIndices[n][i + 2];
			++rDomains[n];
		}
	}
	projector->nnGPreRows = nnrows;
	for (size_t d = 0; d < K->domains.size(); ++d) {
		rows += R->R1.domains[d].ncols;
		cols += L->D2C[d].size();
	}
	for (; n < K->decomposition->neighbors.size(); ++n) {
		for (size_t i = 0; i < projector->rIndices[n].size(); i += projector->rIndices[n][i + 2] + 3) {
			nnrows += projector->rIndices[n][i + 1];
			nncols += projector->rIndices[n][i + 2];
		}
	}

	projector->nnG.resize(rows + nnrows, L->B1.domains.front().nrows, cols + nncols);
	projector->localG.resize(rows, L->B1.domains.front().nrows, 0); // we use rows only

	n = 0;
	esint row = 0, col = 0;
	projector->nnG.rows[row] = IDX;

	auto push_from_neigh = [&] () {
		std::vector<esint> rindices; rindices.reserve(4 * rDomains[n]);
		for (size_t i = 0; i < projector->rIndices[n].size(); i += projector->rIndices[n][i + 2] + 3) {
			rindices.insert(rindices.end(), projector->rIndices[n].data() + i, projector->rIndices[n].data() + i + 3);
			rindices.push_back(col); // store index to nnG matrix
			for (esint r = 0; r < projector->rIndices[n][i + 1]; ++r) {
				projector->nnGMap[projector->rIndices[n][i] + r] = row;
				memcpy(projector->nnG.cols + col, projector->rIndices[n].data() + i + 3, sizeof(esint) * projector->rIndices[n][i + 2]);
				col += projector->rIndices[n][i + 2];
			}
			projector->nnG.rows[++row] = col + IDX;
		}
		projector->rIndices[n].swap(rindices); // columns indices are irrelevant
	};

	for ( ; n < K->decomposition->neighbors.size() && K->decomposition->neighbors[n] < info::mpi::rank; ++n) {
		push_from_neigh();
	}

	esint lr = 0, lc = 0;
	projector->localG.rows[lr] = IDX;
	projector->localG.cols = projector->nnG.cols + col;
	projector->localG.vals = projector->nnG.vals + col;
	for (size_t i = 0; i < projector->sIndices.size(); ++i) {
		projector->sIndices[i].clear(); // probably memory inefficient (due to push_back behavior)
	}
	for (size_t d = 0; d < K->domains.size(); ++d) {
		for (esint r = 0; r < R->R1.domains[d].ncols; ++r) {
			for (size_t lambda = 0; lambda < L->D2C[d].size(); ++lambda, ++lc, ++col) {
				auto neighs = L->L2MPI->cbegin() + L->D2C[d][lambda];
				for (auto nn = neighs->begin(); nn != neighs->end(); ++nn) {
					projector->sIndices[*nn].push_back(lc);
				}
				projector->nnG.cols[col] = L->C2G[L->D2C[d][lambda]] + IDX;
			}
			projector->nnG.rows[++row] = col + IDX;
			projector->localG.rows[++lr] = lc + IDX;
		}
	}

	for ( ; n < K->decomposition->neighbors.size(); ++n) {
		push_from_neigh();
	}

	for (size_t i = 0; i < projector->sIndices.size(); ++i) {
		projector->sBuffer[i].resize(projector->sIndices[i].size());
		size_t size = 0;
		for (size_t j = 0; j < projector->rIndices[i].size(); j += 4) {
			size += projector->rIndices[i][j + 2];
		}
		projector->rBuffer[i].resize(size);
	}
}

template <typename T>
static void _setGGt(Orthogonal<T> *projector)
{
	const Matrix_FETI<Matrix_CSR, T> *K = projector->feti->K;
	const int IDX = _Matrix_CSR_Pattern::Indexing;

	std::vector<esint> ggt(projector->localG.nrows + 1);
	ggt[0] = projector->localG.nrows;

	// TODO: avoid brute force algorithm below
	for (esint r = 0; r < projector->localG.nrows; ++r) {
		auto nonzero = [&] (esint c) {
			const esint *c1    = -IDX + projector->localG.cols + projector->localG.rows[r];
			const esint *c1end = -IDX + projector->localG.cols + projector->localG.rows[r + 1];
			const esint *c2    = -IDX + projector->nnG.cols + projector->nnG.rows[c];
			const esint *c2end = -IDX + projector->nnG.cols + projector->nnG.rows[c + 1];
			while (c1 != c1end && c2 != c2end && *c1 != *c2) { *c1 < *c2 ? ++c1 : ++c2; }
			return !(c1 == c1end || c2 == c2end || *c1 != *c2);
		};

		size_t n = 0;
		esint c = 0;
		auto check_neigh = [&] () {
			for (size_t i = 0; i < projector->rIndices[n].size(); i += 4) {
				// build upper triangle only
				if (r + projector->feti->sinfo.R1offset <= projector->rIndices[n][i] && nonzero(c)) {
					for (esint rr = 0; rr < projector->rIndices[n][i + 1]; ++rr) {
						ggt.push_back(projector->rIndices[n][i] + IDX);
					}
				}
				c += projector->rIndices[n][i + 1];
			}
		};

		for ( ; n < K->decomposition->neighbors.size() && K->decomposition->neighbors[n] < info::mpi::rank; ++n) {
			check_neigh();
		}
		for (esint rr = 0; rr < projector->localG.nrows; ++rr, ++c) {
			// build upper triangle only
			if (r <= c - projector->nnGPreRows && nonzero(c)) {
				ggt.push_back(projector->feti->sinfo.R1offset + rr  + IDX);
			}
		}
		for (; n < K->decomposition->neighbors.size(); ++n) {
			check_neigh();
		}
		ggt[r + 1] = ggt.size() - (projector->localG.nrows + 1) + IDX;
	}
	esint offset = ggt.size() - (projector->localG.nrows + 1);
	Communication::exscan(&offset, nullptr, 1, MPITools::getType(offset).mpitype, MPI_SUM);
	if (info::mpi::rank == 0) {
		offset = 0;
	}
	for (esint r = 0; r < projector->localG.nrows; ++r) {
		ggt[r + 1] += offset;
	}

	// prepare more optimized version (avoid vector resizing)
	if (!Communication::allGatherUnknownSize(ggt)) {
		eslog::error("cannot gather ggt.\n");
	}

	projector->GGt.shape = Matrix_Shape::UPPER;
	projector->GGt.type = Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE;
	projector->GGt.resize(projector->feti->sinfo.R1size, projector->feti->sinfo.R1size, ggt.size() - info::mpi::size - projector->feti->sinfo.R1size);
	projector->GGt.rows[0] = IDX;
	for (size_t i = 0, r = 1, nnz = 0; i < ggt.size();) {
		const esint &rows = ggt[i], cols = ggt[i + rows] - projector->GGt.rows[r - 1];
		memcpy(projector->GGt.rows + r, ggt.data() + i + 1, sizeof(esint) * rows);
		memcpy(projector->GGt.cols + nnz, ggt.data() + i + rows + 1, sizeof(esint) * cols);
		i += 1 + rows + cols; r += rows; nnz += cols;
	}

	math::commit(projector->GGt);
	math::symbolicFactorization(projector->GGt);

	projector->invGGt.shape = Matrix_Shape::FULL;
	projector->invGGt.type = Matrix_Type::REAL_NONSYMMETRIC;
	projector->invGGt.resize(projector->localG.nrows, projector->feti->sinfo.R1size);
}

template <typename T>
static void _set(Orthogonal<T> *projector)
{
	_setG(projector);
	_setGGt(projector);
}

template <typename T>
static void _update(Orthogonal<T> *projector)
{
	const Matrix_FETI<Matrix_CSR, T> *K = projector->feti->K;
	const typename AX_FETI<T>::Regularization *R = projector->feti->regularization;
	const typename AX_FETI<T>::EqualityConstraints *L = projector->feti->equalityConstraints;
	esint GGtOffset = projector->feti->sinfo.R1offset;
	Matrix_CSR<T> &GGt = projector->GGt;
	const int IDX = _Matrix_CSR_Pattern::Indexing;

	esint row = 0, col = 0;
	for (size_t d = 0; d < K->domains.size(); ++d) {
		for (esint r = 0; r < R->R1.domains[d].ncols; ++r, ++row) {
			esint lr = 0;
			for (size_t lambda = 0; lambda < L->D2C[d].size(); ++lambda, ++col) {
				T value = T{};
				for (; lr < L->B1.domains[d].nnz && (esint)lambda + IDX == L->B1.domains[d].rows[lr]; ++lr) { // change B1 to CSR format
					value += R->R1.domains[d].vals[r * R->R1.domains[d].nrows + L->B1.domains[d].cols[lr] - IDX] * L->B1.domains[d].vals[lr];
				}
				projector->localG.vals[col] = value;
			}
		}
	}

	for (size_t n = 0; n < projector->sIndices.size(); ++n) {
		for (size_t i = 0; i < projector->sIndices[n].size(); ++i) {
			projector->sBuffer[n][i] = projector->localG.vals[projector->sIndices[n][i]];
		}
	}

	if (!Communication::exchangeKnownSize(projector->sBuffer, projector->rBuffer, K->decomposition->neighbors)) {
		eslog::error("cannot exchange nnG values.\n");
	}

	for (size_t n = 0; n < projector->rIndices.size(); ++n) {
		memcpy(projector->nnG.vals + projector->rIndices[n][3], projector->rBuffer[n].data(), sizeof(T) * projector->rBuffer[n].size());
	}

	esint cmin = GGtOffset + IDX;
	esint cmax = cmin + projector->localG.nrows;
	for (esint r = 0; r < projector->localG.nrows; ++r) {
		auto mult = [] (const Orthogonal<T> *projector, const esint &r, const esint &c, T &result) {
			const esint *c1    = -IDX + projector->localG.cols + projector->localG.rows[r];
			const esint *c1end = -IDX + projector->localG.cols + projector->localG.rows[r + 1];
			const esint *c2    = -IDX + projector->nnG.cols + projector->nnG.rows[c];
			const esint *c2end = -IDX + projector->nnG.cols + projector->nnG.rows[c + 1];
			result = T{};
			while (c1 != c1end && c2 != c2end) {
				if (*c1 == *c2) {
					result += projector->localG.vals[c1 - projector->localG.cols] * projector->nnG.vals[c2 - projector->nnG.cols];
					++c1; ++c2;
				} else {
					*c1 < *c2 ? ++c1 : ++c2;
				}
			}
		};

		for (esint c = GGt.rows[r + GGtOffset]; c < GGt.rows[r + GGtOffset + 1]; ++c) {
			if (GGt.cols[c - IDX] < cmin || cmax <= GGt.cols[c - IDX]) {
				mult(projector, r, projector->nnGMap[GGt.cols[c - IDX] - IDX], GGt.vals[c - IDX]);
			} else {
				mult(projector, r, GGt.cols[c - IDX] - cmin + projector->nnGPreRows, GGt.vals[c - IDX]);
			}
		}
	}

	if (!Communication::allGatherInplace(GGt.vals, projector->GGt.rows[GGtOffset] - IDX, projector->GGt.rows[GGtOffset + projector->localG.nrows] - projector->GGt.rows[GGtOffset])) {
		eslog::error("cannot gather GGt values.\n");
	}

	Matrix_Dense<T> eye;
	eye.resize(projector->localG.nrows, projector->feti->sinfo.R1size);
	math::set(eye, T{});
	for (esint r = 0; r < projector->localG.nrows; ++r) {
		eye.vals[r * projector->feti->sinfo.R1size + projector->feti->sinfo.R1offset + r] = T{1};
	}

	math::numericalFactorization(projector->GGt);
	math::solve(projector->GGt, eye, projector->invGGt);

	if (info::ecf->output.print_matrices) {
		eslog::storedata(" STORE: feti/{nnG, localG, GGt}\n");
		math::store(projector->nnG, utils::filename(utils::debugDirectory(*projector->feti->step) + "/feti", "nnG").c_str());
		math::store(projector->localG, utils::filename(utils::debugDirectory(*projector->feti->step) + "/feti", "localG").c_str());
		math::store(projector->GGt, utils::filename(utils::debugDirectory(*projector->feti->step) + "/feti", "GGt").c_str());
		math::store(projector->invGGt, utils::filename(utils::debugDirectory(*projector->feti->step) + "/feti", "invGGt").c_str());
	}
}

template <> Orthogonal<double>::Orthogonal(AX_FETI<double> *feti): Projector(feti) { _set<double>(this); }
template <> Orthogonal<std::complex<double> >::Orthogonal(AX_FETI<std::complex<double> > *feti): Projector(feti) { _set<std::complex<double> >(this); }

template <> void Orthogonal<double>::info() { _info<double>(this); }
template <> void Orthogonal<std::complex<double> >::info() { _info<std::complex<double> >(this); }

template <> void Orthogonal<double>::update() { _update<double>(this); }
template <> void Orthogonal<std::complex<double> >::update() { _update<std::complex<double> >(this); }

}

