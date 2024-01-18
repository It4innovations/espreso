
#include "tfeti.orthogonal.symmetric.h"
#include "basis/containers/serializededata.h"
#include "basis/utilities/sysutils.h"
#include "basis/utilities/utils.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/envinfo.h"
#include "esinfo/eslog.hpp"

#include <vector>
#include <unordered_map>

namespace espreso {

template<typename T>
TFETIOrthogonalSymmetric<T>::TFETIOrthogonalSymmetric(FETI<T> &feti)
: Projector<T>(feti)
{
	const typename FETI<T>::Regularization &R = feti.regularization;

	e.resize();
	Gx.resize();
	iGGtGx.resize(feti.sinfo.R1size);

	domainOffset = feti.K.decomposition->dbegin;

	dinfo.reserve(R.R1.size());
	for (size_t d = 0, koffset = feti.sinfo.R1offset; d < R.R1.size(); ++d) {
		dinfo.push_back(DomainInfo((esint)(domainOffset + d), koffset, R.R1[d].nrows));
		koffset += R.R1[d].nrows;
	}

	_computeDualGraph();
	_setG();
	_setGGt();
}

template<typename T>
TFETIOrthogonalSymmetric<T>::~TFETIOrthogonalSymmetric()
{

}

template<typename T>
void TFETIOrthogonalSymmetric<T>::info()
{
	esint nnz = 2 * (GGt.nnz - GGt.nrows) + GGt.nrows;

	eslog::info(" = ORTHOGONAL PROJECTOR PROPERTIES                                                           = \n");
	eslog::info(" =   GGT ROWS                                                                      %9d = \n", GGt.nrows);
	eslog::info(" =   GGT NNZ                                                                       %9d = \n", nnz);
//	eslog::info(" =   GGT FACTORS NNZ                                                               %9d = \n", GGtSolver.nnzL);
	if (feti.configuration.exhaustive_info) {
		// PPt = eye
	}
	eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
}

template<typename T>
void TFETIOrthogonalSymmetric<T>::update(const step::Step &step)
{
	typename FETI<T>::Regularization &R = feti.regularization;

	#pragma omp parallel for
	for (size_t d = 0; d < dinfo.size(); ++d) {
		math::orthonormalize(R.R1[d]);
		Vector_Dense<T> _e;
		_e.size = R.R1[d].nrows;
		_e.vals = e.vals + dinfo[d].koffset;
		math::blas::apply(_e, T{1}, R.R1[d], T{0}, feti.f.domains[d]);
	}
	e.synchronize();
	eslog::checkpointln("FETI: COMPUTE DUAL RHS [e]");

	_updateG();
	_updateGGt();
	_print(step);
}

template<typename T>
void TFETIOrthogonalSymmetric<T>::apply(const Vector_Dual<T> &x, Vector_Dual<T> &y)
{
	x.copyToWithoutHalo(y);
	_applyG(x, Gx);
	_applyInvGGt(Gx, iGGtGx);
	_applyGt(iGGtGx, T{-1}, y);
}

template<typename T>
void TFETIOrthogonalSymmetric<T>::applyGtInvGGt(const Vector_Kernel<T> &x, Vector_Dual<T> &y)
{
	math::set(y, T{0});
	_applyInvGGt(x, iGGtGx);
	_applyGt(iGGtGx, T{-1}, y);
}

template<typename T>
void TFETIOrthogonalSymmetric<T>::applyRInvGGtG(const Vector_Dual<T> &x, Vector_FETI<Vector_Dense, T> &y)
{
	_applyG(x, Gx);
	_applyInvGGt(Gx, iGGtGx);
	_applyR(iGGtGx, y);
}

template<typename T>
void TFETIOrthogonalSymmetric<T>::_applyG(const Vector_Dual<T> &in, Vector_Kernel<T> &out)
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
void TFETIOrthogonalSymmetric<T>::_applyInvGGt(const Vector_Kernel<T> &in, Vector_Dense<T> &out)
{
	#pragma omp parallel for
	for (int t = 0; t < info::env::threads; ++t) {
		Matrix_Dense<T> a;
		Vector_Dense<T> y;
		a.ncols = invGGt.ncols;
		a.nrows = y.size = Vector_Kernel<T>::distribution[t + 1] - Vector_Kernel<T>::distribution[t];

		a.vals = invGGt.vals + invGGt.ncols * Vector_Kernel<T>::distribution[t];
		y.vals = out.vals + Vector_Kernel<T>::distribution[t];

		math::blas::apply(y, T{1}, a, T{0}, in);
	}
}

template<typename T>
void TFETIOrthogonalSymmetric<T>::_applyGt(const Vector_Dense<T> &in, const T &alpha, Vector_Dual<T> &out)
{
	for (esint r = 0; r < G.nrows; ++r) {
		for (esint c = G.rows[r]; c < G.rows[r + 1]; ++c) {
			out.vals[G.cols[c]] += alpha * G.vals[c] * in.vals[r];
		}
	}
	out.synchronize();
}

template<typename T>
void TFETIOrthogonalSymmetric<T>::_applyR(const Vector_Dense<T> &in, Vector_FETI<Vector_Dense, T> &out)
{
	#pragma omp parallel for
	for (size_t d = 0; d < out.domains.size(); ++d) {
		Vector_Dense<T> y;
		y.size = dinfo[d].kernels;
		y.vals = in.vals + dinfo[d].koffset - feti.sinfo.R1offset;

		math::blas::applyT(out.domains[d], T{1}, feti.regularization.R1[d], T{0}, y);
	}
}

template<typename T>
void TFETIOrthogonalSymmetric<T>::_computeDualGraph()
{
	const DOFsDecomposition *decomposition = feti.decomposition;
	const typename FETI<T>::EqualityConstraints &eq = feti.equalityConstraints;

	dualGraph.resize(dinfo.size());
	for (size_t i = 0; i < eq.cmap.size(); ) {
		esint domains = eq.cmap[i + 1];
		for (esint d1 = 0; d1 < domains; ++d1) {
			esint di1 = eq.cmap[i + 2 + d1] - domainOffset;
			for (esint d2 = 0; d2 < domains; ++d2) {
				esint di2 = eq.cmap[i + 2 + d2] - domainOffset;
				if (decomposition->ismy(eq.cmap[i + 2 + d1]) && decomposition->ismy(eq.cmap[i + 2 + d2])) {
					dualGraph[di1].push_back(dinfo[di2]);
				} else {
					if (decomposition->ismy(eq.cmap[i + 2 + d1]) && !decomposition->ismy(eq.cmap[i + 2 + d2])) {
						dualGraph[di1].push_back(DomainInfo(eq.cmap[i + 2 + d2], 0, 0));
					}
					if (decomposition->ismy(eq.cmap[i + 2 + d2]) && !decomposition->ismy(eq.cmap[i + 2 + d1])) {
						dualGraph[di2].push_back(DomainInfo(eq.cmap[i + 2 + d1], 0, 0));
					}
				}
			}
		}
		i += eq.cmap[i + 1] + 2;
	}
	for (size_t d = 0; d < dinfo.size(); ++d) {
		utils::sortAndRemoveDuplicates(dualGraph[d]);
	}

	std::vector<std::vector<DomainInfo> > sBuffer(decomposition->neighbors.size()), rBuffer(decomposition->neighbors.size());
	for (size_t d = 0; d < dualGraph.size(); ++d) {
		esint last = -1;
		for (size_t di = 0; di < dualGraph[d].size(); ++di) {
			esint n = decomposition->noffset(dualGraph[d][di].domain);
			if (!decomposition->ismy(dualGraph[d][di].domain) && last < n) {
				sBuffer[n].push_back(dinfo[d]);
				last = n;
			}
		}
	}

	if (!Communication::exchangeUnknownSize(sBuffer, rBuffer, decomposition->neighbors)) {
		eslog::error("cannot exchange dual graph info\n");
	}

	std::unordered_map<esint, DomainInfo> other;
	for (size_t n = 0; n < rBuffer.size(); ++n) {
		for (size_t i = 0; i < rBuffer[n].size(); ++i) {
			other[rBuffer[n][i].domain] = rBuffer[n][i];
		}
	}

	for (size_t d = 0; d < dinfo.size(); ++d) {
		for (size_t i = 0; i < dualGraph[d].size(); ++i) {
			if (!decomposition->ismy(dualGraph[d][i].domain)) {
				dualGraph[d][i] = other[dualGraph[d][i].domain];
			}
		}
	}

	downinfo.resize(decomposition->neighbors.size());
	std::vector<esint> cOffset(dinfo.size());
	for (size_t i = 0, offset = 0; i < eq.cmap.size(); ) {
		esint lambdas = eq.cmap[i];
		esint domains = eq.cmap[i + 1];
		esint last = -1;
		for (esint d = 0; d < domains; ++d) {
			esint di = eq.cmap[i + 2 + d] - domainOffset;
			if (eq.cmap[i + 2 + d] < decomposition->dbegin) {
				esint n = decomposition->noffset(eq.cmap[i + 2 + d]);
				if (last < n) {
					for (esint d2 = d; d2 < domains; ++d2) {
						if (decomposition->ismy(eq.cmap[i + 2 + d2])) {
							esint di2 = eq.cmap[i + 2 + d2] - domainOffset;
							downinfo[n][eq.cmap[i + 2 + d2]] = dinfo[eq.cmap[i + 2 + d2] - domainOffset];
							downinfo[n][eq.cmap[i + 2 + d2]].cindices.push_back({cOffset[di2], lambdas});
							downinfo[n][eq.cmap[i + 2 + d2]].ncols += lambdas;
						}
					}
					last = n;
				}
			}
			if (decomposition->dend <= eq.cmap[i + 2 + d]) {
				upinfo[eq.cmap[i + 2 + d]] = other[eq.cmap[i + 2 + d]];
				upinfo[eq.cmap[i + 2 + d]].cindices.push_back({(esint)offset, lambdas});
				upinfo[eq.cmap[i + 2 + d]].ncols += lambdas;
			}
			if (decomposition->ismy(eq.cmap[i + 2 + d])) {
				cOffset[di] += lambdas;
			}
		}
		i += eq.cmap[i + 1] + 2;
		offset += lambdas;
	}
}

template<typename T>
void TFETIOrthogonalSymmetric<T>::_setG()
{
	// G is stored with 0-based in indexing
	const typename FETI<T>::EqualityConstraints &eq = feti.equalityConstraints;

	esint Grows = 0, Gnnz = 0;
	for (size_t d = 0; d < dinfo.size(); ++d) {
		Grows += dinfo[d].kernels;
		Gnnz += dinfo[d].kernels * eq.domain[d].D2C.size();
	}
	esint Gtrows = Grows, Gtnnz = Gnnz;
	for (auto di = upinfo.cbegin(); di != upinfo.cend(); ++di) {
		Gtrows += di->second.kernels;
		Gtnnz += di->second.kernels * di->second.ncols;
	}

	Gt.resize(Gtrows, feti.sinfo.lambdasLocal, Gtnnz);
	Gt.rows[0] = 0;
	size_t ri = 0;
	for (size_t d = 0; d < dinfo.size(); ++d) {
		for (esint kr = 0; kr < dinfo[d].kernels; ++kr, ++ri) {
			Gt.rows[ri + 1] = Gt.rows[ri] + eq.domain[d].B1.nrows;
			for (esint c = 0; c < eq.domain[d].B1.nrows; ++c) {
				Gt.cols[Gt.rows[ri] + c] = eq.domain[d].D2C[c];
			}
		}
	}
	for (auto di = upinfo.begin(); di != upinfo.end(); ++di) {
		NeighborDomainInfo &ndi = di->second;
		ndi.koffset = ri;
		for (esint kr = 0; kr < ndi.kernels; ++kr, ++ri) {
			Gt.rows[ri + 1] = Gt.rows[ri] + ndi.ncols;
			for (size_t ci = 0, c = 0; ci < ndi.cindices.size(); ++ci) {
				for (esint cc = 0; cc < ndi.cindices[ci].count; ++cc, ++c) {
					Gt.cols[Gt.rows[ri] + c] = ndi.cindices[ci].offset + cc;
				}
			}
		}
	}


	G.shallowCopy(Gt);
	G.nrows = Grows;
	G.nnz = Gnnz;
	eslog::checkpointln("FETI: SET G");
}

template<typename T>
void TFETIOrthogonalSymmetric<T>::_updateG()
{
	// G is stored with 0-based in indexing
	const typename FETI<T>::Regularization &R = feti.regularization;
	const typename FETI<T>::EqualityConstraints &L = feti.equalityConstraints;

	for (size_t d = 0, r = 0; d < dinfo.size(); ++d) {
		for (esint kr = 0; kr < dinfo[d].kernels; ++kr, ++r) {
			for (esint c = 0; c < L.domain[d].B1.nrows; ++c) {
				G.vals[G.rows[r] + c] = 0;
				for (esint i = L.domain[d].B1.rows[c]; i < L.domain[d].B1.rows[c + 1]; ++i) {
					G.vals[G.rows[r] + c] -= R.R1[d].vals[R.R1[d].ncols * kr + L.domain[d].B1.cols[i]] * L.domain[d].B1.vals[i];
				}
			}
		}
	}

	const DOFsDecomposition *decomposition = feti.decomposition;

	std::vector<std::vector<T> > sBuffer(decomposition->neighbors.size()), rBuffer(decomposition->neighbors.size());
	for (size_t n = 0; n < downinfo.size(); ++n) {
		for (auto di = downinfo[n].cbegin(); di != downinfo[n].cend(); ++di) {
			const NeighborDomainInfo &ndi = di->second;
			for (esint kr = 0; kr < ndi.kernels; ++kr) {
				esint roffset = G.rows[ndi.koffset - feti.sinfo.R1offset + kr];
				for (size_t cc = 0; cc < ndi.cindices.size(); ++cc) {
					for (esint c = 0; c < ndi.cindices[cc].count; ++c) {
						sBuffer[n].push_back(G.vals[roffset + ndi.cindices[cc].offset + c]);
					}
				}
			}
		}
	}

	if (!Communication::exchangeUnknownSize(sBuffer, rBuffer, decomposition->neighbors)) {
		eslog::error("cannot exchange neighbor's G\n");
	}

	for (size_t n = 0, offset = G.rows[G.nrows]; n < rBuffer.size(); ++n) {
		std::copy(rBuffer[n].begin(), rBuffer[n].end(), Gt.vals + offset);
		offset += rBuffer[n].size();
	}

	eslog::checkpointln("FETI: UPDATE G");
}

template<typename T>
void TFETIOrthogonalSymmetric<T>::_setGGt()
{
	const int IDX = Indexing::CSR;

	GGtDataOffset = 0;
	for (size_t d = 0; d < dinfo.size(); ++d) {
		for (esint kr = 0; kr < dinfo[d].kernels; ++kr) {
			for (size_t i = 0; i < dualGraph[d].size(); ++i) {
				for (esint kc = 0; kc < dualGraph[d][i].kernels; ++kc) {
					if (dinfo[d].koffset + kr <= dualGraph[d][i].koffset + kc) {
						++GGtDataOffset;
					}
				}
			}
		}
	}
	GGtDataSize = GGtDataOffset;
	GGtNnz = Communication::exscan(GGtDataOffset);

	GGt.resize(feti.sinfo.R1totalSize, feti.sinfo.R1totalSize, GGtNnz);
	GGt.shape = Matrix_Shape::UPPER;
	GGt.type = Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE;
	GGt.rows[0] = IDX;
	GGt.rows[feti.sinfo.R1offset] = GGtDataOffset + IDX;
	for (size_t d = 0; d < dinfo.size(); ++d) {
		for (esint kr = 0; kr < dinfo[d].kernels; ++kr) {
			GGt.rows[dinfo[d].koffset + kr + 1] = GGt.rows[dinfo[d].koffset + kr];
			for (size_t i = 0, c = GGt.rows[dinfo[d].koffset + kr] - IDX; i < dualGraph[d].size(); ++i) {
				for (esint kc = 0; kc < dualGraph[d][i].kernels; ++kc) {
					if (dinfo[d].koffset + kr <= dualGraph[d][i].koffset + kc) {
						GGt.cols[c++] = dualGraph[d][i].koffset + kc + IDX;
						++GGt.rows[dinfo[d].koffset + kr + 1];
					}
				}
			}
		}
	}

	if (!Communication::allGatherInplace(GGt.rows, feti.sinfo.R1offset + 1, G.nrows)) {
		eslog::error("cannot gather GGt rows.\n");
	}
	if (!Communication::allGatherInplace(GGt.cols, GGtDataOffset, GGtDataSize)) {
		eslog::error("cannot gather GGt cols.\n");
	}

	invGGt.resize(G.nrows, feti.sinfo.R1totalSize);
	eslog::checkpointln("FETI: SET GGT");
}

template<typename T>
void TFETIOrthogonalSymmetric<T>::_updateGGt()
{
	const int IDX = Indexing::CSR;

	for (size_t d = 0; d < dinfo.size(); ++d) {
		for (esint kr = 0; kr < dinfo[d].kernels; ++kr) {
			for (size_t i = 0, c = GGt.rows[dinfo[d].koffset + kr] - IDX; i < dualGraph[d].size(); ++i) {
				for (esint kc = 0; kc < dualGraph[d][i].kernels; ++kc) {
					if (dinfo[d].koffset + kr <= dualGraph[d][i].koffset + kc) {
						GGt.vals[c] = 0;
						esint k1, k2, ke1, ke2;
						k1  = G.rows[dinfo[d].koffset - feti.sinfo.R1offset + kr];
						ke1 = G.rows[dinfo[d].koffset - feti.sinfo.R1offset + kr + 1];
						if (dualGraph[d][i].koffset - feti.sinfo.R1offset < G.nrows) {
							k2  = G.rows[dualGraph[d][i].koffset - feti.sinfo.R1offset + kc];
							ke2 = G.rows[dualGraph[d][i].koffset - feti.sinfo.R1offset + kc + 1];
						} else {
							k2  = Gt.rows[upinfo[dualGraph[d][i].domain].koffset + kc];
							ke2 = Gt.rows[upinfo[dualGraph[d][i].domain].koffset + kc + 1];
						}
						while (k1 < ke1 && k2 < ke2) {
							while (k1 < ke1 && Gt.cols[k1] < Gt.cols[k2]) { ++k1; };
							while (k2 < ke2 && Gt.cols[k2] < Gt.cols[k1]) { ++k2; };
							if (k1 < ke1 && k2 < ke2 && Gt.cols[k1] == Gt.cols[k2]) {
								GGt.vals[c] += Gt.vals[k1++] * Gt.vals[k2++];
							}
						}
						++c;
					}
				}
			}
		}
	}


	if (!Communication::allGatherInplace(GGt.vals, GGtDataOffset, GGtDataSize)) {
		eslog::error("cannot gather GGt vals.\n");
	}
	eslog::checkpointln("FETI: GATHER GGT VALUES");

	DirectSolver<Matrix_CSR, T> GGtSolver;
	GGtSolver.commit(GGt);
	GGtSolver.symbolicFactorization();
	GGtSolver.numericalFactorization();
	eslog::checkpointln("FETI: GGT FACTORIZATION");

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
void TFETIOrthogonalSymmetric<T>::_print(const step::Step &step)
{
	if (info::ecf->output.print_matrices) {
		eslog::storedata(" STORE: feti/projector/{R_orth, G, e, GGt, invGGt}\n");
		for (size_t di = 0; di < feti.regularization.R1.size(); ++di) {
			math::store(feti.regularization.R1[di], utils::filename(utils::debugDirectory(step) + "/feti/projector", (std::string("R_orth") + std::to_string(di)).c_str()).c_str());
		}
		math::store(G, utils::filename(utils::debugDirectory(step) + "/feti/projector", "G").c_str());
		math::store(Gt, utils::filename(utils::debugDirectory(step) + "/feti/projector", "Gt").c_str());
		math::store(e, utils::filename(utils::debugDirectory(step) + "/feti/projector", "e").c_str());
		math::store(GGt, utils::filename(utils::debugDirectory(step) + "/feti/projector", "GGt").c_str());
		math::store(invGGt, utils::filename(utils::debugDirectory(step) + "/feti/projector", "invGGt").c_str());
	}
}

template struct TFETIOrthogonalSymmetric<double>;
template struct TFETIOrthogonalSymmetric<std::complex<double> >;

}
