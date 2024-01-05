

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

	_computeDualGraph();
	_setG();
	_setGGt();
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
	eslog::info(" =   GGT ROWS                                                                      %9d = \n", GGt.nrows);
	eslog::info(" =   GGT NNZ                                                                       %9d = \n", nnz);
//	eslog::info(" =   GGT FACTORS NNZ                                                               %9d = \n", GGtSolver.nnzL);
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
		_e.vals = e.vals + d;
		Matrix_Dense<T> _Rt;
		_Rt.nrows = R.R1.domains[d].ncols;
		_Rt.ncols = R.R1.domains[d].nrows;
		_Rt.vals = R.R1.domains[d].vals;
		math::apply(_e, T{1}, _Rt, T{0}, feti.f.domains[d]);
	}
	e.synchronize();
	eslog::checkpointln("FETI: COMPUTE DUAL RHS [e]");

	_updateG();
	_updateGGt();

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
		y.vals = in.vals + d;

		math::apply(out.domains[d], T{1}, feti.regularization.R1.domains[d], T{0}, y);
	}
}

template<typename T>
void OrthogonalTFETISymmetric<T>::_computeDualGraph()
{
	const Matrix_FETI<Matrix_CSR, T> &K = feti.K;
	const typename FETI<T>::EqualityConstraints &eq = feti.equalityConstraints;

	dualGraph.resize(K.domains.size());
	for (size_t i = 0; i < eq.cmap.size(); ) {
		esint domains = eq.cmap[i + 1];
		for (esint d1 = 0; d1 < domains; ++d1) {
			for (esint d2 = 0; d2 < domains; ++d2) {
				dualGraph[eq.cmap[i + 2 + d1] - K.decomposition->dbegin].push_back(eq.cmap[i + 2 + d2]);
			}
		}
		i += eq.cmap[i + 1] + 2;
	}

	for (size_t d = 0; d < K.domains.size(); ++d) {
		utils::sortAndRemoveDuplicates(dualGraph[d]);
	}
}

template<typename T>
void OrthogonalTFETISymmetric<T>::_setG()
{
	// G is stored with 0-based in indexing
	const Matrix_FETI<Matrix_CSR, T> &K = feti.K;
	const typename FETI<T>::Regularization &R = feti.regularization;
	const typename FETI<T>::EqualityConstraints &eq = feti.equalityConstraints;

	esint Grows = 0, Gnnz = 0;
	for (size_t d = 0; d < K.domains.size(); ++d) {
		Grows += R.R1.domains[d].ncols;
		Gnnz += R.R1.domains[d].ncols * eq.domain[d].D2C.size();
	}

	Gt.resize(Grows, feti.sinfo.lambdasLocal, Gnnz);
	Gt.rows[0] = 0;
	for (size_t d = 0, ri = 0; d < K.domains.size(); ++d) {
		for (esint r = 0; r < R.R1.domains[d].ncols; ++r, ++ri) {
			Gt.rows[ri + 1] = Gt.rows[ri] + eq.domain[d].B1.nrows;
			for (esint c = 0; c < eq.domain[d].B1.nrows; ++c) {
				Gt.cols[Gt.rows[ri] + c] = eq.domain[d].D2C[c];
			}
		}
	}

	G.shallowCopy(Gt);
	G.nrows = Grows;
	G.nnz = Gnnz;
	eslog::checkpointln("FETI: SET G");
}

template<typename T>
void OrthogonalTFETISymmetric<T>::_updateG()
{
	// G is stored with 0-based in indexing
	const Matrix_FETI<Matrix_CSR, T> &K = feti.K;
	const typename FETI<T>::Regularization &R = feti.regularization;
	const typename FETI<T>::EqualityConstraints &L = feti.equalityConstraints;

	for (size_t d = 0, r = 0; d < K.domains.size(); ++d) {
		for (esint k = 0; k < R.R1.domains[d].ncols; ++k, ++r) {
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
void OrthogonalTFETISymmetric<T>::_setGGt()
{
	const Matrix_FETI<Matrix_CSR, T> &K = feti.K;
	const int IDX = Indexing::CSR;

	GGtDataOffset = 0;
	for (size_t d = 0; d < dualGraph.size(); ++d) {
		for (size_t i = 0; i < dualGraph[d].size(); ++i) {
			if (K.decomposition->dbegin + (esint)d <= dualGraph[d][i]) {
				++GGtDataOffset;
			}
		}
	}
	GGtOffset = G.nrows;
	GGtSize = Communication::exscan(GGtOffset);
	GGtNnz = Communication::exscan(GGtDataOffset);

	GGt.resize(GGtSize, GGtSize, GGtNnz);
	GGt.shape = Matrix_Shape::UPPER;
	GGt.type = Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE;
	GGt.rows[0] = IDX;
	GGt.rows[GGtOffset] = GGtDataOffset + IDX;
	for (size_t d = 0; d < dualGraph.size(); ++d) {
		GGt.rows[GGtOffset + d + 1] = GGt.rows[GGtOffset + d];
		for (size_t i = 0, c = GGt.rows[GGtOffset + d] - IDX; i < dualGraph[d].size(); ++i) {
			if (K.decomposition->dbegin + (esint)d <= dualGraph[d][i]) {
				GGt.cols[c++] = dualGraph[d][i] + IDX;
				++GGt.rows[GGtOffset + d + 1];
			}
		}
	}

	if (!Communication::allGatherInplace(GGt.rows, feti.sinfo.R1offset + 1, G.nrows)) {
		eslog::error("cannot gather GGt rows.\n");
	}
	if (!Communication::allGatherInplace(GGt.cols, GGtOffset, GGtSize)) {
		eslog::error("cannot gather GGt cols.\n");
	}

	invGGt.resize(G.nrows, GGtSize);
	eslog::checkpointln("FETI: SET GGT");
}

template<typename T>
void OrthogonalTFETISymmetric<T>::_updateGGt()
{
	const Matrix_FETI<Matrix_CSR, T> &K = feti.K;
	const int IDX = Indexing::CSR;

	for (size_t d = 0; d < dualGraph.size(); ++d) {
		for (size_t i = 0, c = GGt.rows[GGtOffset + d] - IDX; i < dualGraph[d].size(); ++i) {
			if (K.decomposition->dbegin + (esint)d <= dualGraph[d][i]) {
				GGt.vals[++c - 1] = 0;
				esint k1 = G.rows[d], ke1 = G.rows[d + 1];
				esint k2 = G.rows[dualGraph[d][i]], ke2 = G.rows[dualGraph[d][i] + 1];
				while (k1 < ke1 && k2 < ke2) {
					while (k1 < ke1 && G.cols[k1] < G.cols[k2]) { ++k1; };
					while (k2 < ke2 && G.cols[k2] < G.cols[k1]) { ++k2; };
					if (k1 < ke1 && k2 < ke2 && G.cols[k1] == G.cols[k2]) {
						GGt.vals[c - 1] += G.vals[k1++] * G.vals[k2++];
					}
				}
			}
		}
	}


	if (!Communication::allGatherInplace(GGt.vals, GGtOffset, GGtSize)) {
		eslog::error("cannot gather GGt vals.\n");
	}
	eslog::checkpointln("FETI: GATHER GGT VALUES");

	DirectSolver<T, Matrix_CSR> GGtSolver;
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

