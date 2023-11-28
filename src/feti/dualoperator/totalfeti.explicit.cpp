
#include "totalfeti.explicit.h"
#include "feti/common/applyB.h"

#include "basis/utilities/sysutils.h"
#include "esinfo/eslog.hpp"
#include "esinfo/ecfinfo.h"
#include "wrappers/mpi/communication.h"

namespace espreso {

template <typename T>
TotalFETIExplicit<T>::TotalFETIExplicit(FETI<T> &feti)
: DualOperator<T>(feti)
{

}

template <typename T>
TotalFETIExplicit<T>::~TotalFETIExplicit()
{

}

template <typename T>
void TotalFETIExplicit<T>::info()
{
//	DualOperatorInfo sum, min, max;
	size_t minF = INT32_MAX, maxF = 0, sumF = 0;
	for (size_t d = 0; d < F.size(); ++d) {
		minF = std::min(minF, F[d].nrows * F[d].ncols * sizeof(double));
		maxF = std::max(maxF, F[d].nrows * F[d].ncols * sizeof(double));
		sumF += F[d].nrows * F[d].ncols * sizeof(double);
	}

//	TotalFETIImplicit<T>::reduceInfo(sum, min, max);
	Communication::allReduce(&minF, nullptr, 1, MPITools::getType<size_t>().mpitype, MPI_MIN);
	Communication::allReduce(&maxF, nullptr, 1, MPITools::getType<size_t>().mpitype, MPI_MAX);
	Communication::allReduce(&sumF, nullptr, 1, MPITools::getType<size_t>().mpitype, MPI_SUM);

	eslog::info(" = EXPLICIT TOTAL FETI OPERATOR                                                              = \n");
//	TotalFETIImplicit<T>::printInfo(sum, min, max);
	eslog::info(" =   F MEMORY [MB]                                            %8.2f <%8.2f - %8.2f> = \n", (double)sumF / F.size() / 1024. / 1024., minF / 1024. / 1024., maxF / 1024. / 1024.);
	eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
}

template <typename T>
void TotalFETIExplicit<T>::set(const step::Step &step)
{
	const typename FETI<T>::EqualityConstraints &L = feti.equalityConstraints;
//	sparsity = feti.configuration.partial_dual ? DirectSolver<T, Matrix_CSR>::VectorSparsity::SPARSE_RHS | DirectSolver<T, Matrix_CSR>::VectorSparsity::SPARSE_SOLUTION : DirectSolver<T, Matrix_CSR>::VectorSparsity::DENSE;

	Kplus.resize(feti.K.domains.size());
	d.resize();
	Btx.resize(feti.K.domains.size());
	KplusBtx.resize(feti.K.domains.size());
	KSolver.resize(feti.K.domains.size());

	#pragma omp parallel for
	for (size_t di = 0; di < feti.K.domains.size(); ++di) {
		Kplus[di].type = feti.K.domains[di].type;
		Kplus[di].shape = feti.K.domains[di].shape;
		math::combine(Kplus[di], feti.K.domains[di], feti.regularization.RegMat.domains[di]);
		Btx[di].resize(feti.K.domains[di].nrows);
		KplusBtx[di].resize(feti.K.domains[di].nrows);
		math::set(Btx[di], T{0});
	}
	eslog::checkpointln("FETI: SET TOTAL-FETI OPERATOR");

	#pragma omp parallel for
	for (size_t di = 0; di < feti.K.domains.size(); ++di) {
		KSolver[di].commit(Kplus[di]);

		esint suffix = 0;
//		if (sparsity != DirectSolver<T, Matrix_CSR>::VectorSparsity::DENSE) {
//			suffix = *std::min_element(L.domain[di].B1.cols, L.domain[di].B1.cols + L.domain[di].B1.nnz);
//		}

		KSolver[di].symbolicFactorization(suffix);
	}
	eslog::checkpointln("FETI: TFETI SYMBOLIC FACTORIZATION");

	F.resize(feti.K.domains.size());
	in.resize(feti.K.domains.size());
	out.resize(feti.K.domains.size());

	#pragma omp parallel for
	for (size_t d = 0; d < feti.K.domains.size(); ++d) {
		F[d].resize(L.domain[d].B1.nrows, L.domain[d].B1.nrows);
		in[d].resize(L.domain[d].B1.nrows);
		out[d].resize(L.domain[d].B1.nrows);
	}
	eslog::checkpointln("FETI: TFETI SET EXPLICIT OPERATOR");
}

template <typename T>
void TotalFETIExplicit<T>::update(const step::Step &step)
{
	#pragma omp parallel for
	for (size_t di = 0; di < feti.K.domains.size(); ++di) {
		math::sumCombined(Kplus[di], T{1}, feti.K.domains[di], feti.regularization.RegMat.domains[di]);
	}
	eslog::checkpointln("FETI: UPDATE TOTAL-FETI OPERATOR");
	#pragma omp parallel for
	for (size_t di = 0; di < feti.K.domains.size(); ++di) {
		KSolver[di].numericalFactorization();
	}
	eslog::checkpointln("FETI: TFETI NUMERICAL FACTORIZATION");

	#pragma omp parallel for
	for (size_t di = 0; di < feti.K.domains.size(); ++di) {
		KSolver[di].solve(feti.f.domains[di], KplusBtx[di]);
	}
	applyB(feti, KplusBtx, d);
	d.add(T{-1}, feti.equalityConstraints.c);
	eslog::checkpointln("FETI: COMPUTE DUAL RHS [d]");

	const typename FETI<T>::EqualityConstraints &L = feti.equalityConstraints;

	bool perLambda = false;

	if (perLambda) {
		#pragma omp parallel for
		for (size_t d = 0; d < feti.K.domains.size(); ++d) {
			Vector_Dense<T> KplusBt, Bt;
			KplusBt.resize(L.domain[d].B1.ncols);
			Bt.resize(L.domain[d].B1.ncols);

			for (esint r = 0; r < L.domain[d].B1.nrows; ++r) {
				math::set(Bt, T{0});
				for (esint c = L.domain[d].B1.rows[r]; c < L.domain[d].B1.rows[r + 1]; ++c) {
					Bt.vals[L.domain[d].B1.cols[c]] = L.domain[d].B1.vals[c];
				}
				KSolver[d].solve(Bt, KplusBt);

				for (esint fr = 0; fr < L.domain[d].B1.nrows; ++fr) {
					F[d].vals[fr * L.domain[d].B1.nrows + r] = 0;
					for (esint fc = L.domain[d].B1.rows[fr]; fc < L.domain[d].B1.rows[fr + 1]; ++fc) {
						F[d].vals[fr * L.domain[d].B1.nrows + r] += L.domain[d].B1.vals[fc] * KplusBt.vals[L.domain[d].B1.cols[fc]];
					}
				}
			}
		}
	} else {
		#pragma omp parallel for
		for (size_t d = 0; d < feti.K.domains.size(); ++d) {
			Matrix_Dense<T> KplusBt, Bt;
			KplusBt.resize(L.domain[d].B1.nrows, L.domain[d].B1.ncols);
			Bt.resize(L.domain[d].B1.nrows, L.domain[d].B1.ncols);

			for (esint r = 0; r < L.domain[d].B1.nrows; ++r) {
				esint cc = 0;
				for (esint c = L.domain[d].B1.rows[r]; c < L.domain[d].B1.rows[r + 1]; ++c, ++cc) {
					while (cc < L.domain[d].B1.cols[c]) {
						Bt.vals[r * Bt.ncols + cc] = 0;
						++cc;
					}
					Bt.vals[r * Bt.ncols + cc] = L.domain[d].B1.vals[c];
				}
				while (cc < Bt.ncols) {
					Bt.vals[r * Bt.ncols + cc] = 0;
					++cc;
				}
			}
			KSolver[d].solve(Bt, KplusBt);
			for (esint r = 0; r < L.domain[d].B1.nrows; ++r) {
				for (esint lr = 0; lr < L.domain[d].B1.nrows; ++lr) {
					F[d].vals[lr * L.domain[d].B1.nrows + r] = 0;
					for (esint lc = L.domain[d].B1.rows[lr]; lc < L.domain[d].B1.rows[lr + 1]; ++lc) {
						F[d].vals[lr * F[d].ncols + r] += L.domain[d].B1.vals[lc] * KplusBt.vals[r * KplusBt.ncols + L.domain[d].B1.cols[lc]];
					}
				}
			}
		}
	}

	eslog::checkpointln("FETI: TFETI (B * K+ * B') ASSEMBLED");

	print(step);
}

template <typename T>
void TotalFETIExplicit<T>::apply(const Vector_Dual<T> &x, Vector_Dual<T> &y)
{
	#pragma omp parallel for
	for (size_t d = 0; d < feti.K.domains.size(); ++d) {
		extractDomain(feti, d, x, in[d]);
		math::apply(out[d], T{1}, F[d], T{0}, in[d]);
	}
	insertDomains(feti, out, y);
}

template <typename T>
void TotalFETIExplicit<T>::toPrimal(const Vector_Dual<T> &x, Vector_FETI<Vector_Dense, T> &y)
{
	#pragma omp parallel for
	for (size_t di = 0; di < feti.K.domains.size(); ++di) {
		applyBt(feti, di, x, Btx[di]);
		math::copy(KplusBtx[di], feti.f.domains[di]);
		math::add(KplusBtx[di], T{-1}, Btx[di]);
		KSolver[di].solve(KplusBtx[di], y.domains[di]);
	}
}

template <typename T>
void TotalFETIExplicit<T>::print(const step::Step &step)
{
	if (info::ecf->output.print_matrices) {
		eslog::storedata(" STORE: feti/dualop/{Kplus, F}\n");
		for (size_t di = 0; di < feti.K.domains.size(); ++di) {
			math::store(Kplus[di], utils::filename(utils::debugDirectory(step) + "/feti/dualop", (std::string("Kplus") + std::to_string(di)).c_str()).c_str());
			math::store(F[di], utils::filename(utils::debugDirectory(step) + "/feti/dualop", (std::string("F") + std::to_string(di)).c_str()).c_str());
		}
		math::store(d, utils::filename(utils::debugDirectory(step) + "/feti/dualop", "d").c_str());
	}
}

template class TotalFETIExplicit<double>;
template class TotalFETIExplicit<std::complex<double> >;

}
