
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
//	sparsity = feti.configuration.partial_dual ? DirectSolver<T, Matrix_CSR>::VectorSparsity::SPARSE_RHS | DirectSolver<T, Matrix_CSR>::VectorSparsity::SPARSE_SOLUTION : DirectSolver<T, Matrix_CSR>::VectorSparsity::DENSE;

	Kplus.resize(feti.K.size());
	d.resize(feti.lambdas.size);
	Btx.resize(feti.K.size());
	KplusBtx.resize(feti.K.size());
	KSolver.resize(feti.K.size());

	#pragma omp parallel for
	for (size_t di = 0; di < feti.K.size(); ++di) {
		Kplus[di].type = feti.K[di].type;
		Kplus[di].shape = feti.K[di].shape;
		math::combine(Kplus[di], feti.K[di], feti.RegMat[di]);
		Btx[di].resize(feti.K[di].nrows);
		KplusBtx[di].resize(feti.K[di].nrows);
		math::set(Btx[di], T{0});
	}
	eslog::checkpointln("FETI: SET TOTAL-FETI OPERATOR");

	#pragma omp parallel for
	for (size_t di = 0; di < feti.K.size(); ++di) {
		KSolver[di].commit(Kplus[di]);

		esint suffix = 0;
//		if (sparsity != DirectSolver<T, Matrix_CSR>::VectorSparsity::DENSE) {
//			suffix = *std::min_element(feti.B1[di].cols, feti.B1[di].cols + feti.B1[di].nnz);
//		}

		KSolver[di].symbolicFactorization(suffix);
	}
	eslog::checkpointln("FETI: TFETI SYMBOLIC FACTORIZATION");

	F.resize(feti.K.size());
	in.resize(feti.K.size());
	out.resize(feti.K.size());

	#pragma omp parallel for
	for (size_t d = 0; d < feti.K.size(); ++d) {
		F[d].resize(feti.B1[d].nrows, feti.B1[d].nrows);
		in[d].resize(feti.B1[d].nrows);
		out[d].resize(feti.B1[d].nrows);
	}
	eslog::checkpointln("FETI: TFETI SET EXPLICIT OPERATOR");
}

template <typename T>
void TotalFETIExplicit<T>::update(const step::Step &step)
{
	if (feti.updated.K) {
		#pragma omp parallel for
		for (size_t di = 0; di < feti.K.size(); ++di) {
			math::sumCombined(Kplus[di], T{1}, feti.K[di], feti.RegMat[di]);
		}
		eslog::checkpointln("FETI: UPDATE TOTAL-FETI OPERATOR");
		#pragma omp parallel for
		for (size_t di = 0; di < feti.K.size(); ++di) {
			KSolver[di].numericalFactorization();
		}
		eslog::checkpointln("FETI: TFETI NUMERICAL FACTORIZATION");
	}

	#pragma omp parallel for
	for (size_t di = 0; di < feti.K.size(); ++di) {
		KSolver[di].solve(feti.f[di], KplusBtx[di]);
	}
	applyB(feti, KplusBtx, d);
	d.add(T{-1}, feti.c);
	eslog::checkpointln("FETI: COMPUTE DUAL RHS [d]");

	bool perLambda = false;

	if (feti.updated.K) {
		if (perLambda) {
			#pragma omp parallel for
			for (size_t d = 0; d < feti.K.size(); ++d) {
				Vector_Dense<T> KplusBt, Bt;
				KplusBt.resize(feti.B1[d].ncols);
				Bt.resize(feti.B1[d].ncols);

				for (esint r = 0; r < feti.B1[d].nrows; ++r) {
					math::set(Bt, T{0});
					for (esint c = feti.B1[d].rows[r]; c < feti.B1[d].rows[r + 1]; ++c) {
						Bt.vals[feti.B1[d].cols[c]] = feti.B1[d].vals[c];
					}
					KSolver[d].solve(Bt, KplusBt);

					for (esint fr = 0; fr < feti.B1[d].nrows; ++fr) {
						F[d].vals[fr * feti.B1[d].nrows + r] = 0;
						for (esint fc = feti.B1[d].rows[fr]; fc < feti.B1[d].rows[fr + 1]; ++fc) {
							F[d].vals[fr * feti.B1[d].nrows + r] += feti.B1[d].vals[fc] * KplusBt.vals[feti.B1[d].cols[fc]];
						}
					}
				}
			}
		} else {
			#pragma omp parallel for
			for (size_t d = 0; d < feti.K.size(); ++d) {
				Matrix_Dense<T> KplusBt, Bt;
				KplusBt.resize(feti.B1[d].nrows, feti.B1[d].ncols);
				Bt.resize(feti.B1[d].nrows, feti.B1[d].ncols);

				for (esint r = 0; r < feti.B1[d].nrows; ++r) {
					esint cc = 0;
					for (esint c = feti.B1[d].rows[r]; c < feti.B1[d].rows[r + 1]; ++c, ++cc) {
						while (cc < feti.B1[d].cols[c]) {
							Bt.vals[r * Bt.ncols + cc] = 0;
							++cc;
						}
						Bt.vals[r * Bt.ncols + cc] = feti.B1[d].vals[c];
					}
					while (cc < Bt.ncols) {
						Bt.vals[r * Bt.ncols + cc] = 0;
						++cc;
					}
				}
				KSolver[d].solve(Bt, KplusBt);
				for (esint r = 0; r < feti.B1[d].nrows; ++r) {
					for (esint lr = 0; lr < feti.B1[d].nrows; ++lr) {
						F[d].vals[lr * feti.B1[d].nrows + r] = 0;
						for (esint lc = feti.B1[d].rows[lr]; lc < feti.B1[d].rows[lr + 1]; ++lc) {
							F[d].vals[lr * F[d].ncols + r] += feti.B1[d].vals[lc] * KplusBt.vals[r * KplusBt.ncols + feti.B1[d].cols[lc]];
						}
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
	for (size_t d = 0; d < feti.K.size(); ++d) {
		extractDomain(feti, d, x, in[d]);
		math::blas::apply(out[d], T{1}, F[d], T{0}, in[d]);
	}
	insertDomains(feti, out, y);
}

template <typename T>
void TotalFETIExplicit<T>::toPrimal(const Vector_Dual<T> &x, std::vector<Vector_Dense<T> > &y)
{
	#pragma omp parallel for
	for (size_t di = 0; di < feti.K.size(); ++di) {
		applyBt(feti, di, x, Btx[di]);
		math::copy(KplusBtx[di], feti.f[di]);
		math::add(KplusBtx[di], T{-1}, Btx[di]);
		KSolver[di].solve(KplusBtx[di], y[di]);
	}
}

template <typename T>
void TotalFETIExplicit<T>::print(const step::Step &step)
{
	if (info::ecf->output.print_matrices) {
		eslog::storedata(" STORE: feti/dualop/{Kplus, F}\n");
		for (size_t di = 0; di < feti.K.size(); ++di) {
			math::store(Kplus[di], utils::filename(utils::debugDirectory(step) + "/feti/dualop", (std::string("Kplus") + std::to_string(di)).c_str()).c_str());
			math::store(F[di], utils::filename(utils::debugDirectory(step) + "/feti/dualop", (std::string("F") + std::to_string(di)).c_str()).c_str());
		}
		math::store(d, utils::filename(utils::debugDirectory(step) + "/feti/dualop", "d").c_str());
	}
}

template class TotalFETIExplicit<double>;
template class TotalFETIExplicit<std::complex<double> >;

}
