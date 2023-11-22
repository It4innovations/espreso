
#include "totalfeti.explicit.h"
#include "feti/common/applyB.h"

#include "basis/utilities/sysutils.h"
#include "esinfo/eslog.hpp"
#include "esinfo/ecfinfo.h"

namespace espreso {

template <typename T>
TotalFETIExplicit<T>::TotalFETIExplicit(FETI<T> &feti)
: TotalFETIImplicit<T>(feti)
{

}

template <typename T>
TotalFETIExplicit<T>::~TotalFETIExplicit()
{

}

template <typename T>
void TotalFETIExplicit<T>::info()
{
	DualOperatorInfo sum, min, max;
	size_t minF = INT32_MAX, maxF = 0, sumF = 0;
	for (size_t d = 0; d < F.size(); ++d) {
		minF = std::min(minF, F[d].nrows * F[d].ncols * sizeof(double));
		maxF = std::max(maxF, F[d].nrows * F[d].ncols * sizeof(double));
		sumF += F[d].nrows * F[d].ncols * sizeof(double);
	}

	TotalFETIImplicit<T>::reduceInfo(sum, min, max);
	Communication::allReduce(&minF, nullptr, 1, MPITools::getType<size_t>().mpitype, MPI_MIN);
	Communication::allReduce(&maxF, nullptr, 1, MPITools::getType<size_t>().mpitype, MPI_MAX);
	Communication::allReduce(&sumF, nullptr, 1, MPITools::getType<size_t>().mpitype, MPI_SUM);

	eslog::info(" = EXPLICIT TOTAL FETI OPERATOR                                                              = \n");
	TotalFETIImplicit<T>::printInfo(sum, min, max);
	eslog::info(" =   F MEMORY [MB]                                            %8.2f <%8.2f - %8.2f> = \n", (double)sumF / F.size() / 1024. / 1024., minF / 1024. / 1024., maxF / 1024. / 1024.);
	eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
}

template <typename T>
void TotalFETIExplicit<T>::set(const step::Step &step)
{
	TotalFETIImplicit<T>::set(step);

	F.resize(feti.K.domains.size());
	in.resize(feti.K.domains.size());
	out.resize(feti.K.domains.size());

	const typename FETI<T>::EqualityConstraints &L = feti.equalityConstraints;
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
	TotalFETIImplicit<T>::update(step);

	const typename FETI<T>::EqualityConstraints &L = feti.equalityConstraints;
	
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
			KSolver[d].solve(Bt, KplusBt, sparsity);

			for (esint fr = 0; fr < L.domain[d].B1.nrows; ++fr) {
				F[d].vals[fr * L.domain[d].B1.nrows + r] = 0;
				for (esint fc = L.domain[d].B1.rows[fr]; fc < L.domain[d].B1.rows[fr + 1]; ++fc) {
					F[d].vals[fr * L.domain[d].B1.nrows + r] += L.domain[d].B1.vals[fc] * KplusBt.vals[L.domain[d].B1.cols[fc]];
				}
			}
		}
	}

	eslog::checkpointln("FETI: TFETI (B * K+ * B') ASSEMBLED");

	if (info::ecf->output.print_matrices) {
		eslog::storedata(" STORE: feti/dualop/{F}\n");
		for (size_t d = 0; d < feti.K.domains.size(); ++d) {
			math::store(F[d], utils::filename(utils::debugDirectory(step) + "/feti/dualop", (std::string("F") + std::to_string(d)).c_str()).c_str());
		}
	}
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
	TotalFETIImplicit<T>::toPrimal(x, y);
}

template class TotalFETIExplicit<double>;
template class TotalFETIExplicit<std::complex<double> >;

}
