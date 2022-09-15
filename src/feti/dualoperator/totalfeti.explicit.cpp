
#include "totalfeti.explicit.h"
#include "feti/common/applyB.h"

#include "basis/utilities/sysutils.h"
#include "esinfo/eslog.hpp"
#include "esinfo/ecfinfo.h"

namespace espreso {

template class TotalFETIExplicit<double>;
template class TotalFETIExplicit<std::complex<double> >;

template <typename T>
TotalFETIExplicit<T>::TotalFETIExplicit(FETI<T> *feti)
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
	TotalFETIImplicit<T>::reduceInfo(sum, min, max);

	eslog::info(" = EXPLICIT TOTAL FETI OPERATOR                                                              = \n");
	TotalFETIImplicit<T>::printInfo(sum, min, max);
	eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
}

template <typename T>
void TotalFETIExplicit<T>::set()
{
	TotalFETIImplicit<T>::set();

	this->F.resize(this->feti->K->domains.size());
	this->in.resize(this->feti->K->domains.size());
	this->out.resize(this->feti->K->domains.size());

	const typename FETI<T>::EqualityConstraints *L = this->feti->equalityConstraints;
	#pragma omp parallel for
	for (size_t d = 0; d < this->feti->K->domains.size(); ++d) {
		this->F[d].resize(L->domain[d].B1.nrows, L->domain[d].B1.nrows);
		this->in[d].resize(L->domain[d].B1.nrows);
		this->out[d].resize(L->domain[d].B1.nrows);
	}
	eslog::checkpointln("FETI: TFETI SET EXPLICIT OPERATOR");
}

template <typename T>
void TotalFETIExplicit<T>::update()
{
	TotalFETIImplicit<T>::update();

	const typename FETI<T>::EqualityConstraints *L = this->feti->equalityConstraints;
	
	#pragma omp parallel for
	for (size_t d = 0; d < this->feti->K->domains.size(); ++d) {
		Vector_Dense<T> KplusBt, Bt;
		KplusBt.resize(L->domain[d].B1.ncols);
		Bt.resize(L->domain[d].B1.ncols);

		for (esint r = 0; r < L->domain[d].B1.nrows; ++r) {
			math::set(Bt, T{0});
			for (esint c = L->domain[d].B1.rows[r]; c < L->domain[d].B1.rows[r + 1]; ++c) {
				Bt.vals[L->domain[d].B1.cols[c]] = L->domain[d].B1.vals[c];
			}
			math::solve(this->Kplus[d], Bt, KplusBt, this->sparsity);

			for (esint fr = 0; fr < L->domain[d].B1.nrows; ++fr) {
				this->F[d].vals[fr * L->domain[d].B1.nrows + r] = 0;
				for (esint fc = L->domain[d].B1.rows[fr]; fc < L->domain[d].B1.rows[fr + 1]; ++fc) {
					this->F[d].vals[fr * L->domain[d].B1.nrows + r] += L->domain[d].B1.vals[fc] * KplusBt.vals[L->domain[d].B1.cols[fc]];
				}
			}
		}
	}

	eslog::checkpointln("FETI: TFETI (B * K+ * B') ASSEMBLED");

	if (info::ecf->output.print_matrices) {
		eslog::storedata(" STORE: feti/dualop/{F}\n");
		for (size_t d = 0; d < this->feti->K->domains.size(); ++d) {
			math::store(this->F[d], utils::filename(utils::debugDirectory(*this->feti->step) + "/feti/dualop", (std::string("F") + std::to_string(d)).c_str()).c_str());
		}
	}
}

template <typename T>
void TotalFETIExplicit<T>::apply(const Vector_Dual<T> &x, Vector_Dual<T> &y)
{
	#pragma omp parallel for
	for (size_t d = 0; d < this->feti->K->domains.size(); ++d) {
		extractDomain(this->feti, d, x, this->in[d]);
		math::apply(this->out[d], T{1}, this->F[d], T{0}, this->in[d]);
	}
	insertDomains(this->feti, this->out, y);
}

template <typename T>
void TotalFETIExplicit<T>::toPrimal(const Vector_Dual<T> &x, Vector_FETI<Vector_Dense, T> &y)
{
	TotalFETIImplicit<T>::toPrimal(x, y);
}

}
