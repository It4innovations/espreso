
#include "totalfeti.h"

#include "basis/utilities/sysutils.h"
#include "esinfo/eslog.hpp"
#include "esinfo/ecfinfo.h"

namespace espreso {

template <typename T> static void _applyBt(TotalFETI<T> *dual, size_t d, const Vector_Dual<T> &in, Vector_Dense<T> &out);
template <typename T> static void _applyB(TotalFETI<T> *dual, const std::vector<Vector_Dense<T> > &in, Vector_Dual<T> &out);
template <typename T> static void _print(TotalFETI<T> *dual);

template <typename T>
static void _info(TotalFETI<T> *dual)
{
	eslog::info(" = TOTAL FETI OPERATOR PROPERTIES                                                            = \n");
	eslog::info(" =   DOMAINS                                                                       %9d = \n", dual->feti->sinfo.domains);
	eslog::info(" =   DUAL SIZE                                                                     %9d = \n", dual->feti->sinfo.lambdasTotal);
	if (dual->feti->configuration.exhaustive_info) {
		// power method to Eigen values
		// B * Bt = eye
		// pseudo inverse
	}
	eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
}

/*
 * prepare buffers and call symbolic factorization that is independent on the Kplus values
 */
template <typename T>
static void _set(TotalFETI<T> *dual)
{
	dual->Kplus.resize(dual->feti->K->domains.size());
	dual->d.resize();
	dual->Btx.resize(dual->feti->K->domains.size());
	dual->KplusBtx.resize(dual->feti->K->domains.size());

	#pragma omp parallel for
	for (size_t d = 0; d < dual->feti->K->domains.size(); ++d) {
		dual->Kplus[d].type = dual->feti->K->domains[d].type;
		dual->Kplus[d].shape = dual->feti->K->domains[d].shape;
		math::combine(dual->Kplus[d], dual->feti->K->domains[d], dual->feti->regularization->RegMat.domains[d]);
		dual->Btx[d].resize(dual->feti->K->domains[d].nrows);
		dual->KplusBtx[d].resize(dual->feti->K->domains[d].nrows);
		math::set(dual->Btx[d], T{0});
	}
	eslog::checkpointln("FETI: SET TOTAL-FETI OPERATOR");

	#pragma omp parallel for
	for (size_t d = 0; d < dual->feti->K->domains.size(); ++d) {
		math::commit(dual->Kplus[d]);
		math::symbolicFactorization(dual->Kplus[d]);
	}
	eslog::checkpointln("FETI: TFETI SYMBOLIC FACTORIZATION");
}

template <typename T>
static void _free(TotalFETI<T> *dual)
{
	#pragma omp parallel for
	for (size_t d = 0; d < dual->feti->K->domains.size(); ++d) {
		math::free(dual->Kplus[d]);
	}
}

template <typename T>
static void _update(TotalFETI<T> *dual)
{
	#pragma omp parallel for
	for (size_t d = 0; d < dual->feti->K->domains.size(); ++d) {
		math::sumCombined(dual->Kplus[d], T{1}, dual->feti->K->domains[d], dual->feti->regularization->RegMat.domains[d]);
	}
	eslog::checkpointln("FETI: UPDATE TOTAL-FETI OPERATOR");
	#pragma omp parallel for
	for (size_t d = 0; d < dual->feti->K->domains.size(); ++d) {
		math::numericalFactorization(dual->Kplus[d]);
	}
	eslog::checkpointln("FETI: TFETI NUMERICAL FACTORIZATION");

	#pragma omp parallel for
	for (size_t d = 0; d < dual->feti->K->domains.size(); ++d) {
		math::solve(dual->Kplus[d], dual->feti->f->domains[d], dual->KplusBtx[d]);
	}
	_applyB(dual, dual->KplusBtx, dual->d);
	dual->d.synchronize();
	dual->d.add(T{-1}, dual->feti->equalityConstraints->c);
	eslog::checkpointln("FETI: COMPUTE DUAL RHS [d]");

	_print(dual);
}

template <typename T>
static void _apply(TotalFETI<T> *dual, const Vector_Dual<T> &x, Vector_Dual<T> &y)
{
	#pragma omp parallel for
	for (size_t d = 0; d < dual->feti->K->domains.size(); ++d) {
		_applyBt(dual, d, x, dual->Btx[d]);
		math::solve(dual->Kplus[d], dual->Btx[d], dual->KplusBtx[d]);
	}
	_applyB(dual, dual->KplusBtx, y);
}

template <typename T>
static void _applyBt(TotalFETI<T> *dual, size_t d, const Vector_Dual<T> &in, Vector_Dense<T> &out)
{
	const typename AX_FETI<T>::EqualityConstraints::Domain &L = dual->feti->equalityConstraints->domain[d];

	// compare performance of the loop and set function
	//	for (esint i = 0; i < L->B1.domains[d].nnz; ++i) {
	//		out.vals[L->B1.domains[d].cols[i]] = 0;
	//	}
	math::set(out, T{0});

	// check performance with loop unrolling?
	for (esint r = 0; r < L.B1.nrows; ++r) {
		for (esint c = L.B1.rows[r]; c < L.B1.rows[r + 1]; ++c) {
			out.vals[L.B1.cols[c]] += L.B1.vals[c] * in.vals[L.D2C[r]];
		}
	}
}

// TODO: threaded implementation + more efficient 'beta' scale
template <typename T>
static void _applyB(TotalFETI<T> *dual, const std::vector<Vector_Dense<T> > &in, Vector_Dual<T> &out)
{
	const typename AX_FETI<T>::EqualityConstraints *L = dual->feti->equalityConstraints;

	math::set(out, T{0});
	for (size_t d = 0; d < L->domain.size(); ++d) {
		for (esint r = 0; r < L->domain[d].B1.nrows; ++r) {
			for (esint c = L->domain[d].B1.rows[r]; c < L->domain[d].B1.rows[r + 1]; ++c) {
				out.vals[L->domain[d].D2C[r]] += L->domain[d].B1.vals[c] * in[d].vals[L->domain[d].B1.cols[c]];
			}
		}
	}
	out.synchronize();
}

template <typename T>
static void _print(TotalFETI<T> *dual)
{
	if (info::ecf->output.print_matrices) {
		eslog::storedata(" STORE: feti/dualop/{Kplus}\n");
		for (size_t d = 0; d < dual->feti->K->domains.size(); ++d) {
			math::store(dual->Kplus[d], utils::filename(utils::debugDirectory(*dual->feti->step) + "/feti/dualop", (std::string("Kplus") + std::to_string(d)).c_str()).c_str());
		}
		math::store(dual->d, utils::filename(utils::debugDirectory(*dual->feti->step) + "/feti/dualop", "d").c_str());
	}
}

template <> TotalFETI<double>::TotalFETI(AX_FETI<double> *feti): DualOperator(feti) { _set<double>(this); }
template <> TotalFETI<std::complex<double> >::TotalFETI(AX_FETI<std::complex<double> > *feti): DualOperator(feti) { _set<std::complex<double> >(this); }

template <> TotalFETI<double>::~TotalFETI() { _free<double>(this); }
template <> TotalFETI<std::complex<double> >::~TotalFETI() { _free<std::complex<double> >(this); }

template <> void TotalFETI<double>::info() { _info<double>(this); }
template <> void TotalFETI<std::complex<double> >::info() { _info<std::complex<double> >(this); }

template <> void TotalFETI<double>::update() { _update<double>(this); }
template <> void TotalFETI<std::complex<double> >::update() { _update<std::complex<double> >(this); }

template <> void TotalFETI<double>::apply(const Vector_Dual<double> &x, Vector_Dual<double> &y) { _apply(this, x, y); }
template <> void TotalFETI<std::complex<double> >::apply(const Vector_Dual<std::complex<double> > &x, Vector_Dual<std::complex<double> > &y) { _apply(this, x, y); }

}
