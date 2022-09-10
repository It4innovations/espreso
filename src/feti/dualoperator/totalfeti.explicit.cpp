
#include "totalfeti.explicit.h"
#include "feti/common/applyB.h"

#include "basis/utilities/sysutils.h"
#include "esinfo/eslog.hpp"
#include "esinfo/ecfinfo.h"
#include "math/physics/math.physics.copy.hpp"

namespace espreso {

template <typename T> static void _print(TotalFETIExplicit<T> *dual);

template <typename T>
static void _info(TotalFETIExplicit<T> *dual)
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
static void _set(TotalFETIExplicit<T> *dual)
{
	dual->F.resize(dual->feti->K->domains.size());

	const typename FETI<T>::EqualityConstraints *L = dual->feti->equalityConstraints;
	#pragma omp parallel for
	for (size_t d = 0; d < dual->feti->K->domains.size(); ++d) {
		dual->F[d].resize(L->domain[d].B1.nrows, L->domain[d].B1.nrows);
	}
	eslog::checkpointln("FETI: TFETI SET EXPLICIT OPERATOR");
}

template <typename T>
static void _free(TotalFETIExplicit<T> *dual)
{

}

template <typename T>
static void _update(TotalFETIExplicit<T> *dual)
{
	const typename FETI<T>::EqualityConstraints *L = dual->feti->equalityConstraints;

	#pragma omp parallel for
	for (size_t d = 0; d < dual->feti->K->domains.size(); ++d) {
		math::solve(dual->Kplus[d], dual->feti->f->domains[d], dual->KplusBtx[d], math::VectorSparsity::SPARSE_RHS);
	}

	
	// #pragma omp parallel for
	for (size_t d = 0; d < dual->feti->K->domains.size(); ++d) {
		double start = eslog::time(), _start;

		_start = eslog::time();
		Matrix_Dense<T> KplusBt, Bt;
		KplusBt.resize(L->domain[d].B1.nrows, L->domain[d].B1.ncols);
		Bt.resize(L->domain[d].B1.nrows, L->domain[d].B1.ncols);
		printf("resize: %fs\n", eslog::time() - _start);

		_start = eslog::time();
		math::copy(Bt, L->domain[d].B1);
		printf("copy: %fs\n", eslog::time() - _start);

		_start = eslog::time();
		math::solve(dual->Kplus[d], Bt, KplusBt, dual->sparsity);
		printf("solve: %fs\n", eslog::time() - _start);

		_start = eslog::time();
		math::set(dual->F[d], T{0});
		printf("setF: %fs\n", eslog::time() - _start);

		_start = eslog::time();
		for (esint r = 0; r < L->domain[d].B1.nrows; ++r) {
			for (esint c = 0; c < L->domain[d].B1.nrows; ++c) {
				for (esint k = L->domain[d].B1.rows[r]; k < L->domain[d].B1.rows[r + 1]; ++k) {
					// printf("%d / %d::%d / %d::%d / %d\n", r * L->domain[d].B1.nrows + c, dual->F[d].nnz, k, L->domain[d].B1.nnz, c * L->domain[d].B1.ncols + L->domain[d].B1.cols[k], KplusBt.nnz);
					dual->F[d].vals[r * L->domain[d].B1.nrows + c] += L->domain[d].B1.vals[k] * KplusBt.vals[c * L->domain[d].B1.ncols + L->domain[d].B1.cols[k]];
				}
			}
		}
		printf("mult: %fs\n", eslog::time() - _start);
		printf("explicit B * K+ * Bt: %fs\n", eslog::time() - start);
	}
	

	eslog::checkpointln("FETI: TFETI EXPLICIT DUAL ASSEMBLED");
	_print(dual);
}

template <typename T>
static void _apply(TotalFETIExplicit<T> *dual, const Vector_Dual<T> &x, Vector_Dual<T> &y)
{

}

template <typename T>
static void _toPrimal(TotalFETIExplicit<T> *dual, const Vector_Dual<T> &x, Vector_FETI<Vector_Dense, T> &y)
{
	#pragma omp parallel for
	for (size_t d = 0; d < dual->feti->K->domains.size(); ++d) {
		applyBt(dual->feti, d, x, dual->Btx[d]);
		math::copy(dual->KplusBtx[d], dual->feti->f->domains[d]);
		math::add(dual->KplusBtx[d], T{-1}, dual->Btx[d]);
		math::solve(dual->Kplus[d], dual->KplusBtx[d], y.domains[d]);
	}
}

template <typename T>
static void _print(TotalFETIExplicit<T> *dual)
{
	if (info::ecf->output.print_matrices) {
		eslog::storedata(" STORE: feti/dualop/{F}\n");
		for (size_t d = 0; d < dual->feti->K->domains.size(); ++d) {
			math::store(dual->Kplus[d], utils::filename(utils::debugDirectory(*dual->feti->step) + "/feti/dualop", (std::string("F") + std::to_string(d)).c_str()).c_str());
		}
	}
}

template <> TotalFETIExplicit<double>::TotalFETIExplicit(FETI<double> *feti): TotalFETIImplicit(feti) { _set<double>(this); }
template <> TotalFETIExplicit<std::complex<double> >::TotalFETIExplicit(FETI<std::complex<double> > *feti): TotalFETIImplicit(feti) { _set<std::complex<double> >(this); }

template <> TotalFETIExplicit<double>::~TotalFETIExplicit() { _free<double>(this); }
template <> TotalFETIExplicit<std::complex<double> >::~TotalFETIExplicit() { _free<std::complex<double> >(this); }

template <> void TotalFETIExplicit<double>::info() { _info<double>(this); }
template <> void TotalFETIExplicit<std::complex<double> >::info() { _info<std::complex<double> >(this); }

template <> void TotalFETIExplicit<double>::update() { TotalFETIImplicit::update(); _update<double>(this); }
template <> void TotalFETIExplicit<std::complex<double> >::update() { TotalFETIImplicit::update(); _update<std::complex<double> >(this); }

template <> void TotalFETIExplicit<double>::apply(const Vector_Dual<double> &x, Vector_Dual<double> &y) { TotalFETIImplicit::apply(x, y); }
template <> void TotalFETIExplicit<std::complex<double> >::apply(const Vector_Dual<std::complex<double> > &x, Vector_Dual<std::complex<double> > &y) { TotalFETIImplicit::apply(x, y); }

template <> void TotalFETIExplicit<double>::toPrimal(const Vector_Dual<double> &x, Vector_FETI<Vector_Dense, double> &y) { _toPrimal(this, x, y); }
template <> void TotalFETIExplicit<std::complex<double> >::toPrimal(const Vector_Dual<std::complex<double> > &x, Vector_FETI<Vector_Dense, std::complex<double> > &y) { _toPrimal(this, x, y); }

}
