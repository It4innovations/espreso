
#include "totalfeti.h"
#include "feti/common/applyB.h"

#include "basis/utilities/sysutils.h"
#include "esinfo/eslog.hpp"
#include "esinfo/ecfinfo.h"

namespace espreso {

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
	const typename FETI<T>::EqualityConstraints *L = dual->feti->equalityConstraints;
	dual->sparsity = dual->feti->configuration.restricted_dual ? math::VectorSparsity::SPARSE_RHS | math::VectorSparsity::SPARSE_SOLUTION : math::VectorSparsity::DENSE;

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
		math::initSolver(dual->Kplus[d]);

		esint suffix = 0;
		if (dual->sparsity != math::VectorSparsity::DENSE) {
			suffix = *std::min_element(L->domain[d].B1.cols, L->domain[d].B1.cols + L->domain[d].B1.nnz);
		}

		math::symbolicFactorization(dual->Kplus[d], suffix);
	}
	eslog::checkpointln("FETI: TFETI SYMBOLIC FACTORIZATION");
}

template <typename T>
static void _free(TotalFETI<T> *dual)
{
	#pragma omp parallel for
	for (size_t d = 0; d < dual->feti->K->domains.size(); ++d) {
		math::freeSolver(dual->Kplus[d]);
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
		math::solve(dual->Kplus[d], dual->feti->f->domains[d], dual->KplusBtx[d], math::VectorSparsity::SPARSE_RHS);
	}
	applyB(dual->feti, dual->KplusBtx, dual->d);
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
		applyBt(dual->feti, d, x, dual->Btx[d]);
		math::solve(dual->Kplus[d], dual->Btx[d], dual->KplusBtx[d], dual->sparsity);
	}
	applyB(dual->feti, dual->KplusBtx, y);
}

template <typename T>
static void _toPrimal(TotalFETI<T> *dual, const Vector_Dual<T> &x, Vector_FETI<Vector_Dense, T> &y)
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

template <> TotalFETI<double>::TotalFETI(FETI<double> *feti): DualOperator(feti) { _set<double>(this); }
template <> TotalFETI<std::complex<double> >::TotalFETI(FETI<std::complex<double> > *feti): DualOperator(feti) { _set<std::complex<double> >(this); }

template <> TotalFETI<double>::~TotalFETI() { _free<double>(this); }
template <> TotalFETI<std::complex<double> >::~TotalFETI() { _free<std::complex<double> >(this); }

template <> void TotalFETI<double>::info() { _info<double>(this); }
template <> void TotalFETI<std::complex<double> >::info() { _info<std::complex<double> >(this); }

template <> void TotalFETI<double>::update() { _update<double>(this); }
template <> void TotalFETI<std::complex<double> >::update() { _update<std::complex<double> >(this); }

template <> void TotalFETI<double>::apply(const Vector_Dual<double> &x, Vector_Dual<double> &y) { _apply(this, x, y); }
template <> void TotalFETI<std::complex<double> >::apply(const Vector_Dual<std::complex<double> > &x, Vector_Dual<std::complex<double> > &y) { _apply(this, x, y); }

template <> void TotalFETI<double>::toPrimal(const Vector_Dual<double> &x, Vector_FETI<Vector_Dense, double> &y) { _toPrimal(this, x, y); }
template <> void TotalFETI<std::complex<double> >::toPrimal(const Vector_Dual<std::complex<double> > &x, Vector_FETI<Vector_Dense, std::complex<double> > &y) { _toPrimal(this, x, y); }

}
