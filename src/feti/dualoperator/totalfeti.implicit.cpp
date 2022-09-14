
#include "totalfeti.implicit.h"
#include "feti/common/applyB.h"

#include "basis/utilities/sysutils.h"
#include "esinfo/eslog.hpp"
#include "esinfo/ecfinfo.h"

namespace espreso {

template class TotalFETIImplicit<double>;
template class TotalFETIImplicit<std::complex<double> >;

template <typename T>
TotalFETIImplicit<T>::TotalFETIImplicit(FETI<T> *feti)
: DualOperator<T>(feti), sparsity(math::VectorSparsity::DENSE)
{

}

template <typename T>
TotalFETIImplicit<T>::~TotalFETIImplicit()
{
	#pragma omp parallel for
	for (size_t d = 0; d < this->feti->K->domains.size(); ++d) {
		math::freeSolver(this->Kplus[d]);
	}
}

template <typename T>
void TotalFETIImplicit<T>::reduceInfo(math::SolverInfo &sum, math::SolverInfo &min, math::SolverInfo &max)
{
	sum.nnzA = sum.nnzL = sum.rows = 0;
	min.nnzA = min.nnzL = min.rows = INT32_MAX;
	max.nnzA = max.nnzL = max.rows = 0;
	for (size_t d = 0; d < this->feti->K->domains.size(); ++d) {
		auto info = math::getSolverInfo(this->Kplus[d]);
		min.rows = std::min(min.rows, info.rows);
		min.nnzA = std::min(min.nnzA, info.nnzA);
		min.nnzL = std::min(min.nnzL, info.nnzL);
		max.rows = std::max(max.rows, info.rows);
		max.nnzA = std::max(max.nnzA, info.nnzA);
		max.nnzL = std::max(max.nnzL, info.nnzL);
		sum.rows += info.rows;
		sum.nnzA += info.nnzA;
		sum.nnzL += info.nnzL;
	}

	Communication::allReduce(&min.rows, nullptr, 3, MPI_INT, MPI_MIN);
	Communication::allReduce(&max.rows, nullptr, 3, MPI_INT, MPI_MAX);
	Communication::allReduce(&sum.rows, nullptr, 3, MPI_INT, MPI_SUM);
}

template <typename T>
void TotalFETIImplicit<T>::printInfo(math::SolverInfo &sum, math::SolverInfo &min, math::SolverInfo &max)
{
	eslog::info(" =   DOMAINS                                                                       %9d = \n", this->feti->sinfo.domains);
	eslog::info(" =   DUAL SIZE                                                                     %9d = \n", this->feti->sinfo.lambdasTotal);
	eslog::info(" =   K+ ROWS                                                  %8.0f <%8d - %8d> = \n", (double)sum.rows / this->feti->sinfo.domains, min.rows, max.rows);
	eslog::info(" =   K+ NNZ                                                   %8.0f <%8d - %8d> = \n", (double)sum.nnzA / this->feti->sinfo.domains, min.nnzA, max.nnzA);
	eslog::info(" =   K+ FACTORS NNZ                                           %8.0f <%8d - %8d> = \n", (double)sum.nnzL / this->feti->sinfo.domains, min.nnzL, max.nnzL);
	if (this->sparsity != math::VectorSparsity::DENSE) {
		eslog::info(" =   K+ FACTORIZATION                                                        RESPECT SURFACE = \n");
	}
	if (this->feti->configuration.exhaustive_info) {
		// power method to Eigen values
		// B * Bt = eye
		// pseudo inverse
	}
}

template <typename T>
void TotalFETIImplicit<T>::info()
{
	math::SolverInfo sum, min, max;
	reduceInfo(sum, min, max);

	eslog::info(" = IMPLICIT TOTAL FETI OPERATOR                                                              = \n");
	printInfo(sum, min, max);
	eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
}

/*
 * prepare buffers and call symbolic factorization that is independent on the Kplus values
 */
template <typename T>
void TotalFETIImplicit<T>::set()
{
	const typename FETI<T>::EqualityConstraints *L = this->feti->equalityConstraints;
	this->sparsity = this->feti->configuration.partial_dual ? math::VectorSparsity::SPARSE_RHS | math::VectorSparsity::SPARSE_SOLUTION : math::VectorSparsity::DENSE;

	this->Kplus.resize(this->feti->K->domains.size());
	this->d.resize();
	this->Btx.resize(this->feti->K->domains.size());
	this->KplusBtx.resize(this->feti->K->domains.size());

	#pragma omp parallel for
	for (size_t d = 0; d < this->feti->K->domains.size(); ++d) {
		this->Kplus[d].type = this->feti->K->domains[d].type;
		this->Kplus[d].shape = this->feti->K->domains[d].shape;
		math::combine(this->Kplus[d], this->feti->K->domains[d], this->feti->regularization->RegMat.domains[d]);
		this->Btx[d].resize(this->feti->K->domains[d].nrows);
		this->KplusBtx[d].resize(this->feti->K->domains[d].nrows);
		math::set(this->Btx[d], T{0});
	}
	eslog::checkpointln("FETI: SET TOTAL-FETI OPERATOR");

	#pragma omp parallel for
	for (size_t d = 0; d < this->feti->K->domains.size(); ++d) {
		math::initSolver(this->Kplus[d]);

		esint suffix = 0;
		if (this->sparsity != math::VectorSparsity::DENSE) {
			suffix = *std::min_element(L->domain[d].B1.cols, L->domain[d].B1.cols + L->domain[d].B1.nnz);
		}

		math::symbolicFactorization(this->Kplus[d], suffix);
	}
	eslog::checkpointln("FETI: TFETI SYMBOLIC FACTORIZATION");
}

template <typename T>
void TotalFETIImplicit<T>::update()
{
	#pragma omp parallel for
	for (size_t d = 0; d < this->feti->K->domains.size(); ++d) {
		math::sumCombined(this->Kplus[d], T{1}, this->feti->K->domains[d], this->feti->regularization->RegMat.domains[d]);
	}
	eslog::checkpointln("FETI: UPDATE TOTAL-FETI OPERATOR");
	#pragma omp parallel for
	for (size_t d = 0; d < this->feti->K->domains.size(); ++d) {
		math::numericalFactorization(this->Kplus[d]);
	}
	eslog::checkpointln("FETI: TFETI NUMERICAL FACTORIZATION");

	#pragma omp parallel for
	for (size_t d = 0; d < this->feti->K->domains.size(); ++d) {
		math::solve(this->Kplus[d], this->feti->f->domains[d], this->KplusBtx[d], math::VectorSparsity::SPARSE_RHS);
	}
	applyB(this->feti, this->KplusBtx, this->d);
	this->d.synchronize();
	this->d.add(T{-1}, this->feti->equalityConstraints->c);
	eslog::checkpointln("FETI: COMPUTE DUAL RHS [d]");

	printMatrices();
}


template <typename T>
void TotalFETIImplicit<T>::apply(const Vector_Dual<T> &x, Vector_Dual<T> &y)
{
	#pragma omp parallel for
	for (size_t d = 0; d < this->feti->K->domains.size(); ++d) {
		applyBt(this->feti, d, x, this->Btx[d]);
		math::solve(this->Kplus[d], this->Btx[d], this->KplusBtx[d], this->sparsity);
	}
	applyB(this->feti, this->KplusBtx, y);
}

template <typename T>
void TotalFETIImplicit<T>::toPrimal(const Vector_Dual<T> &x, Vector_FETI<Vector_Dense, T> &y)
{
	#pragma omp parallel for
	for (size_t d = 0; d < this->feti->K->domains.size(); ++d) {
		applyBt(this->feti, d, x, this->Btx[d]);
		math::copy(this->KplusBtx[d], this->feti->f->domains[d]);
		math::add(this->KplusBtx[d], T{-1}, this->Btx[d]);
		math::solve(this->Kplus[d], this->KplusBtx[d], y.domains[d]);
	}
}

template <typename T>
void TotalFETIImplicit<T>::printMatrices()
{
	if (info::ecf->output.print_matrices) {
		eslog::storedata(" STORE: feti/dualop/{Kplus}\n");
		for (size_t d = 0; d < this->feti->K->domains.size(); ++d) {
			math::store(this->Kplus[d], utils::filename(utils::debugDirectory(*this->feti->step) + "/feti/dualop", (std::string("Kplus") + std::to_string(d)).c_str()).c_str());
		}
		math::store(this->d, utils::filename(utils::debugDirectory(*this->feti->step) + "/feti/dualop", "d").c_str());
	}
}

}