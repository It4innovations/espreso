
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
: DualOperator<T>(feti), sparsity(DirectSolver<T, Matrix_CSR>::VectorSparsity::DENSE)
{

}

template <typename T>
TotalFETIImplicit<T>::~TotalFETIImplicit()
{

}

template <typename T>
void TotalFETIImplicit<T>::reduceInfo(DualOperatorInfo &sum, DualOperatorInfo &min, DualOperatorInfo &max)
{
	const typename FETI<T>::EqualityConstraints *L = this->feti->equalityConstraints;

	sum.nnzA = sum.nnzL = sum.memoryL = sum.rows = sum.dualA = sum.surfaceA = 0;
	min.nnzA = min.nnzL = min.memoryL = min.rows = min.dualA = min.surfaceA = INT32_MAX;
	max.nnzA = max.nnzL = max.memoryL = max.rows = max.dualA = max.surfaceA = 0;
	for (size_t d = 0; d < this->feti->K->domains.size(); ++d) {
		size_t dualA = L->domain[d].B1.nrows;
		size_t surfaceA = this->Kplus[d].nrows - *std::min_element(L->domain[d].B1.cols, L->domain[d].B1.cols + L->domain[d].B1.nnz);
		min.rows = std::min(min.rows, KSolver[d].rows);
		min.nnzA = std::min(min.nnzA, KSolver[d].nnzA);
		min.nnzL = std::min(min.nnzL, KSolver[d].nnzL);
		min.memoryL = std::min(min.memoryL, KSolver[d].memoryL);
		min.dualA = std::min(min.dualA, dualA);
		min.surfaceA = std::min(min.surfaceA, surfaceA);
		max.rows = std::max(max.rows, KSolver[d].rows);
		max.nnzA = std::max(max.nnzA, KSolver[d].nnzA);
		max.nnzL = std::max(max.nnzL, KSolver[d].nnzL);
		max.memoryL = std::max(max.memoryL, KSolver[d].memoryL);
		max.dualA = std::max(max.dualA, dualA);
		max.surfaceA = std::max(max.surfaceA, surfaceA);
		sum.rows += KSolver[d].rows;
		sum.nnzA += KSolver[d].nnzA;
		sum.nnzL += KSolver[d].nnzL;
		sum.memoryL += KSolver[d].memoryL;
		sum.dualA += dualA;
		sum.surfaceA += surfaceA;
	}

	Communication::allReduce(&min, nullptr, 6, MPITools::getType<size_t>().mpitype, MPI_MIN);
	Communication::allReduce(&max, nullptr, 6, MPITools::getType<size_t>().mpitype, MPI_MAX);
	Communication::allReduce(&sum, nullptr, 6, MPITools::getType<size_t>().mpitype, MPI_SUM);
}

template <typename T>
void TotalFETIImplicit<T>::printInfo(DualOperatorInfo &sum, DualOperatorInfo &min, DualOperatorInfo &max)
{
	eslog::info(" =   DOMAINS TOTAL                                                                 %9d = \n", this->feti->sinfo.domains);
	eslog::info(" =   DUAL SIZE                                                                     %9d = \n", this->feti->sinfo.lambdasTotal);
	eslog::info(" =   B1 ROWS                                                  %8.0f <%8d - %8d> = \n", (double)sum.dualA / this->feti->sinfo.domains, min.dualA, max.dualA);
	eslog::info(" =   K+ SURFACE                                               %8.0f <%8d - %8d> = \n", (double)sum.surfaceA / this->feti->sinfo.domains, min.surfaceA, max.surfaceA);
	eslog::info(" =   K+ ROWS                                                  %8.0f <%8d - %8d> = \n", (double)sum.rows / this->feti->sinfo.domains, min.rows, max.rows);
	eslog::info(" =   K+ NNZ                                                   %8.0f <%8d - %8d> = \n", (double)sum.nnzA / this->feti->sinfo.domains, min.nnzA, max.nnzA);
	eslog::info(" =   K+ FACTORS NNZ                                           %8.0f <%8d - %8d> = \n", (double)sum.nnzL / this->feti->sinfo.domains, min.nnzL, max.nnzL);
	eslog::info(" =   K+ SOLVER MEMORY [MB]                                    %8.2f <%8.2f - %8.2f> = \n", (double)sum.memoryL / this->feti->sinfo.domains / 1024. / 1024., min.memoryL / 1024. / 1024., max.memoryL / 1024. / 1024.);
	if (this->sparsity != DirectSolver<T, Matrix_CSR>::VectorSparsity::DENSE) {
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
	DualOperatorInfo sum, min, max;
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
	this->sparsity = this->feti->configuration.partial_dual ? DirectSolver<T, Matrix_CSR>::VectorSparsity::SPARSE_RHS | DirectSolver<T, Matrix_CSR>::VectorSparsity::SPARSE_SOLUTION : DirectSolver<T, Matrix_CSR>::VectorSparsity::DENSE;

	this->Kplus.resize(this->feti->K->domains.size());
	this->d.resize();
	this->Btx.resize(this->feti->K->domains.size());
	this->KplusBtx.resize(this->feti->K->domains.size());
	this->KSolver.resize(this->feti->K->domains.size());

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
		KSolver[d].commit(this->Kplus[d]);

		esint suffix = 0;
		if (this->sparsity != DirectSolver<T, Matrix_CSR>::VectorSparsity::DENSE) {
			suffix = *std::min_element(L->domain[d].B1.cols, L->domain[d].B1.cols + L->domain[d].B1.nnz);
		}

		KSolver[d].symbolicFactorization(suffix);
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
		KSolver[d].numericalFactorization();
	}
	eslog::checkpointln("FETI: TFETI NUMERICAL FACTORIZATION");

	#pragma omp parallel for
	for (size_t d = 0; d < this->feti->K->domains.size(); ++d) {
		KSolver[d].solve(this->feti->f->domains[d], this->KplusBtx[d], sparsity);
	}
	applyB(this->feti, this->KplusBtx, this->d);
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
		KSolver[d].solve(this->Btx[d], this->KplusBtx[d], sparsity);
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
		KSolver[d].solve(this->KplusBtx[d], y.domains[d]);
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
