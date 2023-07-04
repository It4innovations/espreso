
#include "totalfeti.implicit.h"
#include "feti/common/applyB.h"

#include "basis/utilities/sysutils.h"
#include "esinfo/eslog.hpp"
#include "esinfo/ecfinfo.h"

namespace espreso {

template class TotalFETIImplicit<double>;
template class TotalFETIImplicit<std::complex<double> >;

template <typename T>
TotalFETIImplicit<T>::TotalFETIImplicit(FETI<T> &feti)
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
	const typename FETI<T>::EqualityConstraints &L = feti.equalityConstraints;

	sum.nnzA = sum.nnzL = sum.memoryL = sum.rows = sum.dualA = sum.surfaceA = 0;
	min.nnzA = min.nnzL = min.memoryL = min.rows = min.dualA = min.surfaceA = INT32_MAX;
	max.nnzA = max.nnzL = max.memoryL = max.rows = max.dualA = max.surfaceA = 0;
	for (size_t d = 0; d < feti.K.domains.size(); ++d) {
		size_t dualA = L.domain[d].B1.nrows;
		size_t surfaceA = Kplus[d].nrows - *std::min_element(L.domain[d].B1.cols, L.domain[d].B1.cols + L.domain[d].B1.nnz);
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
	eslog::info(" =   DOMAINS TOTAL                                                                 %9d = \n", feti.sinfo.domains);
	eslog::info(" =   DUAL SIZE                                                                     %9d = \n", feti.sinfo.lambdasTotal);
	eslog::info(" =   B1 ROWS                                                  %8.0f <%8d - %8d> = \n", (double)sum.dualA / feti.sinfo.domains, min.dualA, max.dualA);
	eslog::info(" =   K+ SURFACE                                               %8.0f <%8d - %8d> = \n", (double)sum.surfaceA / feti.sinfo.domains, min.surfaceA, max.surfaceA);
	eslog::info(" =   K+ ROWS                                                  %8.0f <%8d - %8d> = \n", (double)sum.rows / feti.sinfo.domains, min.rows, max.rows);
	eslog::info(" =   K+ NNZ                                                   %8.0f <%8d - %8d> = \n", (double)sum.nnzA / feti.sinfo.domains, min.nnzA, max.nnzA);
	eslog::info(" =   K+ FACTORS NNZ                                           %8.0f <%8d - %8d> = \n", (double)sum.nnzL / feti.sinfo.domains, min.nnzL, max.nnzL);
	eslog::info(" =   K+ SOLVER MEMORY [MB]                                    %8.2f <%8.2f - %8.2f> = \n", (double)sum.memoryL / feti.sinfo.domains / 1024. / 1024., min.memoryL / 1024. / 1024., max.memoryL / 1024. / 1024.);
	if (sparsity != DirectSolver<T, Matrix_CSR>::VectorSparsity::DENSE) {
		eslog::info(" =   K+ FACTORIZATION                                                        RESPECT SURFACE = \n");
	}
	if (feti.configuration.exhaustive_info) {
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
void TotalFETIImplicit<T>::set(const step::Step &step)
{
	const typename FETI<T>::EqualityConstraints &L = feti.equalityConstraints;
	sparsity = feti.configuration.partial_dual ? DirectSolver<T, Matrix_CSR>::VectorSparsity::SPARSE_RHS | DirectSolver<T, Matrix_CSR>::VectorSparsity::SPARSE_SOLUTION : DirectSolver<T, Matrix_CSR>::VectorSparsity::DENSE;

	Kplus.resize(feti.K.domains.size());
	d.resize();
	Btx.resize(feti.K.domains.size());
	KplusBtx.resize(feti.K.domains.size());
	KSolver.resize(feti.K.domains.size());

	#pragma omp parallel for
	for (size_t d = 0; d < feti.K.domains.size(); ++d) {
		Kplus[d].type = feti.K.domains[d].type;
		Kplus[d].shape = feti.K.domains[d].shape;
		math::combine(Kplus[d], feti.K.domains[d], feti.regularization.RegMat.domains[d]);
		Btx[d].resize(feti.K.domains[d].nrows);
		KplusBtx[d].resize(feti.K.domains[d].nrows);
		math::set(Btx[d], T{0});
	}
	eslog::checkpointln("FETI: SET TOTAL-FETI OPERATOR");

	#pragma omp parallel for
	for (size_t d = 0; d < feti.K.domains.size(); ++d) {
		KSolver[d].commit(Kplus[d]);

		esint suffix = 0;
		if (sparsity != DirectSolver<T, Matrix_CSR>::VectorSparsity::DENSE) {
			suffix = *std::min_element(L.domain[d].B1.cols, L.domain[d].B1.cols + L.domain[d].B1.nnz);
		}

		KSolver[d].symbolicFactorization(suffix);
	}
	eslog::checkpointln("FETI: TFETI SYMBOLIC FACTORIZATION");
}

template <typename T>
void TotalFETIImplicit<T>::update(const step::Step &step)
{
	#pragma omp parallel for
	for (size_t d = 0; d < feti.K.domains.size(); ++d) {
		math::sumCombined(Kplus[d], T{1}, feti.K.domains[d], feti.regularization.RegMat.domains[d]);
	}
	eslog::checkpointln("FETI: UPDATE TOTAL-FETI OPERATOR");
	#pragma omp parallel for
	for (size_t d = 0; d < feti.K.domains.size(); ++d) {
		KSolver[d].numericalFactorization();
	}
	eslog::checkpointln("FETI: TFETI NUMERICAL FACTORIZATION");

	#pragma omp parallel for
	for (size_t d = 0; d < feti.K.domains.size(); ++d) {
		KSolver[d].solve(feti.f.domains[d], KplusBtx[d], sparsity);
	}
	applyB(feti, KplusBtx, d);
	d.add(T{-1}, feti.equalityConstraints.c);
	eslog::checkpointln("FETI: COMPUTE DUAL RHS [d]");

	print(step);
}


template <typename T>
void TotalFETIImplicit<T>::apply(const Vector_Dual<T> &x, Vector_Dual<T> &y)
{
	#pragma omp parallel for
	for (size_t d = 0; d < feti.K.domains.size(); ++d) {
		applyBt(feti, d, x, Btx[d]);
		KSolver[d].solve(Btx[d], KplusBtx[d], sparsity);
	}
	applyB(feti, KplusBtx, y);
}

template <typename T>
void TotalFETIImplicit<T>::toPrimal(const Vector_Dual<T> &x, Vector_FETI<Vector_Dense, T> &y)
{
	#pragma omp parallel for
	for (size_t d = 0; d < feti.K.domains.size(); ++d) {
		applyBt(feti, d, x, Btx[d]);
		math::copy(KplusBtx[d], feti.f.domains[d]);
		math::add(KplusBtx[d], T{-1}, Btx[d]);
		KSolver[d].solve(KplusBtx[d], y.domains[d]);
	}
}

template <typename T>
void TotalFETIImplicit<T>::print(const step::Step &step)
{
	if (info::ecf->output.print_matrices) {
		eslog::storedata(" STORE: feti/dualop/{Kplus}\n");
		for (size_t d = 0; d < feti.K.domains.size(); ++d) {
			math::store(Kplus[d], utils::filename(utils::debugDirectory(step) + "/feti/dualop", (std::string("Kplus") + std::to_string(d)).c_str()).c_str());
		}
		math::store(d, utils::filename(utils::debugDirectory(step) + "/feti/dualop", "d").c_str());
	}
}

}
