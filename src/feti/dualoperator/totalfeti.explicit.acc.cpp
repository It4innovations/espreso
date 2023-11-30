
#include "totalfeti.explicit.acc.h"
#include "feti/common/applyB.h"

#include "basis/utilities/sysutils.h"
#include "esinfo/eslog.hpp"
#include "esinfo/ecfinfo.h"
#include "esinfo/mpiinfo.h"

namespace espreso {

template <typename T>
TotalFETIExplicitAcc<T>::TotalFETIExplicitAcc(FETI<T> &feti)
: DualOperator<T>(feti), acc(info::mpi::rank)
{

}

template <typename T>
TotalFETIExplicitAcc<T>::~TotalFETIExplicitAcc()
{

}

template <typename T>
void TotalFETIExplicitAcc<T>::info()
{
//	DualOperatorInfo sum, min, max;
//	size_t minF = INT32_MAX, maxF = 0, sumF = 0;
//	for (size_t d = 0; d < F.size(); ++d) {
//		minF = std::min(minF, F[d].nrows * F[d].ncols * sizeof(double));
//		maxF = std::max(maxF, F[d].nrows * F[d].ncols * sizeof(double));
//		sumF += F[d].nrows * F[d].ncols * sizeof(double);
//	}
//
////	TotalFETIImplicit<T>::reduceInfo(sum, min, max);
//	Communication::allReduce(&minF, nullptr, 1, MPITools::getType<size_t>().mpitype, MPI_MIN);
//	Communication::allReduce(&maxF, nullptr, 1, MPITools::getType<size_t>().mpitype, MPI_MAX);
//	Communication::allReduce(&sumF, nullptr, 1, MPITools::getType<size_t>().mpitype, MPI_SUM);
//
//	eslog::info(" = EXPLICIT TOTAL FETI OPERATOR ON GPU                                                       = \n");
////	TotalFETIImplicit<T>::printInfo(sum, min, max);
//	eslog::info(" =   F MEMORY [MB]                                            %8.2f <%8.2f - %8.2f> = \n", (double)sumF / F.size() / 1024. / 1024., minF / 1024. / 1024., maxF / 1024. / 1024.);
//	eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
}

template <typename T>
void TotalFETIExplicitAcc<T>::set(const step::Step &step)
{
	const typename FETI<T>::EqualityConstraints &L = feti.equalityConstraints;

	Kplus.resize(feti.K.domains.size());
	Btx.resize(feti.K.domains.size());
	KplusBtx.resize(feti.K.domains.size());
	KSolver.resize(feti.K.domains.size());
	B1.resize(feti.K.domains.size());
	D2C.resize(feti.K.domains.size());
	d.resize();

	#pragma omp parallel for
	for (size_t di = 0; di < feti.K.domains.size(); ++di) {
		Kplus[di].type = feti.K.domains[di].type;
		Kplus[di].shape = feti.K.domains[di].shape;
		math::combine(Kplus[di], feti.K.domains[di], feti.regularization.RegMat.domains[di]);
		Btx[di].resize(feti.K.domains[di].nrows);
		KplusBtx[di].resize(feti.K.domains[di].nrows);
		math::set(Btx[di], T{0});
		B1[di].shallowCopy(L.domain[di].B1);
		D2C[di] = L.domain[di].D2C;
	}
	eslog::checkpointln("FETI: SET TOTAL-FETI OPERATOR");

	#pragma omp parallel for
	for (size_t di = 0; di < feti.K.domains.size(); ++di) {
		KSolver[di].commit(Kplus[di]);
		KSolver[di].symbolicFactorization(0);
	}
	eslog::checkpointln("FETI: TFETI SYMBOLIC FACTORIZATION");

	acc.set(Kplus, B1);
	eslog::checkpointln("FETI: TFETI SET ACC");
}

template <typename T>
void TotalFETIExplicitAcc<T>::update(const step::Step &step)
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

	acc.update(Kplus);
	eslog::checkpointln("FETI: TFETI UPDATE ACC");

	#pragma omp parallel for
	for (size_t di = 0; di < feti.K.domains.size(); ++di) {
		KSolver[di].solve(feti.f.domains[di], KplusBtx[di]);
	}
	applyB(feti, KplusBtx, d);
	d.add(T{-1}, feti.equalityConstraints.c);
	eslog::checkpointln("FETI: COMPUTE DUAL RHS [d]");

//	print(step);
}

template <typename T>
void TotalFETIExplicitAcc<T>::apply(const Vector_Dual<T> &x, Vector_Dual<T> &y)
{
	acc.apply(x, y, D2C);
	y.synchronize();
}

template <typename T>
void TotalFETIExplicitAcc<T>::toPrimal(const Vector_Dual<T> &x, Vector_FETI<Vector_Dense, T> &y)
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
		eslog::storedata(" STORE: feti/dualop/{Kplus}\n");
		for (size_t di = 0; di < feti.K.domains.size(); ++di) {
			math::store(Kplus[di], utils::filename(utils::debugDirectory(step) + "/feti/dualop", (std::string("Kplus") + std::to_string(di)).c_str()).c_str());
		}
		math::store(d, utils::filename(utils::debugDirectory(step) + "/feti/dualop", "d").c_str());
	}
}

template class TotalFETIExplicitAcc<double>;
template class TotalFETIExplicitAcc<std::complex<double> >;

}
