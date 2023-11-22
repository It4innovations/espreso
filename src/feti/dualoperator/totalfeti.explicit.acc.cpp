
#include "totalfeti.explicit.acc.h"
#include "feti/common/applyB.h"

#include "basis/utilities/sysutils.h"
#include "esinfo/eslog.hpp"
#include "esinfo/ecfinfo.h"

namespace espreso {

template <typename T>
TotalFETIExplicitAcc<T>::TotalFETIExplicitAcc(FETI<T> &feti)
: TotalFETIExplicit<T>(feti)
{

}

template <typename T>
TotalFETIExplicitAcc<T>::~TotalFETIExplicitAcc()
{

}

template <typename T>
void TotalFETIExplicitAcc<T>::info()
{
	DualOperatorInfo sum, min, max;
	TotalFETIImplicit<T>::reduceInfo(sum, min, max);

	eslog::info(" = ACCELERATED TOTAL FETI OPERATOR                                                           = \n");
	TotalFETIImplicit<T>::printInfo(sum, min, max);
	eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
}

template <typename T>
void TotalFETIExplicitAcc<T>::set(const step::Step &step)
{
	TotalFETIExplicit<T>::set(step);

	L.resize(feti.K.domains.size());
	U.resize(feti.K.domains.size());
	p.resize(feti.K.domains.size());
}

template <typename T>
void TotalFETIExplicitAcc<T>::update(const step::Step &step)
{
	TotalFETIExplicit<T>::update(step);

	#pragma omp parallel for
	for (size_t d = 0; d < feti.K.domains.size(); ++d) {
		KSolver[d].getFactors(L[d], U[d], p[d]);
		switch (Kplus[d].type) {
		case Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE:
		case Matrix_Type::COMPLEX_HERMITIAN_POSITIVE_DEFINITE:
			// TODO: copy factors to ACC, solve, ...
			// results should be equal to F[d]
			break;
		default:
			// TODO: implement non-symmetric case
			break;
		}
	}
}

template <typename T>
void TotalFETIExplicitAcc<T>::apply(const Vector_Dual<T> &x, Vector_Dual<T> &y)
{
	TotalFETIExplicit<T>::apply(x, y);
}

template <typename T>
void TotalFETIExplicitAcc<T>::toPrimal(const Vector_Dual<T> &x, Vector_FETI<Vector_Dense, T> &y)
{
	TotalFETIExplicit<T>::toPrimal(x, y);
}

template class TotalFETIExplicitAcc<double>;
template class TotalFETIExplicitAcc<std::complex<double> >;

}
