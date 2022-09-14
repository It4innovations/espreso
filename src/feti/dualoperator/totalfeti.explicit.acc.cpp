
#include "totalfeti.explicit.acc.h"
#include "feti/common/applyB.h"

#include "basis/utilities/sysutils.h"
#include "esinfo/eslog.hpp"
#include "esinfo/ecfinfo.h"

namespace espreso {

template class TotalFETIExplicitAcc<double>;
template class TotalFETIExplicitAcc<std::complex<double> >;

template <typename T>
TotalFETIExplicitAcc<T>::TotalFETIExplicitAcc(FETI<T> *feti)
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
	math::SolverInfo sum, min, max;
	TotalFETIImplicit<T>::reduceInfo(sum, min, max);

	eslog::info(" = ACCELERATED TOTAL FETI OPERATOR                                                           = \n");
	TotalFETIImplicit<T>::printInfo(sum, min, max);
	eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
}

template <typename T>
void TotalFETIExplicitAcc<T>::set()
{
	TotalFETIExplicit<T>::set();

	L.resize(this->feti->K->domains.size());
	U.resize(this->feti->K->domains.size());
	p.resize(this->feti->K->domains.size());
}

template <typename T>
void TotalFETIExplicitAcc<T>::update()
{
	TotalFETIExplicit<T>::update();

	const typename FETI<T>::EqualityConstraints *L = this->feti->equalityConstraints;

	#pragma omp parallel for
	for (size_t d = 0; d < this->feti->K->domains.size(); ++d) {
		math::getFactors(this->Kplus[d], this->L[d], this->U[d], this->p[d]);
		switch (this->Kplus[d].type) {
		case Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE:
		case Matrix_Type::COMPLEX_HERMITIAN_POSITIVE_DEFINITE:
			// TODO: copy factors to ACC, solve, ...
			// results should be equal to this->F[d]
			break;
		default:
			// TODO: implement non-symmetric case
			break;
		}

		math::freeFactor(this->L[d]); // ?
//		math::freeFactor(this->U[d]);
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

}
