
#include "lumped.h"
#include "feti/common/applyB.h"

#include "esinfo/eslog.hpp"

namespace espreso {

template <typename T>
Lumped<T>::Lumped(FETI<T> &feti)
: Preconditioner<T>(feti)
{
	Btx.resize(feti.K.size());
	KBtx.resize(feti.K.size());
	KSpBlas.resize(feti.K.size());

	#pragma omp parallel for
	for (size_t d = 0; d < feti.K.size(); ++d) {
		KSpBlas[d].insert(feti.K[d]);
		Btx[d].resize(feti.K[d].nrows);
		KBtx[d].resize(feti.K[d].nrows);
	}

	eslog::checkpointln("FETI: SET LUMPED PRECONDITIONER");
}

template <typename T>
Lumped<T>::~Lumped()
{

}

template <typename T>
void Lumped<T>::info()
{
	if (feti.configuration.exhaustive_info) {
		eslog::info(" = LUMPED PRECONDITIONER PROPERTIES                                                          = \n");

		eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
	}
}

template <typename T>
void Lumped<T>::update(const step::Step &step)
{

}

template <typename T> void
Lumped<T>::apply(const Vector_Dual<T> &x, Vector_Dual<T> &y)
{
	#pragma omp parallel for
	for (size_t d = 0; d < feti.K.size(); ++d) {
		applyBt(feti, d, x, Btx[d]);
		KSpBlas[d].apply(KBtx[d], T{1}, T{0}, Btx[d]);
	}
	applyB(feti, KBtx, y);
}

template struct Lumped<double>;
template struct Lumped<std::complex<double> >;

}
