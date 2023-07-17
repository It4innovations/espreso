
#include "lumped.h"
#include "feti/common/applyB.h"

#include "esinfo/eslog.hpp"

namespace espreso {

template struct Lumped<double>;
template struct Lumped<std::complex<double> >;

template <typename T>
Lumped<T>::Lumped(FETI<T> &feti)
: Preconditioner<T>(feti)
{
	Btx.resize(feti.K.domains.size());
	KBtx.resize(feti.K.domains.size());
	KSpBlas.resize(feti.K.domains.size());

	#pragma omp parallel for
	for (size_t d = 0; d < feti.K.domains.size(); ++d) {
		KSpBlas[d].insert(feti.K.domains[d]);
		Btx[d].resize(feti.K.domains[d].nrows);
		KBtx[d].resize(feti.K.domains[d].nrows);
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
	for (size_t d = 0; d < feti.K.domains.size(); ++d) {
		applyBt(feti, d, x, Btx[d]);
		KSpBlas[d].apply(KBtx[d], T{1}, T{0}, Btx[d]);
	}
	applyB(feti, KBtx, y);
}

}
