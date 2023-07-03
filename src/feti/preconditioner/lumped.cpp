
#include "lumped.h"
#include "feti/common/applyB.h"

#include "esinfo/eslog.hpp"

namespace espreso {

template struct Lumped<double>;
template struct Lumped<std::complex<double> >;

template <typename T>
Lumped<T>::Lumped(FETI<T> *feti)
: Preconditioner<T>(feti), K(feti->K)
{
	Btx.resize(K->domains.size());
	KBtx.resize(K->domains.size());
	KSpBlas.resize(K->domains.size());

	#pragma omp parallel for
	for (size_t d = 0; d < K->domains.size(); ++d) {
		KSpBlas[d].commit(K->domains[d]);
		Btx[d].resize(K->domains[d].nrows);
		KBtx[d].resize(K->domains[d].nrows);
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
	if (this->feti->configuration.exhaustive_info) {
		eslog::info(" = LUMPED PRECONDITIONER PROPERTIES                                                          = \n");

		eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
	}
}

template <typename T>
void Lumped<T>::update()
{

}

template <typename T> void
Lumped<T>::apply(const Vector_Dual<T> &x, Vector_Dual<T> &y)
{
	#pragma omp parallel for
	for (size_t d = 0; d < K->domains.size(); ++d) {
		applyBt(this->feti, d, x, Btx[d]);
		KSpBlas[d].apply(KBtx[d], T{1}, T{0}, Btx[d]);
	}
	applyB(this->feti, KBtx, y);
}

}
