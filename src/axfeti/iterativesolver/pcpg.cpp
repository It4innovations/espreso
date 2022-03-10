
#include "pcpg.h"

namespace espreso {

template <typename T>
static void _info(PCPG<T> *solver)
{

}

template <typename T>
static void _set(PCPG<T> *solver)
{

}

template <typename T>
static void _update(PCPG<T> *solver)
{

}

template <> PCPG<double>::PCPG(AX_FETI<double> *feti): IterativeSolver<double>(feti) { _set<double>(this); }
template <> PCPG<std::complex<double> >::PCPG(AX_FETI<std::complex<double> > *feti): IterativeSolver<std::complex<double> >(feti) { _set<std::complex<double> >(this); }

template <> void PCPG<double>::info() { _info<double>(this); }
template <> void PCPG<std::complex<double> >::info() { _info<std::complex<double> >(this); }

template <> void PCPG<double>::update() { _update<double>(this); }
template <> void PCPG<std::complex<double> >::update() { _update<std::complex<double> >(this); }

}

