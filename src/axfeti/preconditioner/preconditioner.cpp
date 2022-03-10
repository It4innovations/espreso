
#include "preconditioner.h"

#include "emptypreconditioner.h"

namespace espreso {

template <typename T>
static Preconditioner<T>* _set(AX_FETI<T> *feti)
{
	switch (feti->configuration.preconditioner) {
	case FETIConfiguration::PRECONDITIONER::NONE:
		eslog::info(" = PRECONDITIONER                                                                       NONE = \n");
		return new EmptyPreconditioner<T>(feti);
		break;
	case FETIConfiguration::PRECONDITIONER::LUMPED:
	case FETIConfiguration::PRECONDITIONER::WEIGHT_FUNCTION:
	case FETIConfiguration::PRECONDITIONER::DIRICHLET:
	case FETIConfiguration::PRECONDITIONER::SUPER_DIRICHLET:
	case FETIConfiguration::PRECONDITIONER::MAGIC:
	default: return nullptr;
	}
}

template <> Preconditioner<double>* Preconditioner<double>::set(AX_FETI<double> *feti) { return _set<double>(feti); }
template <> Preconditioner<std::complex<double> >* Preconditioner<std::complex<double> >::set(AX_FETI<std::complex<double> > *feti) { return _set<std::complex<double> >(feti); }

}
