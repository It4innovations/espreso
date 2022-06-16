
#include "preconditioner.h"

#include "dirichlet.h"
#include "emptypreconditioner.h"
#include "lumped.h"
#include "weightfunction.h"

namespace espreso {

template <typename T>
static Preconditioner<T>* _set(FETI<T> *feti)
{
	switch (feti->configuration.preconditioner) {
	case FETIConfiguration::PRECONDITIONER::NONE:
		eslog::info(" = PRECONDITIONER                                                                       NONE = \n");
		return new EmptyPreconditioner<T>(feti);
	case FETIConfiguration::PRECONDITIONER::LUMPED:
		eslog::info(" = PRECONDITIONER                                                                     LUMPED = \n");
		return new Lumped<T>(feti);
	case FETIConfiguration::PRECONDITIONER::WEIGHT_FUNCTION:
		eslog::info(" = PRECONDITIONER                                                            WEIGHT FUNCTION = \n");
		return new WeightFunction<T>(feti);
	case FETIConfiguration::PRECONDITIONER::DIRICHLET:
		eslog::info(" = PRECONDITIONER                                                                  DIRICHLET = \n");
		return new Dirichlet<T>(feti);
	case FETIConfiguration::PRECONDITIONER::SUPER_DIRICHLET:
	case FETIConfiguration::PRECONDITIONER::MAGIC:
	default: return nullptr;
	}
}

template <> Preconditioner<double>* Preconditioner<double>::set(FETI<double> *feti) { return _set<double>(feti); }
template <> Preconditioner<std::complex<double> >* Preconditioner<std::complex<double> >::set(FETI<std::complex<double> > *feti) { return _set<std::complex<double> >(feti); }

}
