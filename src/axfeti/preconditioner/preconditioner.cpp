
#include "preconditioner.h"

#include "dirichlet.h"
#include "emptypreconditioner.h"
#include "lumped.h"
#include "weightfunction.h"

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
		eslog::info(" = PRECONDITIONER                                                                     LUMPED = \n");
		return new Lumped<T>(feti);
		break;
	case FETIConfiguration::PRECONDITIONER::WEIGHT_FUNCTION:
		eslog::info(" = PRECONDITIONER                                                            WEIGHT FUNCTION = \n");
		return new WeightFunction<T>(feti);
		break;
	case FETIConfiguration::PRECONDITIONER::DIRICHLET:
		eslog::info(" = PRECONDITIONER                                                                  DIRICHLET = \n");
		return new Dirichlet<T>(feti);
		break;
	case FETIConfiguration::PRECONDITIONER::SUPER_DIRICHLET:
	case FETIConfiguration::PRECONDITIONER::MAGIC:
	default: return nullptr;
	}
}

template <> Preconditioner<double>* Preconditioner<double>::set(AX_FETI<double> *feti) { return _set<double>(feti); }
template <> Preconditioner<std::complex<double> >* Preconditioner<std::complex<double> >::set(AX_FETI<std::complex<double> > *feti) { return _set<std::complex<double> >(feti); }

}
