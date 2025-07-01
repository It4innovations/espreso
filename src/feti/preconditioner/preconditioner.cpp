
#include "preconditioner.h"

#include "emptypreconditioner.h"
#include "lumped.h"
#include "dirichlet.h"
#include "dirichlet.generalsc.cpu.h"
#include "esinfo/eslog.h"

#include <complex>

namespace espreso {

template <typename T>
Preconditioner<T>* Preconditioner<T>::create(FETI<T> &feti, const step::Step &step)
{
    switch (feti.configuration.preconditioner) {
    case FETIConfiguration::PRECONDITIONER::NONE:
        eslog::info(" = PRECONDITIONER                                                                       NONE = \n");
        return new EmptyPreconditioner<T>(feti);
    case FETIConfiguration::PRECONDITIONER::LUMPED:
        eslog::info(" = PRECONDITIONER                                                                     LUMPED = \n");
        return new Lumped<T>(feti);
    case FETIConfiguration::PRECONDITIONER::DIRICHLET:
        eslog::info(" = PRECONDITIONER                                                                  DIRICHLET = \n");
        return new Dirichlet<T>(feti);
    case FETIConfiguration::PRECONDITIONER::DIRICHLET_GENERALSC_CPU:
        eslog::info(" = PRECONDITIONER                             DIRICHLET GENERALIZED SC IMPLEMENTATION ON CPU = \n");
        return new DirichletGeneralScCpu<T,int>(feti);
    default: return nullptr;
    }
}

template struct Preconditioner<double>;
template struct Preconditioner<std::complex<double> >;

}
