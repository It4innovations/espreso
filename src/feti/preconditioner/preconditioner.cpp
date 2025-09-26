
#include "preconditioner.h"

#include "emptypreconditioner.h"
#include "lumped.h"
#include "dirichlet.h"
#include "dirichlet.implicit.h"
#include "dirichlet.generalschur.h"
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
    case FETIConfiguration::PRECONDITIONER::DIRICHLET_IMPLICIT:
        eslog::info(" = PRECONDITIONER                                                         IMPLICIT DIRICHLET = \n");
        return new DirichletImplicit<T>(feti);
    case FETIConfiguration::PRECONDITIONER::DIRICHLET_GENERALSCHUR_CPU:
        eslog::info(" = PRECONDITIONER                          DIRICHLET GENERALIZED SCHUR IMPLEMENTATION ON CPU = \n");
        return new DirichletGeneralSchur<T,int>(feti, 'C');
    case FETIConfiguration::PRECONDITIONER::DIRICHLET_GENERALSCHUR_GPU:
        eslog::info(" = PRECONDITIONER                          DIRICHLET GENERALIZED SCHUR IMPLEMENTATION ON GPU = \n");
        if (!gpu::mgm::is_linked()) {
            eslog::globalerror("PRECONDITIONER::DIRICHLET_GENERALSCHUR_GPU cannot be created. GPU support is not built.\n");
        }
        if (!gpu::mgm::is_available()) {
            eslog::globalerror("PRECONDITIONER::DIRICHLET_GENERALSCHUR_GPU cannot be created. No GPUs detected or other error occured.\n");
        }
        return new DirichletGeneralSchur<T,int>(feti, 'G');
    default: return nullptr;
    }
}

template struct Preconditioner<double>;
template struct Preconditioner<std::complex<double> >;

}
