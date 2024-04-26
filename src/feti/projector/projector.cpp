
#include "projector.h"
#include "hfeti.orthogonal.symmetric.h"
#include "hfeti.orthogonal.symmetric.withfactors.h"
#include "tfeti.orthogonal.symmetric.h"
#include "tfeti.orthogonal.symmetric.withfactors.h"
#include "tfeti.conjugate.symmetric.h"
#include "esinfo/eslog.h"

namespace espreso {

template <typename T>
Projector<T>* Projector<T>::set(FETI<T> &feti, const step::Step &step)
{
    if (feti.configuration.iterative_solver == FETIConfiguration::ITERATIVE_SOLVER::SMALBE) {
        if (feti.configuration.projector == FETIConfiguration::PROJECTOR::ORTHOGONAL) {
            feti.configuration.projector = FETIConfiguration::PROJECTOR::ORTHOGONAL_WITH_FACTORS;
        }
    }

    switch (feti.configuration.method) {
    case FETIConfiguration::METHOD::TOTAL_FETI: {
        switch (feti.configuration.projector) {
        case FETIConfiguration::PROJECTOR::ORTHOGONAL:
        case FETIConfiguration::PROJECTOR::ORTHOGONAL_FULL:
            eslog::info(" = PROJECTOR                                                             EXPLICIT ORTHOGONAL = \n");
            return new TFETIOrthogonalSymmetric<T>(feti);
        case FETIConfiguration::PROJECTOR::ORTHOGONAL_WITH_FACTORS:
        case FETIConfiguration::PROJECTOR::ORTHOGONAL_FULL_WITH_FACTORS:
            eslog::info(" = PROJECTOR                                                EXPLICIT ORTHOGONAL WITH FACTORS = \n");
            return new TFETIOrthogonalSymmetricWithFactors<T>(feti);
        case FETIConfiguration::PROJECTOR::CONJUGATE:
            eslog::info(" = PROJECTOR                                                              EXPLICIT CONJUGATE = \n");
            return new TFETIConjugateSymmetric<T>(feti);
        default: return nullptr;
        }
    } break;
    case FETIConfiguration::METHOD::HYBRID_FETI: {
        switch (feti.configuration.projector) {
        case FETIConfiguration::PROJECTOR::ORTHOGONAL:
            eslog::info(" = PROJECTOR                                                      HYBRID EXPLICIT ORTHOGONAL = \n");
            return new HFETIOrthogonalSymmetric<T>(feti);
        case FETIConfiguration::PROJECTOR::ORTHOGONAL_WITH_FACTORS:
            eslog::info(" = PROJECTOR                                         HYBRID EXPLICIT ORTHOGONAL WITH FACTORS = \n");
            return new HFETIOrthogonalSymmetricWithFactors<T>(feti);
//        case FETIConfiguration::PROJECTOR::CONJUGATE:
//            eslog::info(" = PROJECTOR                                                        HYBRID EXPLICIT CONJUGATE = \n");
//            return new TFETIConjugateSymmetric<T>(feti);
        case FETIConfiguration::PROJECTOR::ORTHOGONAL_FULL:
            eslog::info(" = PROJECTOR                                                             EXPLICIT ORTHOGONAL = \n");
            return new TFETIOrthogonalSymmetric<T>(feti);
        case FETIConfiguration::PROJECTOR::ORTHOGONAL_FULL_WITH_FACTORS:
            eslog::info(" = PROJECTOR                                                EXPLICIT ORTHOGONAL WITH FACTORS = \n");
            return new TFETIOrthogonalSymmetricWithFactors<T>(feti);
        default: return nullptr;
        }
    } break;
    }
    return nullptr;
}

template struct Projector<double>;
template struct Projector<std::complex<double> >;

}
