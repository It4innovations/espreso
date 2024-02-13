
#include "projector.h"
#include "tfeti.orthogonal.symmetric.h"
#include "tfeti.conjugate.symmetric.h"
#include "esinfo/eslog.h"

namespace espreso {

template <typename T>
Projector<T>* Projector<T>::set(FETI<T> &feti, const step::Step &step)
{
    switch (feti.configuration.projector) {
    case FETIConfiguration::PROJECTOR::ORTHOGONAL:
        eslog::info(" = PROJECTOR                                                                      ORTHOGONAL = \n");
        return new TFETIOrthogonalSymmetric<T>(feti);
    case FETIConfiguration::PROJECTOR::CONJUGATE:
        eslog::info(" = PROJECTOR                                                                       CONJUGATE = \n");
        return new TFETIConjugateSymmetric<T>(feti);
    default: return nullptr;
    }
}

template struct Projector<double>;
template struct Projector<std::complex<double> >;

}
