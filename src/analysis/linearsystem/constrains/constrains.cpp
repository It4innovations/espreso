
#include "constrains.h"
#include "equalityconstrains.h"
#include "fixedwall.h"

#include "esinfo/ecfinfo.h"
#include "math/math.h"

namespace espreso {

template <typename T>
void Constrains<T>::set(const step::Step &step, FETI<T> &feti, const Vector_Distributed<Vector_Sparse, T> &dirichlet)
{
    eq.set(step, feti, dirichlet);
    switch (info::ecf->physics) {
    case PhysicsConfiguration::TYPE::HEAT_TRANSFER:
        mortar.set(step, feti, 1);
        break;
    case PhysicsConfiguration::TYPE::STRUCTURAL_MECHANICS:
        mortar.set(step, feti, 3);
        fw.set(step, feti, dirichlet);
        break;
    }

    if (feti.configuration.method == FETIConfiguration::METHOD::HYBRID_FETI) {
        cfg.set(step, feti);
    }
}

template <typename T>
void Constrains<T>::update(const step::Step &step, FETI<T> &feti, const Vector_Distributed<Vector_Sparse, T> &dirichlet)
{
    eq.update(step, feti, dirichlet);
    switch (info::ecf->physics) {
    case PhysicsConfiguration::TYPE::HEAT_TRANSFER:
        mortar.update(step, feti);
        break;
    case PhysicsConfiguration::TYPE::STRUCTURAL_MECHANICS:
        mortar.update(step, feti);
        fw.update(step, feti, dirichlet);
        break;
    }

    if (feti.configuration.method == FETIConfiguration::METHOD::HYBRID_FETI) {
        cfg.update(step, feti);
    }
}

template struct Constrains<double>;

}

