
#include "regularization.h"

#include "esinfo/ecfinfo.h"
#include "math/math.h"

namespace espreso {

template <typename T>
void Regularization<T>::empty(FETI<T> &feti)
{
    #pragma omp parallel for
    for (size_t d = 0; d < feti.K.size(); ++d) {
        feti.R1[d].resize(0, feti.K[d].nrows);
        feti.R1[d].type = Matrix_Type::REAL_NONSYMMETRIC;
        feti.R1[d].shape = Matrix_Shape::FULL;

        feti.R2[d].resize(0, feti.K[d].nrows);
        feti.R2[d].type = Matrix_Type::REAL_NONSYMMETRIC;
        feti.R2[d].shape = Matrix_Shape::FULL;

        feti.RegMat[d].resize(feti.K[d].nrows, feti.K[d].ncols, 0);
        feti.RegMat[d].type = feti.K[d].type;
        feti.RegMat[d].shape = feti.K[d].shape;
        for (esint r = 0; r <= feti.RegMat[d].nrows; ++r) {
            feti.RegMat[d].rows[r] = 0;
        }
    }
}

template <typename T>
void Regularization<T>::algebraic(FETI<T> &feti, int defect, int sc_size)
{
    if (feti.configuration.projector == FETIConfiguration::PROJECTOR::CONJUGATE) {
        #pragma omp parallel for
        for (size_t d = 0; d < feti.assembledK.size(); ++d) {
            math::getKernel(feti.assembledK[d], feti.R1[d], feti.RegMat[d], defect, sc_size);
            feti.RegMat[d].resize(feti.K[d].nrows, feti.K[d].nrows, 0);
        }
    } else {
        #pragma omp parallel for
        for (size_t d = 0; d < feti.K.size(); ++d) {
            math::getKernel(feti.K[d], feti.R1[d], feti.RegMat[d], defect, sc_size);
        }
    }
}

template <typename T>
template <typename Settings, typename Configuration>
void Regularization<T>::analyze(FETI<T> &feti, Settings &settings, Configuration &configuration)
{
    regMat = R1 = R2 = onSurface = false;
    if (configuration.type == LoadStepSolverConfiguration::TYPE::STEADY_STATE) {
        regMat = R1 = true;
    }
    if (configuration.type == LoadStepSolverConfiguration::TYPE::TRANSIENT && feti.configuration.projector == FETIConfiguration::PROJECTOR::CONJUGATE) {
        R1 = true;
    }

    for (auto disc = settings.discretization.cbegin(); disc != settings.discretization.cend(); ++disc) {
        if (disc->second == PhysicsConfiguration::DISCRETIZATION::BEM) {
            onSurface = true;
            break;
        }
    }
}

template <typename T>
void Regularization<T>::set(const step::Step &step, FETI<T> &feti)
{
    feti.R1.resize(feti.K.size());
    feti.R2.resize(feti.K.size());
    feti.RegMat.resize(feti.K.size());
    empty(feti);

    switch (info::ecf->physics) {
    case PhysicsConfiguration::TYPE::HEAT_TRANSFER:
        analyze(feti, info::ecf->heat_transfer, info::ecf->heat_transfer.load_steps_settings.at(step.loadstep + 1));
        set(feti, info::ecf->heat_transfer.load_steps_settings.at(step.loadstep + 1));
        break;
    case PhysicsConfiguration::TYPE::STRUCTURAL_MECHANICS:
        analyze(feti, info::ecf->structural_mechanics, info::ecf->structural_mechanics.load_steps_settings.at(step.loadstep + 1));
        set(feti, info::ecf->structural_mechanics.load_steps_settings.at(step.loadstep + 1));
        break;
    }
}

template <typename T>
void Regularization<T>::update(const step::Step &step, FETI<T> &feti, Vector_Distributed<Vector_Dense, T> *solution)
{
    switch (info::ecf->physics) {
    case PhysicsConfiguration::TYPE::HEAT_TRANSFER:
        update(feti, info::ecf->heat_transfer.load_steps_settings.at(step.loadstep + 1));
        break;
    case PhysicsConfiguration::TYPE::STRUCTURAL_MECHANICS:
        update(feti, info::ecf->structural_mechanics.load_steps_settings.at(step.loadstep + 1), solution);
        break;
    }
}

template struct Regularization<double>;

}
