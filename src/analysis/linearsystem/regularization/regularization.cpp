
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
void Regularization<T>::orthonormalize(FETI<T> &feti)
{
    if (!feti.updated.K) {
        return;
    }
    if (feti.configuration.regularization == FETIConfiguration::REGULARIZATION::ALGEBRAIC) {
        // automatically orthogonal
        return;
    }

    bool cluster = (feti.configuration.method == FETIConfiguration::METHOD::TOTAL_FETI) || (feti.configuration.projector_opt & FETIConfiguration::PROJECTOR_OPT::FULL);

    std::vector<int> offset;
    Matrix_Dense<T> _R1;

    if (cluster) {
        for (size_t d = 0; d < feti.R1.size(); ++d) {
            offset.push_back(_R1.ncols);
            _R1.nrows = std::max(_R1.nrows, feti.R1[d].nrows);
            _R1.ncols += feti.R1[d].ncols;
        }
        _R1.resize(_R1.nrows, _R1.ncols);
        for (size_t d = 0; d < feti.R1.size(); ++d) {
            for (int r = 0; r < feti.R1[d].nrows; ++r) {
                math::blas::copy(feti.R1[d].ncols, _R1.vals + _R1.ncols * r + offset[d], 1, feti.R1[d].vals + feti.R1[d].ncols * r, 1);
            }
        }
        math::orthonormalize(_R1);
        for (size_t d = 0; d < feti.R1.size(); ++d) {
            for (int r = 0; r < feti.R1[d].nrows; ++r) {
                math::blas::copy(feti.R1[d].ncols, feti.R1[d].vals + feti.R1[d].ncols * r, 1, _R1.vals + _R1.ncols * r + offset[d], 1);
            }
        }
    } else {
        #pragma omp parallel for
        for (size_t d = 0; d < feti.R1.size(); ++d) {
            math::orthonormalize(feti.R1[d]);
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
void Regularization<T>::update(const step::Step &step, FETI<T> &feti)
{
    switch (info::ecf->physics) {
    case PhysicsConfiguration::TYPE::HEAT_TRANSFER:
        update(feti, info::ecf->heat_transfer.load_steps_settings.at(step.loadstep + 1));
        break;
    case PhysicsConfiguration::TYPE::STRUCTURAL_MECHANICS:
        update(feti, info::ecf->structural_mechanics.load_steps_settings.at(step.loadstep + 1));
        break;
    }
}

template struct Regularization<double>;

}
