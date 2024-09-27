
#include "regularization.h"
#include "regularization.heattransfer.h"
#include "regularization.elasticity.h"
#include "regularization.empty.h"

#include "esinfo/ecfinfo.h"
#include "math/math.h"

namespace espreso {

template <typename T>
void Regularization<T>::set(const step::Step &step, FETI<T> &feti)
{
    feti.R1.resize(feti.K.size());
    feti.R2.resize(feti.K.size());
    feti.RegMat.resize(feti.K.size());

    switch (feti.configuration.regularization) {
    case FETIConfiguration::REGULARIZATION::ANALYTIC:
        switch (info::ecf->physics) {
        case PhysicsConfiguration::TYPE::HEAT_TRANSFER:
            switch (info::ecf->heat_transfer.load_steps_settings.at(step.loadstep + 1).type) {
            case HeatTransferLoadStepConfiguration::TYPE::STEADY_STATE: RegularizationHeatTransfer<T>::set(feti); break;
            case HeatTransferLoadStepConfiguration::TYPE::TRANSIENT:    RegularizationEmpty<T>::set(feti); break;
            default: eslog::error("unknown regularization\n"); break;
            } break;
        case PhysicsConfiguration::TYPE::STRUCTURAL_MECHANICS:
            switch (info::ecf->structural_mechanics.load_steps_settings.at(step.loadstep + 1).type) {
            case StructuralMechanicsLoadStepConfiguration::TYPE::STEADY_STATE: RegularizationElasticity<T>::set(feti); break;
            case StructuralMechanicsLoadStepConfiguration::TYPE::TRANSIENT:    RegularizationEmpty<T>::set(feti); break;
            default: eslog::error("unknown regularization\n"); break;
            } break;
        }
        break;
    case FETIConfiguration::REGULARIZATION::ALGEBRAIC:
        // we need matrices K, set only sizes that we know surely
        for (size_t d = 0; d < feti.K.size(); ++d) {
            feti.R1[d].resize(0, feti.K[d].ncols);
            feti.RegMat[d].resize(feti.K[d].nrows, feti.K[d].ncols, 0);
            math::set(feti.RegMat[d].nrows + 1, feti.RegMat[d].rows, 1, Indexing::CSR);
        }
    }
}

template <typename T>
void Regularization<T>::update(const step::Step &step, FETI<T> &feti)
{
    switch (feti.configuration.regularization) {
    case FETIConfiguration::REGULARIZATION::ANALYTIC:
        switch (info::ecf->physics) {
        case PhysicsConfiguration::TYPE::HEAT_TRANSFER:
            switch (info::ecf->heat_transfer.load_steps_settings.at(step.loadstep + 1).type) {
            case HeatTransferLoadStepConfiguration::TYPE::STEADY_STATE: RegularizationHeatTransfer<T>::update(feti); break;
            case HeatTransferLoadStepConfiguration::TYPE::TRANSIENT:    RegularizationEmpty<T>::update(feti); break;
            default: break;
            } break;

        case PhysicsConfiguration::TYPE::STRUCTURAL_MECHANICS:
            switch (info::ecf->structural_mechanics.load_steps_settings.at(step.loadstep + 1).type) {
            case StructuralMechanicsLoadStepConfiguration::TYPE::STEADY_STATE: RegularizationElasticity<T>::update(feti); break;
            case StructuralMechanicsLoadStepConfiguration::TYPE::TRANSIENT:    RegularizationEmpty<T>::update(feti); break;
            default: break;
            }
            if (step.iteration) {
                // TODO: remove rotations and implement conjugate projector
                // TODO: keeping rotations is probably non-general solution that works for the cubic domains only
//                for (size_t d = 0; d < feti.K.size(); ++d) {
//                    math::blas::copy(feti.R1[d].ncols, feti.R1[d].vals + 3 * feti.R1[d].ncols, 1, feti.R1[d].vals + 4 * feti.R1[d].ncols, 1);
//                    feti.R1[d].resize(3, feti.K[d].ncols);
//                }
            }
            break;
        }
        break;
    case FETIConfiguration::REGULARIZATION::ALGEBRAIC:
        switch (info::ecf->physics) {
        case PhysicsConfiguration::TYPE::HEAT_TRANSFER:
            for (size_t d = 0; d < feti.K.size(); ++d) {
                math::getKernel(feti.K[d], feti.R1[d], feti.RegMat[d], 1, (int)feti.configuration.sc_size);
            }
            break;
        case PhysicsConfiguration::TYPE::STRUCTURAL_MECHANICS:
            switch (info::ecf->structural_mechanics.load_steps_settings.at(step.loadstep + 1).type) {
            case HeatTransferLoadStepConfiguration::TYPE::STEADY_STATE:
            case HeatTransferLoadStepConfiguration::TYPE::TRANSIENT:
                for (size_t d = 0; d < feti.K.size(); ++d) {
                    math::getKernel(feti.K[d], feti.R1[d], feti.RegMat[d], 6, (int)feti.configuration.sc_size);
                }
                break;
            case HeatTransferLoadStepConfiguration::TYPE::HARMONIC:
                for (size_t d = 0; d < feti.K.size(); ++d) {
                    math::getKernel(feti.K[d], feti.R1[d], feti.RegMat[d], 12, (int)feti.configuration.sc_size);
                }
                break;
            }
        }
        break;
    }

    if (feti.updated.K && feti.configuration.regularization == FETIConfiguration::REGULARIZATION::ANALYTIC) {
        std::vector<int> offset;
        Matrix_Dense<T> _R1;
        bool fetiProjection =
                (feti.configuration.method == FETIConfiguration::METHOD::TOTAL_FETI) ||
                (feti.configuration.projector_opt & FETIConfiguration::PROJECTOR_OPT::FULL);

        if (fetiProjection) {
            #pragma omp parallel for
            for (size_t d = 0; d < feti.R1.size(); ++d) {
                math::orthonormalize(feti.R1[d]);
            }
        } else {
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
        }
    }
}

template struct Regularization<double>;

}
