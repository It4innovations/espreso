
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
            math::set(feti.RegMat[d].nrows + 1, feti.RegMat[d].rows, 1, feti.RegMat[d].rows[0]);
        }
    } else {
        #pragma omp parallel for
        for (size_t d = 0; d < feti.K.size(); ++d) {
            math::getKernel(feti.K[d], feti.R1[d], feti.RegMat[d], defect, sc_size);
        }
    }
}

template <typename T>
static void get_svd(Matrix_CSR<T> &K, Matrix_Dense<T> &R1, Matrix_Dense<T> &MoorePenroseInv)
{
    T eps = std::numeric_limits<T>::epsilon();
    T tol = std::max(K.nrows, K.ncols) * eps;

    Matrix_Dense<T> dK; dK.resize(K.nrows, K.ncols);
    math::copy(dK, K);

    Vector_Dense<T> s;
    Matrix_Dense<T> U, VT;
    math::lapack::get_svd(dK, s, U, VT);
    T s_max = *std::max_element(s.vals, s.vals + s.size);
    R1.nrows = 0;
    for (int i = 0; i < s.size; ++i) {
        if (s.vals[i] <= tol * s_max) {
            R1.nrows += 1;
        }
    }
    R1.resize(R1.nrows, K.ncols);
    std::copy(VT.vals + (VT.nrows - R1.nrows) * VT.ncols, VT.vals + VT.nnz, R1.vals);
    for (int r = 0; r < U.nrows; ++r) {
        for (int c = 0; c < U.ncols - R1.nrows; ++c) {
            U.vals[r * U.ncols + c] *= 1 / s.vals[c];
        }
        for (int c = U.ncols - R1.nrows; c < U.ncols; ++c) {
            U.vals[r * U.ncols + c] = 0;
        }
    }
    MoorePenroseInv.resize(dK.nrows, dK.ncols);
    math::blas::multiply(T{1}, U, VT, T{0}, MoorePenroseInv);
}

template <typename T>
void Regularization<T>::svd(FETI<T> &feti)
{
    if (feti.configuration.projector == FETIConfiguration::PROJECTOR::CONJUGATE) {
        #pragma omp parallel for
        for (size_t d = 0; d < feti.assembledK.size(); ++d) {
            get_svd(feti.assembledK[d], feti.R1[d], feti.MoorePenroseInv[d]);
        }
    } else {
        #pragma omp parallel for
        for (size_t d = 0; d < feti.K.size(); ++d) {
            get_svd(feti.K[d], feti.R1[d], feti.MoorePenroseInv[d]);
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
    feti.MoorePenroseInv.resize(feti.K.size());
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
