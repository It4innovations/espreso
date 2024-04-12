
#include "dualoperator.h"
#include "totalfeti.implicit.h"
#include "totalfeti.explicit.h"
#include "totalfeti.explicit.acc.h"
#include "feti/projector/projector.h"

#include "esinfo/eslog.hpp"
#include "math/math.h"

#include <climits>

namespace espreso {

template<typename T>
DualOperator<T>* DualOperator<T>::set(FETI<T> &feti, const step::Step &step)
{
    DualOperator<T>* dual = nullptr;
    switch (feti.configuration.dual_operator) {
    case FETIConfiguration::DUAL_OPERATOR::IMPLICIT:
        eslog::info(" = DUAL OPERATOR                                                         IMPLICIT TOTAL FETI = \n");
        dual = new TotalFETIImplicit<T>(feti);
        break;
    case FETIConfiguration::DUAL_OPERATOR::EXPLICIT:
        eslog::info(" = DUAL OPERATOR                                                         EXPLICIT TOTAL FETI = \n");
        dual = new TotalFETIExplicit<T>(feti);
        break;
    case FETIConfiguration::DUAL_OPERATOR::EXPLICIT_GPU:
        if (!gpu::mgm::is_linked()) {
            eslog::globalerror("GPU acceleration is not supported: GPU support in not built.\n");
        }
        if (!DirectSparseSolver<T>::provideFactors()) {
            eslog::globalerror("GPU acceleration is not supported: Third party sparse solver does not provide factors.\n");
        }
        eslog::info(" = DUAL OPERATOR                                                  EXPLICIT TOTAL FETI ON GPU = \n");
        dual = new TotalFETIExplicitAcc<T,int>(feti);
        break;
    }
    dual->set(step);
    return dual;
}

template<typename T>
void DualOperator<T>::getInitVector(Vector_Dual<T> &v)
{
    double norm = 1 / std::sqrt(feti.sinfo.dual_total);
    int i = 0;
    for (int j = 0; j < feti.lambdas.eq_halo; ++i, ++j) {
        v.vals[i] = 0;
    }
    for (int j = 0; j < feti.lambdas.eq_size; ++i, ++j) {
        int mul = (feti.sinfo.eq_offset + j) % 2;
        v.vals[i] = (1 - 2 * mul) * norm;
    }
    for (int j = 0; j < feti.lambdas.nc_halo; ++i, ++j) {
        v.vals[i] = 0;
    }
    for (int j = 0; j < feti.lambdas.nc_size; ++i, ++j) {
        int mul = (feti.sinfo.eq_total + feti.sinfo.nc_offset + j) % 2;
        v.vals[i] = (1 - 2 * mul) * norm;
    }
    v.synchronize();
}

template<typename T>
void DualOperator<T>::estimateMaxEigenValue(double &lambda, int &iterations, double epsilon, int maxIterations)
{
    DualOperator<T> *F = feti.dualOperator;

    // {-1, 1} / norma
    Vector_Dual<T> y, v;
    getInitVector(y);

    lambda = std::sqrt(v.dot());
    double err = std::numeric_limits<T>::max();
    for (iterations = 0; iterations <= maxIterations && epsilon < err; ++iterations) {
        math::copy(y, v);
        math::scale(1 / lambda, y);
        F->apply(y, v);
        double _lambda = lambda;
        lambda = std::sqrt(v.dot());
        err = std::fabs(lambda - _lambda) / std::fabs(lambda);
    }
}

template<typename T>
void DualOperator<T>::estimateMaxProjectedEigenValue(double &lambda, int &iterations, double epsilon, int maxIterations, double rho, double normPFP)
{
    DualOperator<T> *F = feti.dualOperator;
    Projector<double> *P = feti.projector;

    // {-1, 1} / norma
    Vector_Dual<T> v, y, x, z;
    getInitVector(y);
    lambda = std::sqrt(y.dot());
    double err = std::numeric_limits<T>::max();
    for (iterations = 1; iterations <= maxIterations && epsilon < err; ++iterations) {
        math::copy(x, y);
        math::scale(1 / lambda, x);

        // y = H(x)
        P->apply(x, v);
        F->apply(v, z);
        P->apply(z, y);
        double _lambda = lambda;
        lambda = std::sqrt(y.dot());
        err = std::fabs(lambda - _lambda) / std::fabs(lambda);
    }
}

template class DualOperator<double>;

}
