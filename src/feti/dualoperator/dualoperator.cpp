
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
        if (DirectSparseSolver<T>::provideFactors()) {
            eslog::info(" = DUAL OPERATOR                                                  EXPLICIT TOTAL FETI ON GPU = \n");
            dual = new TotalFETIExplicitAcc<T,esint>(feti);
        } else {
            eslog::globalerror("Third party software problem: solver does not provide factors that are required for EXPLICIT TOTAL FETI ON GPU.\n");
        }
        break;
    }
    dual->set(step);
    return dual;
}

template<typename T>
void DualOperator<T>::estimateMaxEigenValue(double &lambda, int &iterations, double epsilon, int maxIterations)
{
    DualOperator<T> *F = feti.dualOperator;

    // {-1, 1} / norma
    Vector_Dual<T> y, v;
    double norm = 1 / std::sqrt(feti.sinfo.lambdasTotal);
    for (esint i = 0; i < feti.lambdas.eq_halo; ++i) {
        v.vals[i] = 0;
    }
    for (esint i = feti.lambdas.eq_halo, j = 0; i < v.size; ++i, ++j) {
        int mul = (feti.sinfo.lambdasOffset + j) % 2;
        v.vals[i] = (1 - 2 * mul) * norm;
    }
    v.synchronize();

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
    double norm = 1 / std::sqrt(feti.sinfo.lambdasTotal);
    for (esint i = 0; i < feti.lambdas.eq_halo; ++i) {
        y.vals[i] = 0;
    }
    for (esint i = feti.lambdas.eq_halo, j = 0; i < y.size; ++i, ++j) {
        int mul = (feti.sinfo.lambdasOffset + j) % 2;
        y.vals[i] = (1 - 2 * mul) * norm;
    }
    y.synchronize();

    lambda = std::sqrt(y.dot());
    double err = std::numeric_limits<T>::max();
    for (iterations = 1; iterations <= maxIterations && epsilon < err; ++iterations) {
        math::copy(x, y);
        math::scale(1 / lambda, x);

        // y = H(x)
        P->apply(x, v);
        F->apply(v, z);
        P->apply(z, y);
//        math::scale(1 / normPFP, y);
//        math::add(y,  rho, x);
//        math::add(y, -rho, v);

        double _lambda = lambda;
        lambda = std::sqrt(y.dot());
        err = std::fabs(lambda - _lambda) / std::fabs(lambda);
    }
}

template class DualOperator<double>;

}
