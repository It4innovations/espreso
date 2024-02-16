
#include "orthopcpg.h"

#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.hpp"
#include "math/math.h"
#include "feti/projector/projector.h"
#include "feti/dualoperator/dualoperator.h"
#include "feti/preconditioner/preconditioner.h"

namespace espreso {

// https://digital.library.unt.edu/ark:/67531/metadc739671/m2/1/high_res_d/792775.pdf
// page 15

// initialization
// l_0: L => Gt * inv(GGt) * e
// r_0: L => d - F * lambda_0
// w_0: L => r_0 - Gt * inv(GGt) * G * r_0 :: (I - Q) * r_0
// p_0: L => w_0

// loop
// gamma_k: 1 => (w_k,w_k) / (p_k, F * p_k)
//   x_k+1: L => x_k + gama_k * p_k
//   r_k+1: L => r_k - gama_k * F * p_k
//   w_k+1: L => r_k+1 - Gt * inv(GGt) * G * r_k+1 :: (I - Q) * r_k+1
//   y_k+1: L => w_k+1 - Gt * inv(GGt) * G * w_k+1 :: (I - Q) * w_k+1

//   p_k+1: L => w_k+1 - SUM{0->i}[ ((w_k+1, F * p_i) / (p_i, F * p_i)) * p_i ]

template <typename T>
void OrthogonalizedPCPG<T>::info()
{
    eslog::info(" = ORTHOGONAL PRECONDITIONED CONJUGATE PROJECTED GRADIENT SETTINGS                           = \n");
    switch (feti.configuration.stopping_criterion) {
    case FETIConfiguration::STOPPING_CRITERION::RELATIVE:
        eslog::info(" =   STOPPING CRITERION                                                             RELATIVE = \n");
        break;
    case FETIConfiguration::STOPPING_CRITERION::ABSOLUTE:
        eslog::info(" =   STOPPING CRITERION                                                             ABSOLUTE = \n");
        break;
    case FETIConfiguration::STOPPING_CRITERION::ARIOLI:
        eslog::info(" =   STOPPING CRITERION                                                               ARIOLI = \n");
        break;
    }
    eslog::info(" =   PRECISION                                                                      %.2e = \n", feti.configuration.precision);
    if (feti.configuration.max_iterations == 0) {
        eslog::info(" =   MAX_ITERATIONS                                                                     AUTO = \n");
    } else {
        eslog::info(" =   MAX_ITERATIONS                                                                  %7d = \n", feti.configuration.max_iterations);
    }
    eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
}

template <typename T>
OrthogonalizedPCPG<T>::OrthogonalizedPCPG(FETI<T> &feti)
: IterativeSolver<T>(feti)
{
    yFp.reserve(pi.initial_space);
    pFp.reserve(pi.initial_space);
}

template <> void OrthogonalizedPCPG<double>::solve(const step::Step &step, IterativeSolverInfo &info)
{
    DualOperator<double> *F = feti.dualOperator;
    Projector<double> *P = feti.projector;
    Preconditioner<double> *S = feti.preconditioner;

    Vector_Dual<double> p, Fp;
    pi.next(p);

    Vector_Dense<double> _yFp;

    P->applyGtInvGGt(P->e, l);             // l = Gt * inv(GGt) * e

    F->apply(l, r);                        // r = d - F * l
    math::scale(-1., r);                   //
    math::add(r, 1., F->d);                //

    P->apply(r, w);                        // w = P * r
    S->apply(w, z);                        // z = S * w
    P->apply(z, y);                        // y = P * z (y = P * S * w)

    math::copy(p, y);                      // p = y
    math::copy(x, l);                      // x = l

    double yw = y.dot(w);
    setInfo(info, feti.configuration, yw);

    eslog::checkpointln("FETI: PCPG INITIALIZATION");
    eslog::startln("ORTHOGONAL PCPG: ITERATIONS STARTED", "orthopcpg");
    while (!info.converged) {
        // gamma = (w, w) / (p, F * p)
        Fpi.next(Fp);
        F->apply(p, Fp);
        eslog::accumulatedln("orthopcpg: apply F");
        pFp.push_back(p.dot(Fp));
        double gamma = yw / pFp.back();
        eslog::accumulatedln("orthopcpg: dot(p, Fp)");

        // x = x + gamma * p
        // r = r - gamma * F * p
        math::add(x,  gamma,  p);
        math::add(r, -gamma, Fp);
        eslog::accumulatedln("orthopcpg: update x, r");

        // w = P * r
        P->apply(r, w);
        eslog::accumulatedln("orthopcpg: apply P");

        // z = S * w
        S->apply(w, z);
        eslog::accumulatedln("orthopcpg: apply S * w");

        // y = P * z
        P->apply(z, y);
        eslog::accumulatedln("orthopcpg: apply P * z");

        // p = w - SUM{0->i}[ ((w, F * p_i) / (p_i, F * p_i)) * p_i ]
        yFp.push_back(0); _yFp.vals = yFp.data(); _yFp.size = yFp.size();
        Fpi.apply(y, _yFp);
        for (size_t i = 0; i < yFp.size(); ++i) {
            _yFp.vals[i] /= -pFp[i];
        }
        pi.next(p);
        pi.applyT(_yFp, p);
        math::add(p, 1., y);
        eslog::accumulatedln("orthopcpg: orthogonalization");

        updateInfo(info, feti.configuration, yw, 0, 0);
        yw = y.dot(w);
        eslog::accumulatedln("orthopcpg: check criteria");
    }
    eslog::endln("orthopcpg: finished");
    eslog::checkpointln("FETI: ORTHOGONAL PCPG ITERATIONS");
    reconstructSolution(x, r);
    eslog::checkpointln("FETI: SOLUTION RECONSTRUCTION");
    eslog::info("       = ----------------------------------------------------------------------------- = \n");
}

template <> void OrthogonalizedPCPG<std::complex<double> >::solve(const step::Step &step, IterativeSolverInfo &info)
{

}

template class OrthogonalizedPCPG<double>;
template class OrthogonalizedPCPG<std::complex<double> >;

}

