
#include "orthocpg.h"

#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.hpp"
#include "math/math.h"
#include "feti/projector/projector.h"
#include "feti/dualoperator/dualoperator.h"

namespace espreso {

// https://digital.library.unt.edu/ark:/67531/metadc739671/m2/1/high_res_d/792775.pdf
// page 14

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

//   p_k+1: L => w_k+1 - SUM{0->i}[ ((w_k+1, F * p_i) / (p_i, F * p_i)) * p_i ]

template <typename T>
OrthogonalizedCPG<T>::OrthogonalizedCPG(FETI<T> &feti)
: IterativeSolver<T>(feti)
{
    wFp.reserve(pi.initial_space);
    pFp.reserve(pi.initial_space);
}

template <typename T>
void OrthogonalizedCPG<T>::info()
{
    eslog::info(" = ORTHOGONAL CONJUGATE PROJECTED GRADIENT SETTINGS                                          = \n");
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

template <> void OrthogonalizedCPG<double>::solve(const step::Step &step, IterativeSolverInfo &info)
{
    DualOperator<double> *F = feti.dualOperator;
    Projector<double> *P = feti.projector;

    Vector_Dual<double> p, Fp;
    pi.next(p);

    Vector_Dense<double> _wFp;

    P->apply_e(P->e, l);             // l = Gt * inv(GGt) * e

    F->apply(l, r);                        // r = d - F * l
    math::scale(-1., r);                   //
    math::add(r, 1., F->d);                //

    P->apply(r, w);                        // w = P * r
    math::copy(p, w);                      // p = w
    math::copy(x, l);                      // x = l

    double ww = w.dot();
    setInfo(info, feti.configuration, ww);

    eslog::checkpointln("FETI: CPG INITIALIZATION");
    eslog::startln("ORTHOGONAL CPG: ITERATIONS STARTED", "orthocpg");
    while (!info.converged) {
        // gamma = (w, w) / (p, F * p)
        Fpi.next(Fp);
        F->apply(p, Fp);
        eslog::accumulatedln("orthocpg: apply F");
        pFp.push_back(p.dot(Fp));
        double gamma = ww / pFp.back();
        eslog::accumulatedln("orthocpg: dot(p, Fp)");

        math::add(x,  gamma,  p);  // x = x + gamma * p
        math::add(r, -gamma, Fp);  // r = r - gamma * F * p
        eslog::accumulatedln("orthocpg: update x, r");

        // w = P * r
        P->apply(r, w);
        eslog::accumulatedln("orthocpg: apply P");

        // p = w - SUM{0->i}[ ((w, F * p_i) / (p_i, F * p_i)) * p_i ]
        wFp.push_back(0); _wFp.vals = wFp.data(); _wFp.size = wFp.size();
        Fpi.apply(w, _wFp);
        for (size_t i = 0; i < wFp.size(); ++i) {
            _wFp.vals[i] /= -pFp[i];
        }
        pi.next(p);
        pi.applyT(_wFp, p);
        math::add(p, 1., w);
        eslog::accumulatedln("orthocpg: orthogonalization");

        updateInfo(info, feti.configuration, ww, 0, 0);
        ww = w.dot();
        eslog::accumulatedln("orthocpg: check criteria");
    }
    eslog::endln("orthocpg: finished");
    eslog::checkpointln("FETI: ORTHOGONAL CPG ITERATIONS");
    reconstructSolution(x, r, step);
    eslog::checkpointln("FETI: SOLUTION RECONSTRUCTION");
    eslog::info("       = ----------------------------------------------------------------------------- = \n");
}

template <> void OrthogonalizedCPG<std::complex<double> >::solve(const step::Step &step, IterativeSolverInfo &info)
{

}

template class OrthogonalizedCPG<double>;
template class OrthogonalizedCPG<std::complex<double> >;

}


