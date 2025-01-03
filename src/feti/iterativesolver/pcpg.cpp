
#include "pcpg.h"

#include "basis/utilities/sysutils.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.hpp"
#include "math/math.h"
#include "feti/projector/projector.h"
#include "feti/dualoperator/dualoperator.h"
#include "feti/preconditioner/preconditioner.h"

namespace espreso {

// initialization
// l_0: L => Gt * inv(GGt) * e
// r_0: L => d - F * lambda_0
// w_0: L => r_0 - Gt * inv(GGt) * G * r_0 :: (I - Q) * r_0
// y_0: L -> r_0 - Gt * inv(GGt) * G * S * w_0 :: (I - Q) * S * w_0
// p_0: L => w_0

// loop
// gamma_k: 1 => (y_k,w_k) / (p_k, F * p_k)
//   x_k+1: L => x_k + gama_k * p_k
//   r_k+1: L => r_k - gama_k * F * p_k
//   w_k+1: L => r_k+1 - Gt * inv(GGt) * G * r_k+1 :: (I - Q) * r_k+1
//   y_k+1: L => (I - Q) * S * w_k+1
//  beta_k: 1 => (y_k+1,w_k+1) / (p_k,w_k)
//   p_k+1: L => y_k+1 + beta_k * p_k

template <typename T>
void PCPG<T>::info()
{
    eslog::info(" = PRECONDITIONED CONJUGATE PROJECTED GRADIENT SETTINGS                                      = \n");
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
PCPG<T>::PCPG(FETI<T> &feti)
: IterativeSolver<T>(feti)
{

}

template <typename T>
static void _print(const char *name, const IterativeSolverInfo &info, const step::Step &step, const Vector_Dual<T> &v)
{
    if (info::ecf->output.print_matrices > 1) {
        eslog::storedata(" STORE: feti/pcpg/{%s%s}\n", name, std::to_string(info.iterations).c_str());
        math::store(v, utils::filename(utils::debugDirectory(step) + "/feti/pcpg", std::string(name) + std::to_string(info.iterations)).c_str());
    }
}

template <> void PCPG<double>::solve(const step::Step &step, IterativeSolverInfo &info)
{
    DualOperator<double> *F = feti.dualOperator;
    Projector<double> *P = feti.projector;
    Preconditioner<double> *S = feti.preconditioner;

    if (feti.configuration.projector == FETIConfiguration::PROJECTOR::CONJUGATE) {
        P->apply_GtintGGtG(F->d, l);       // l = Gt * inv(GFGt) * G * d
    } else {
        P->apply_GtintGGt(P->e, l);        // l = Gt * inv(GGt) * e
    }

    F->apply(l, r);                        // r = d - F * l
    math::scale(-1., r);                   //
    math::add(r, 1., F->d);                //

    P->applyT(r, w);                       // w = P * r
    S->apply(w, z);                        // z = S * w
    P->apply(z, y);                        // y = P * z (y = P * S * w)

    math::copy(p, y);                      // p = w
    math::copy(x, l);                      // x = l

    _print("p", info, step, p);
    _print("x", info, step, x);
    _print("r", info, step, r);
    _print("z", info, step, x);

    double yw = y.dot(w);
    setInfo(info, feti.configuration, yw);

    eslog::checkpointln("FETI: CPG INITIALIZATION");
    eslog::startln("PCPG: ITERATIONS STARTED", "pcpg");
    while (!info.converged) {
        // gamma = (y, w) / (p, F * p)
        F->apply(p, Fp);
        eslog::accumulatedln("pcpg: apply F");
        double pFp = p.dot(Fp), gamma = yw / pFp;
        eslog::accumulatedln("pcpg: dot(p, Fp)");

        // x = x + gamma * p
        // r = r - gamma * F * p
        math::add(x,  gamma,  p);
        math::add(r, -gamma, Fp);
        eslog::accumulatedln("pcpg: update x, r");

        // w = Pt * r
        P->applyT(r, w);
        eslog::accumulatedln("pcpg: apply P * r");

            // z = S * w
        S->apply(w, z);
        eslog::accumulatedln("pcpg: apply S * w");

        // y = P * z
        P->apply(z, y);
        eslog::accumulatedln("pcpg: apply P * z");


        // beta = (y+1, w+1) / (y, w)
        double _yw = y.dot(w), beta = _yw / yw;
        eslog::accumulatedln("pcpg: dot(y, w)");

        // p = y + beta * p  (y is not used anymore)
        math::add(y, beta, p); y.swap(p);
        eslog::accumulatedln("pcpg: update p");

        updateInfo(info, feti.configuration, yw, 0, 0);
        yw = _yw; // keep yw for the next iteration
        eslog::accumulatedln("pcpg: check criteria");
    }
    eslog::endln("pcpg: finished");
    eslog::checkpointln("FETI: PCPG ITERATIONS");
    reconstructSolution(x, r, step);
    eslog::checkpointln("FETI: SOLUTION RECONSTRUCTION");
    eslog::info("       = ----------------------------------------------------------------------------- = \n");
}

template <> void PCPG<std::complex<double> >::solve(const step::Step &step, IterativeSolverInfo &info)
{

}

template class PCPG<double>;
template class PCPG<std::complex<double> >;

}


