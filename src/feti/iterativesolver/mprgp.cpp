
#include "mprgp.h"

#include "basis/utilities/sysutils.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.hpp"
#include "math/math.h"
#include "feti/projector/projector.h"
#include "feti/preconditioner/preconditioner.h"
#include "feti/dualoperator/dualoperator.h"
#include "wrappers/mpi/communication.h"

namespace espreso {

// https://doi.org/10.1007/b138610

template <typename T>
MPRGP<T>::MPRGP(FETI<T> &feti)
: IterativeSolver<T>(feti)
{
    Vector_Dual<int>::initBuffers();
}

template <typename T>
void MPRGP<T>::info()
{
    eslog::info(" = MODIFIED PROPORTIONING WITH REDUCED GRADIENT PROJECTION SETTING                           = \n");
//    switch (feti.configuration.stopping_criterion) {
//    case FETIConfiguration::STOPPING_CRITERION::RELATIVE:
//        eslog::info(" =   STOPPING CRITERION                                                             RELATIVE = \n");
//        break;
//    case FETIConfiguration::STOPPING_CRITERION::ABSOLUTE:
//        eslog::info(" =   STOPPING CRITERION                                                             ABSOLUTE = \n");
//        break;
//    case FETIConfiguration::STOPPING_CRITERION::ARIOLI:
//        eslog::info(" =   STOPPING CRITERION                                                               ARIOLI = \n");
//        break;
//    }
    eslog::info(" =   PRECISION                                                                      %.2e = \n", feti.configuration.precision);
    eslog::info(" =   MAX_ITERATIONS                                                                  %7d = \n", feti.configuration.max_iterations);
    eslog::info(" =   ALPHAM                                                                         %.2e = \n", feti.configuration.alpham);
    eslog::info(" =   DELTA                                                                          %.2e = \n", feti.configuration.delta);
    eslog::info(" =   GRADPROJ                                                                        %7d = \n", feti.configuration.gradproj);
    eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
}

template <typename T>
static void _print(const char *name, const IterativeSolverInfo &info, const step::Step &step, const Vector_Dual<T> &v)
{
    if (info::ecf->output.print_matrices > 1) {
        eslog::storedata(" STORE: feti/MPRGP/{%s%s}\n", name, std::to_string(info.iterations).c_str());
        math::store(v, utils::filename(utils::debugDirectory(step) + "/feti/MPRGP", std::string(name) + std::to_string(info.iterations)).c_str());
    }
}

template <>
void MPRGP<double>::restrictToFeasibleSet(Vector_Dual<double> &x)
{
    for (int i = feti.lambdas.equalities, j = 0; i < feti.lambdas.size; ++i, ++j) {
        x.vals[i] = std::max(feti.lb.vals[j], std::min(feti.ub.vals[j], x.vals[i]));
    }
}

template <>
void MPRGP<double>::updateFreeAndActiveSet(Vector_Dual<int> &free, Vector_Dual<int> &active, const Vector_Dual<double> &x)
{
    math::set(free, 1);
    math::set(active, 0);
    for (int i = feti.lambdas.equalities, j = 0; i < feti.lambdas.size; ++i, ++j) {
        free.vals[i] = (x.vals[i] - feti.lb.vals[j] > feti.configuration.precision_set && feti.ub.vals[j] - x.vals[i] > feti.configuration.precision_set) ? 1 : 0;
        active.vals[i] = free.vals[i] ? 0 : 1;
    }
}

template <>
void MPRGP<double>::updateReducedGradient(Vector_Dual<double> &g_red, const Vector_Dual<double> &g, const Vector_Dual<double> &x, double alpha)
{
    math::copy(xp, x);
    math::add(xp, -alpha, g);
    restrictToFeasibleSet(xp);
    math::copy(g_red, x);
    math::add(g_red, -1., xp);
    math::scale(1 / alpha, g_red);
    for (int i = feti.lambdas.equalities; i < feti.lambdas.size; ++i) {
        if (g_red.vals[i] * g.vals[i] < 0) {
            g_red.vals[i] = 0;
        }
    }
//    for (int i = 0; i < g_red.size; ++i) {
//        fi_red.vals[i] = free.vals[i] ? g_red.vals[i] : 0;
//        beta_red.vals[i] = active.vals[i] ? g_red.vals[i] : 0;
//    }
}

template <>
void MPRGP<double>::updateStoppingGradient(Vector_Dual<double> &g_stop, Vector_Dual<double> &g, Vector_Dual<double> &x, double alpha)
{
    updateReducedGradient(g_stop, g, x, alpha);
}

template <>
void MPRGP<double>::updateFreeGradient(Vector_Dual<double> &g_free, Vector_Dual<double> &g, Vector_Dual<int> &free)
{
    for (int i = 0; i < feti.lambdas.size; ++i) {
        g_free.vals[i] = free.vals[i] ? g.vals[i] : 0;
    }
}

template <>
void MPRGP<double>::multByFree(Vector_Dual<double> &z, Vector_Dual<int> &free)
{
    for (int i = 0; i < z.size; ++i) {
        z.vals[i] = free.vals[i] * z.vals[i];
    }
}

template <>
double MPRGP<double>::getFeasibleStepLength(Vector_Dual<double> &x, Vector_Dual<double> &p)
{
    double alpha = std::numeric_limits<double>::max();
    for (int i = feti.lambdas.equalities, j = 0; i < feti.lambdas.size; ++i, ++j) {
        if (p.vals[i] > 0) {
            alpha = std::min(alpha, (x.vals[i] - feti.lb.vals[j]) / p.vals[i]);
        }
        if (p.vals[i] < 0) {
            alpha = std::min(alpha, (x.vals[i] - feti.ub.vals[j]) / p.vals[i]);
        }
    }
    Communication::allReduce(&alpha, nullptr, 1, MPITools::getType(alpha).mpitype, MPI_MIN);
    return alpha;
}

template <> void MPRGP<double>::run(const step::Step &step, MPRGPSolverInfo &info, double alpha, std::function<void(Vector_Dual<double> &in, Vector_Dual<double> &out)> H, std::function<bool(const Vector_Dual<double> &x, const Vector_Dual<double> &g_stop)> stop)
{
    Preconditioner<double> *M = feti.preconditioner;

    const char *cg = "   cg", *mixed = "cg-gp", *prop = "   gp", *opt = cg;
    info.n_cg = info.n_mixed = info.n_gproj = 0;

    math::copy(g0, b);
    math::scale(-1., g0);
    H(x0, g);
    math::add(g, 1., g0);
    updateReducedGradient(g_red, g, x0, alpha);
    updateFreeAndActiveSet(free, active, x0);
    updateStoppingGradient(g_stop, g, x0, alpha);
    updateFreeGradient(g_free, g, free);
    math::copy(x, x0); // to feasible?
    M->apply(g_free, z);
    multByFree(z, free);
    math::copy(p, z);

    while (info.iterations++ < feti.configuration.max_iterations && !stop(x, g_stop)) {
        info.time.current = eslog::time();
        // PROPORTIONALITY TEST
        if (2 * feti.configuration.delta * std::max(0., g_stop.dot(g)) <= std::max(0., g_stop.dot(g_free))) {
            H(p, Fp);
            double pFp = p.dot(Fp);
            double alpha_cg = z.dot(g) / pFp;
            double alpha_f = getFeasibleStepLength(x, p);
            if (alpha_cg <= alpha_f) { // Conjugate gradient step
                ++info.n_cg;
                opt = cg;
                math::add(x, -alpha_cg, p);
                math::add(g, -alpha_cg, Fp);
                updateReducedGradient(g_red, g, x, alpha);
                updateFreeAndActiveSet(free, active, x);
                updateFreeGradient(g_free, g, free);
                M->apply(g_free, z);
                multByFree(z, free);
                double gamma_cg = z.dot(Fp) / pFp;
                math::add(z, -gamma_cg, p);
                math::copy(p, z);
            } else { // Mixed step
                ++info.n_mixed;
                opt = mixed;
                math::copy(gg0, g0); math::add(gg0, 1., g);
                double fx = .5 * x.dot(gg0);
                math::copy(nx, x); math::add(nx, -alpha_cg, p);
                restrictToFeasibleSet(nx);
                H(nx, ng);
                math::add(ng, 1., g0);
                math::copy(ngg0, ng);
                math::add(ngg0, 1., g0);
                double fnx = .5 * nx.dot(ngg0);
                if (fnx < fx) {
                    math::copy(x, nx);
                    math::copy(g, ng);
                } else { // Halfstep
                    math::add(x, -alpha_f, p);
                    math::add(g, -alpha_f, Fp);
                    updateReducedGradient(g_red, g, x, alpha);
                    updateFreeAndActiveSet(free, active, x);
                    if (feti.configuration.gradproj) {
                        math::add(x, -alpha, g_red);
                    } else {
                        multByFree(g_red, free);
                        math::add(x, -alpha, g_red);
                    }
                    H(x, g);
                    math::add(g, 1., g0);
                }
                // Gradient projection step
                updateReducedGradient(g_red, g, x, alpha);
                updateFreeAndActiveSet(free, active, x);
                updateFreeGradient(g_free, g, free);
                M->apply(g_free, z);
                multByFree(z, free);
                math::copy(p, z);
            }
        } else {
            ++info.n_gproj;
            opt = prop;
            // Proportioning step
            if (feti.configuration.gradproj) {
                math::add(x, -alpha, g_red);
            } else {
                multByFree(g_red, active);
                math::add(x, -alpha, g_red);
            }
            H(x, g);
            math::add(g, 1., g0);
            updateReducedGradient(g_red, g, x, alpha);
            updateFreeAndActiveSet(free, active, x);
            updateFreeGradient(g_free, g, free);
            M->apply(g_free, z);
            multByFree(z, free);
            math::copy(p, z);
        }
        updateStoppingGradient(g_stop, g, x, alpha);
        if (info.print) {
            eslog::info("       - %9d    %9s             %9.4e        %9.4e      %7.2e - \n", info.iterations, opt, std::sqrt(g_stop.dot()), feti.configuration.precision * std::sqrt(feti.dualOperator->d.dot()), eslog::time() - info.time.current);
        }
    }
}

template <> void MPRGP<double>::solve(const step::Step &step, IterativeSolverInfo &info)
{
    DualOperator<double> *F = feti.dualOperator;
//    Preconditioner<double> *M = feti.preconditioner;
    // Projector<double> *P = feti.projector;

    eslog::info("       = ----------------------------------------------------------------------------- = \n");

    double maxEIG;
    int nIt;
    F->estimateMaxEigenValue(maxEIG, nIt, feti.configuration.power_precision, feti.configuration.power_maxit);
    math::set(x0, 0.);
    math::copy(b, F->d);

    eslog::info("       - ESTIMATED MAX EIGEN VALUE                             %.3e in %4d steps - \n", maxEIG, nIt);
    eslog::info("       - ----------------------------------------------------------------------------- - \n");
    eslog::info("       - ITERATION        STEP                 NORM(G)         TOLERANCE      TIME [s] - \n");
    eslog::info("       - ----------------------------------------------------------------------------- - \n");

    double alpha = feti.configuration.alpham / maxEIG;

    auto F_apply = [&] (Vector_Dual<double> &in, Vector_Dual<double> &out) {
        F->apply(in, out);
    };

    auto stop = [&] (const Vector_Dual<double> &x, const Vector_Dual<double> &g_stop) {
        return std::sqrt(g_stop.dot()) <= feti.configuration.precision * std::sqrt(feti.dualOperator->d.dot());
    };

    eslog::checkpointln("FETI: MPRGP INITIALIZATION");
    eslog::startln("MPRGP: ITERATIONS STARTED", "mprgp");

    MPRGPSolverInfo mprgp_info;
    mprgp_info.print = true;
    run(step, mprgp_info, alpha, F_apply, stop);
    info = mprgp_info;

    info.converged = stop(x, g_stop);
    if (!info.converged && feti.configuration.max_iterations <= info.iterations) {
        info.error = IterativeSolverInfo::ERROR::MAX_ITERATIONS_REACHED;
    }

    eslog::endln("mprgp: finished");
    eslog::checkpointln("FETI: MPRGP ITERATIONS");
    reconstructSolution(x, z, step);
    eslog::checkpointln("FETI: SOLUTION RECONSTRUCTION");
    eslog::info("       = ----------------------------------------------------------------------------- = \n");
}

template <> void MPRGP<std::complex<double> >::run(const step::Step &step, MPRGPSolverInfo &info, double alpha, std::function<void(Vector_Dual<std::complex<double>> &in, Vector_Dual<std::complex<double>> &out)> H, std::function<bool(const Vector_Dual<std::complex<double>> &x, const Vector_Dual<std::complex<double>> &g_stop)> stop)
{

}

template <> void MPRGP<std::complex<double> >::solve(const step::Step &step, IterativeSolverInfo &info)
{

}

template class MPRGP<double>;
template class MPRGP<std::complex<double> >;

}
