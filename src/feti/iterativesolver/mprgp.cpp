
#include "mprgp.h"

#include "basis/utilities/sysutils.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.hpp"
#include "math/math.h"
#include "feti/projector/projector.h"
#include "feti/preconditioner/preconditioner.h"
#include "feti/dualoperator/dualoperator.h"

namespace espreso {

// https://doi.org/10.1007/b138610

template <typename T>
MPRGP<T>::MPRGP(FETI<T> &feti)
: IterativeSolver<T>(feti)
{
    Dual_Buffer<int>::set(feti);
}

template <typename T>
void MPRGP<T>::info()
{
    eslog::info(" = CONJUGATE PROJECTED GRADIENT SETTINGS                                                     = \n");
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
	for (esint i = feti.lambdas.equalities, j = 0; i < feti.lambdas.size; ++i, ++j) {
		x.vals[i] = std::max(feti.lb.vals[j], std::min(feti.ub.vals[j], x.vals[i]));
	}
}

template <>
void MPRGP<double>::updateFreeAndActiveSet(Vector_Dual<int> &free, Vector_Dual<int> &active, const Vector_Dual<double> &x, double epsilon)
{
    math::set(free, 1);
    math::set(active, 0);
    for (int i = feti.lambdas.equalities, j = 0; i < feti.lambdas.size; ++i, ++j) {
        free.vals[i] = (x.vals[i] - feti.lb.vals[j] > epsilon && feti.ub.vals[j] - x.vals[i] > epsilon) ? 1 : 0;
        active.vals[i] = free.vals[i] ? 0 : 1;
    }
    free.synchronize();
    active.synchronize();
}

template <>
void MPRGP<double>::updateReducedGradient(Vector_Dual<double> &g_red, Vector_Dual<double> &fi_red, Vector_Dual<double> &beta_red, const Vector_Dual<double> &g, const Vector_Dual<double> &x, double alpha, double epsilon)
{
    math::copy(xp, x);
    restrictToFeasibleSet(xp);
    math::copy(g_red, x);
    math::add(g_red, -1., xp);
    math::scale(1 / alpha, g_red);
    for (int i = feti.lambdas.equalities; i < feti.lambdas.size; ++i) {
        if (g_red.vals[i] * g.vals[i] < 0) {
            g_red.vals[i] = 0;
        }
    }
    g_red.synchronize();
    updateFreeAndActiveSet(free, active, x, epsilon);
    for (int i = 0; i < g_red.size; ++i) {
        fi_red.vals[i] = free.vals[i] ? g_red.vals[i] : 0;
        beta_red.vals[i] = active.vals[i] ? g_red.vals[i] : 0;
    }
}

template <>
void MPRGP<double>::updateStoppingGradient(Vector_Dual<double> &g_stop, Vector_Dual<double> &fi_stop, Vector_Dual<double> &beta_stop, Vector_Dual<double> &g, Vector_Dual<double> &x, double alpha, double epsilon)
{
    updateReducedGradient(g_stop, fi_stop, beta_stop, g, x, alpha, epsilon);
}

template <>
void MPRGP<double>::updateFreeGradient(Vector_Dual<double> &g_free, Vector_Dual<double> &g, Vector_Dual<int> &free, double alpha, double epsilon)
{
    for (int i = 0; i < feti.lambdas.size; ++i) {
        g_free.vals[i] = free.vals[i] ? g.vals[i] : 0;
    }
    g_free.synchronize();
}

template <>
double MPRGP<double>::getFeasibleStepLength(Vector_Dual<double> &x, Vector_Dual<double> &p)
{
    double alpha = std::numeric_limits<double>::max();
    for (int i = feti.lambdas.equalities, j = 0; i < feti.lambdas.size; ++i, ++j) {
        if (p.vals[i] > 0) {
            alpha = std::min(alpha, x.vals[i] - feti.lb.vals[j] / p.vals[i]);
        }
        if (p.vals[i] < 0) {
            alpha = std::min(alpha, x.vals[i] - feti.ub.vals[j] / p.vals[i]);
        }
    }
    return alpha;
}

template <>
bool MPRGP<double>::stop(const Vector_Dual<double> &x, double epsilon)
{
    // if (it == maxit) { return true; }
    return std::sqrt(x.dot()) <= epsilon * std::sqrt(feti.dualOperator->d.dot());
}

template <> void MPRGP<double>::solve(const step::Step &step, IterativeSolverInfo &info)
{
    DualOperator<double> *F = feti.dualOperator;
    Preconditioner<double> *M = feti.preconditioner;
    // Projector<double> *P = feti.projector;

    double maxEIG;
    int nIt;
    F->estimateMaxEigenValue(0.001, 10, maxEIG, nIt);
    double alpha = 1.95 / maxEIG;
    double delta = 0.25;
    double epsilon = 1e-10;
    bool gradproj = true;
    math::set(x0, 0.);

    math::copy(g0, F->d);
    math::scale(-1., g0);
    F->apply(x0, g);
    math::add(g0, 1., g);
    updateReducedGradient(g_red, fi_red, beta_red, g, x0, alpha, epsilon);
    updateStoppingGradient(g_stop, fi_stop, beta_stop, g, x0, alpha, epsilon);
    updateFreeGradient(g_free, g, free, alpha, epsilon);
    math::copy(x, x0); // to feasible?
    math::copy(z, g_free);
    M->apply(g_free, z);
    // z * free -> z -->> copy of replace by zero
    math::copy(p, z);

    while (!stop(g_stop, epsilon)) {
        // PROPORTIONALITY TEST
        if (2 * delta * std::max(0., g.dot(g_stop)) <= std::max(0., g_stop.dot(g_free))) {
            F->apply(p, Fp);
            double pFp = p.dot(Fp);
            double alpha_cg = p.dot() / pFp;
            double alpha_f = getFeasibleStepLength(x, p);
            if (alpha_cg <= alpha_f) { // Conjugate gradient step
                math::add(x, -alpha_cg, p);
                math::add(g, -alpha_cg, Fp);
                updateReducedGradient(g_red, fi_red, beta_red, g, x, alpha, epsilon);
                M->apply(g_free, z);
                // z * free -> z -->> copy of replace by zero
                double gamma_cg = z.dot(Fp) / pFp;
                math::copy(p, z);
                math::add(p, -gamma_cg, p);
            } else { // Mixed step
                math::copy(gg0, g0); math::add(gg0, 1., g);
                double fx = .5 * x.dot(gg0);
                math::copy(nx, x); math::add(nx, -alpha_cg, p);
                restrictToFeasibleSet(nx);
                F->apply(nx, ng);
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
                    updateReducedGradient(g_red, fi_red, beta_red, g, x, alpha, epsilon);
                    if (gradproj) {
                        math::add(x, -alpha, g_red);
                    } else {
                        math::add(x, -alpha, g_red); // g_red.* free
                    }
                    F->apply(x, g);
                    math::add(g, 1., g0);
                }
                // Gradient projection step
                updateReducedGradient(g_red, fi_red, beta_red, g, x, alpha, epsilon);
                M->apply(g_free, z);
                // z * free -> z -->> copy of replace by zero
                math::copy(p, z);
            }
        } else {
            // Proportioning step
            if (gradproj) {
                math::add(x, -alpha, g_red);
            } else {
                math::add(x, -alpha, g_red); // g_red.* free
            }
            F->apply(x, g);
            math::add(g, 1., g0);
            updateReducedGradient(g_red, fi_red, beta_red, g, x, alpha, epsilon);
            M->apply(g_free, z);
            // z * free -> z -->> copy of replace by zero
            math::copy(p, z);
        }
        updateStoppingGradient(g_stop, fi_stop, beta_stop, g, x, alpha, epsilon);
    }
}

template <> void MPRGP<std::complex<double> >::solve(const step::Step &step, IterativeSolverInfo &info)
{

}

template class MPRGP<double>;
template class MPRGP<std::complex<double> >;

}
