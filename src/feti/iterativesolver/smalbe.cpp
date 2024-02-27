
#include "smalbe.h"

#include "basis/utilities/sysutils.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.hpp"
#include "math/math.h"
#include "feti/projector/projector.h"
#include "feti/dualoperator/dualoperator.h"

namespace espreso {

template <typename T>
SMALBE<T>::SMALBE(FETI<T> &feti)
: IterativeSolver<T>(feti), mprgp(feti)
{

}

template <typename T>
void SMALBE<T>::info()
{
    eslog::info(" = MODIFIED PROPORTIONING WITH REDUCED GRADIENT PROJECTION SETTING                           = \n");
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
        eslog::storedata(" STORE: feti/SMALBE/{%s%s}\n", name, std::to_string(info.iterations).c_str());
        math::store(v, utils::filename(utils::debugDirectory(step) + "/feti/SMALBE", std::string(name) + std::to_string(info.iterations)).c_str());
    }
}

template <> void SMALBE<double>::solve(const step::Step &step, IterativeSolverInfo &info)
{
    DualOperator<double> *F = feti.dualOperator;
    Projector<double> *P = feti.projector;

    math::copy(b, F->d);
    // Homogenization of the equality constraints

    math::set(mprgp.x0, 0.);
    P->apply_invL(P->e, invLce);
    if (std::sqrt(invLce.dot()) > 10 * feti.configuration.precision) {
        P->apply_GtinvU(invLce, x_im);
        F->apply(x_im, Fx_im);
        math::add(b, -1., Fx_im);
        for (esint i = feti.lambdas.equalities, j = 0; i < feti.lambdas.size; ++i, ++j) {
            feti.lb.vals[j] -= x_im.vals[i];
            feti.ub.vals[j] -= x_im.vals[i];
        }
        math::add(mprgp.x0, -1., x_im);
    } else {
        math::set(x_im, 1.);
    }

    // constraints=initialize_slip_bound(); coulomb, later

    P->apply(b, Pb);
    math::copy(b, Pb);

    // Initialization

    int nIt;
    double rho = 0, normPFP = 1, M_const = 1, maxEIG_H;
    F->estimateMaxProjectedEigenValue(maxEIG_H, nIt, 0.0001, 10);
    eslog::info("       - ESTIMATED MAX PROJECTED EIGEN VALUE                   %.3e in %4d steps - \n", maxEIG_H, nIt);
    eslog::info("       - ----------------------------------------------------------------------------- - \n");
    eslog::info("       - ITERATION        STEP                 NORM(G)         TOLERANCE      TIME [s] - \n");
    eslog::info("       - ----------------------------------------------------------------------------- - \n");

    // switch PFP-scaling, TODO
    normPFP = maxEIG_H;
    maxEIG_H = 1;
    math::scale(1 / normPFP, b);

    M_const *= 100; // based on experience of O. Vlach
    double alpha = feti.configuration.alpham / maxEIG_H;
    rho = feti.configuration.rho * maxEIG_H;

    double norm_b = std::sqrt(b.dot());

    math::set(mu, 0.0);
    P->apply_GtinvU(mu, bCtmu);
    math::scale(-1., bCtmu);
    math::add(bCtmu, 1., b);

    math::copy(mprgp.x, mprgp.x0);
    mprgp.restrictToFeasibleSet(mprgp.x);

    auto A_apply = [&] (Vector_Dual<double> &in, Vector_Dual<double> &out) {
        P->apply(in, y);
        F->apply(y, z);
        P->apply(z, out);
        math::scale(1. / normPFP, out);
        math::add(out,  rho, in);
        math::add(out, -rho, y);
    };

    auto stop = [&] (const Vector_Dual<double> &x, const Vector_Dual<double> &g_stop) {
        double norm_g_stop = std::sqrt(g_stop.dot());

        P->apply_invLG(x, Gx);
        double norm_Gx = std::sqrt(Gx.dot());

        return
             norm_g_stop <= std::min(M_const * norm_Gx, feti.configuration.eta * norm_b) ||
            (norm_g_stop <= feti.configuration.precision_in * norm_b &&
            norm_Gx      <= feti.configuration.precision_in * norm_b / maxEIG_H);
    };

    // Hessian multiplication (PAP+rho*Ct*inv(C*Ct)*C)*x
    A_apply(mprgp.x, mprgp.g);
    math::add(mprgp.g, -1., bCtmu);
    mprgp.updateStoppingGradient(mprgp.g_stop, mprgp.g, mprgp.x, alpha);
    mprgp.updateFreeAndActiveSet(mprgp.free, mprgp.active, mprgp.x);

    double Lag0 = -std::numeric_limits<double>::infinity();
    P->apply_invLG(mprgp.x, Gx);
    for (size_t it = 1; it <= feti.configuration.max_iterations; ++it) {
        math::copy(mprgp.b, bCtmu);
        mprgp.run(step, info, alpha, A_apply, stop);
        math::copy(mprgp.x0, mprgp.x);
        mprgp.updateStoppingGradient(mprgp.g_stop, mprgp.g, mprgp.x, alpha);
        eslog::info("       - ----------------------------------------------------------------------------- - \n");
        P->apply_invLG(mprgp.x, Gx);

        double norm_stop = std::sqrt(mprgp.g_stop.dot());
        double norm_Gx = std::sqrt(Gx.dot());
        if (norm_stop <= feti.configuration.precision * norm_b && norm_Gx <= feti.configuration.precision * norm_b / maxEIG_H) {
            break;
        }

        // 2. Update mu and M
        math::copy(gbCtmu, mprgp.g);
        math::add(gbCtmu, -1., bCtmu);
        double Lag1 = .5 * gbCtmu.dot(mprgp.x); // Lag1 = 0.5*(x'*(g - bCtmu));
        math::add(mu, rho, Gx); // mu = mu + rho*Gx;
        math::copy(bCtmu_prev, bCtmu);

        // bCtmu = b - Ct*(Uc\mu);
        P->apply_GtinvU(mu, bCtmu);
        math::scale(-1., bCtmu);
        math::add(bCtmu, 1., b);

        // g = g - bCtmu + bCtmu_prev;
        math::add(mprgp.g, -1., bCtmu);
        math::add(mprgp.g,  1., bCtmu_prev);

        if (Lag1 < Lag0 + rho * Gx.dot() / 2) {
            M_const = feti.configuration.beta * M_const;
        }
        Lag0 = Lag1;
    }

    math::add(mprgp.x, 1., x_im);
    // rbm = Uc\( Lc\(C*(-b_+A(x)))-mu*normPAP ); rbm = -rbm(ip);
    F->apply(mprgp.x, b);
    math::add(b, -1., F->d);
    P->apply_invLG(b, Gx);
    math::add(Gx, -normPFP, mu);
    P->apply_invU(Gx, mu); // rbm = mu
    math::scale(-1., mu);
    // permute -mu

    reconstructSolution(mprgp.x, mu);
}

template <> void SMALBE<std::complex<double> >::solve(const step::Step &step, IterativeSolverInfo &info)
{

}

template class SMALBE<double>;
template class SMALBE<std::complex<double> >;

}
