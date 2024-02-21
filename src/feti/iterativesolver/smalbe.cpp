
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

// function y=H(x)
// % Hessian multiplication (PAP+rho*Ct*inv(C*Ct)*C)*x
// v = P(x);
// y = A(v);
// y = P(y);
// y = y/normPAP; % Scale PAP
// y = y + rho*(x-v);
// end

//function flag=minimization_good_enough_stop(~)
//  norm_b = norm(b);
//  flag = (norm(stopping_gradient)<=options.epsilon*norm_b && norm(Cx)<=options.epsilon*norm_b/maxeig);
//end



template <> void SMALBE<double>::solve(const step::Step &step, IterativeSolverInfo &info)
{
    DualOperator<double> *F = feti.dualOperator;
    Preconditioner<double> *M = feti.preconditioner;
    Projector<double> *P = feti.projector;

    math::copy(b, F->d);
    math::copy(b_, F->d);
    // Homogenization of the equality constraints

    math::set(mprgp.x0, 0.);
    if (std::sqrt(P->e.dot()) > 10 * feti.configuration.precision) {
        P->apply_GtinvU(P->e, x_im);
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

    // Initialization

    int nIt;
    double rho = 0, normPFP = 1, M_const = 1, maxEIG_H;
    F->estimateMaxProjectedEigenValue(maxEIG_H, nIt, 0.001, 10);
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
    math::set(mprgp.x0, 0.);

    math::copy(mprgp.x, mprgp.x0);
    mprgp.restrictToFeasibleSet(mprgp.x);

    // Hessian multiplication (PAP+rho*Ct*inv(C*Ct)*C)*x
    P->apply(mprgp.x, y);
    F->apply(y, z);
    P->apply(z, mprgp.g);
    math::scale(1. / normPFP, mprgp.g);
    math::add(mprgp.g,  rho, mprgp.x);
    math::add(mprgp.g, -rho, z);
    math::add(mprgp.g, -1., bCtmu);
    mprgp.updateStoppingGradient(mprgp.g_stop, mprgp.g, mprgp.x, alpha, feti.configuration.precision);
    // Gx = Lg\(G*x);
    mprgp.updateFreeAndActiveSet(mprgp.free, mprgp.active, mprgp.x, feti.configuration.precision);

    for (size_t it = 1; it <= feti.configuration.max_iterations; ++it) {
        mprgp.solve(step, info);
        mprgp.updateStoppingGradient(mprgp.g_stop, mprgp.g, mprgp.x, alpha, feti.configuration.precision);
        // Gx = Lg\(G*x);

        double norm_b = std::sqrt(b.dot());
        double norm_stop = std::sqrt(mprgp.g_stop.dot());
        double norm_Cx = 0;
        if (norm_stop <= feti.configuration.precision * norm_b && norm_Cx <= feti.configuration.precision * norm_b / maxEIG_H) {
            break;
        }

        // 2. Update mu and M
        math::copy(gbCtmu, mprgp.g);
        math::add(gbCtmu, -1., bCtmu);
        double Lag1 = gbCtmu.dot(mprgp.x); // Lag1 = 0.5*(x'*(g - bCtmu));
        math::add(mu, rho, Cx); // mu = mu + rho*Cx;
        math::copy(bCtmu_prev, bCtmu);
        // bCtmu = b - Ct*(Uc\mu);
        P->apply_GtinvU(mu, bCtmu);
        math::scale(-1., bCtmu);
        math::add(bCtmu, 1., b);
        // g = g - bCtmu + bCtmu_prev;
        math::add(mprgp.g, -1., bCtmu);
        math::add(mprgp.g, -1., bCtmu_prev);

        // if Lag1 < Lag0 + rho*norm(Gx)^2/2
          // M = options.beta*M; mpgp_options.M = M;
        // end
        // Lag0 = Lag1;
    }

    // x = x+x_im;
    math::add(mprgp.x, 1., x_im);
    // rbm = Uc\( Lc\(C*(-b_+A(x)))-mu*normPAP ); rbm = -rbm(ip);
    reconstructSolution(mprgp.x, mprgp.g);
}

template <> void SMALBE<std::complex<double> >::solve(const step::Step &step, IterativeSolverInfo &info)
{

}

template class SMALBE<double>;
template class SMALBE<std::complex<double> >;

}
