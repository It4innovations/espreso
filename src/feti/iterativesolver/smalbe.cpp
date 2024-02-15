
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
: IterativeSolver<T>(feti)
{

}

template <typename T>
void SMALBE<T>::info()
{
    eslog::info(" = MODIFIED PROPORTIONING WITH REDUCED GRADIENT PROJECTION SETTING                           = \n");
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

//    if (std::sqrt(P->e.dot()) > 10 * feti.configuration.precision) {
//        x_im = Ct*(Uc\P->e);
//        b = b - F->apply(x_im);
//        constraints.sbc.lb = constraints.sbc.lb - x_im(constraints.indsbc); // indsbc -> inequalities
//        constraints.sbc.ub = constraints.sbc.ub - x_im(constraints.indsbc);
//        x0 = x0 - x_im;
//    } else {
//        x_im = zeros(n,1);
//    }

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
    // bCtmu = b - Ct*(Uc\mu)
    // x0 = 0
    // x = project_to_feasible_set(x0,constraints,options.prec);
    // g = H(x) - bCtmu
    // stopping_gradient = compute_stopping_gradient(options.stopping_gradient,g,x,constraints,alpha,options.prec);
    // Gx = Lg\(G*x);
    // free = compute_free_and_active_sets(x,constraints,options.prec);

    while (1) {
        // zavolat mprgp
        // [x,g,mpgp_outputs] = mprgp(@H,bCtmu,x,g,mpgp_options,@Hplus,constraints,C,Ct,Lc,Uc,fem,data_master);
        //
        // stopping_gradient = compute_stopping_gradient(options.stopping_gradient,g,x,constraints,alpha,options.prec);
        // Gx = Lg\(G*x);

        // mit_flag=out_of_maxiter_stop; min_flag=minimization_good_enough_stop;
        // if mit_flag || (min_flag && bch_flag); break; end;

        // 2. Update mu and M
        // Lag1 = 0.5*(x'*(g - bCtmu));
        // mu = mu + rho*Cx; bCtmu_prev=bCtmu; bCtmu = b - Ct*(Uc\mu);
        // g = g - bCtmu + bCtmu_prev;
        // if Lag1 < Lag0 + rho*norm(Gx)^2/2
          // M = options.beta*M; mpgp_options.M = M;
        // end
        // Lag0 = Lag1;
    }

    // x = x+x_im;
    // rbm = Uc\( Lc\(C*(-b_+A(x)))-mu*normPAP ); rbm = -rbm(ip);
    // reconstructSolution(x, g);
}

template <> void SMALBE<std::complex<double> >::solve(const step::Step &step, IterativeSolverInfo &info)
{

}

template class SMALBE<double>;
template class SMALBE<std::complex<double> >;

}
