
#include "iterativesolver.h"

#include "pcpg.h"
#include "orthopcpg.h"
#include "smalbe.h"
#include "mprgp.h"
#include "math/math.h"
#include "feti/projector/projector.h"
#include "feti/dualoperator/dualoperator.h"

#include "basis/utilities/sysutils.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.hpp"

#include <limits>

namespace espreso {

template <typename T>
IterativeSolver<T>* IterativeSolver<T>::create(FETI<T> &feti, const step::Step &step)
{
    if (feti.lambdas.equalities != feti.lambdas.size) {
        if (
                feti.configuration.iterative_solver != FETIConfiguration::ITERATIVE_SOLVER::SMALBE &&
                feti.configuration.iterative_solver != FETIConfiguration::ITERATIVE_SOLVER::MPRGP) {
            eslog::globalerror("FETI solver error: use SMALBE or MPRGP for solving system with inequalities.\n");
        }
    }

    switch (feti.configuration.iterative_solver) {
    case FETIConfiguration::ITERATIVE_SOLVER::PCG:           return new PCPG<T>(feti);
    case FETIConfiguration::ITERATIVE_SOLVER::orthogonalPCG: return new OrthogonalizedPCPG<T>(feti);
    case FETIConfiguration::ITERATIVE_SOLVER::SMALBE:        return new SMALBE<T>(feti);
    case FETIConfiguration::ITERATIVE_SOLVER::MPRGP:         return new MPRGP<T>(feti);
    default: return nullptr;
    }
}

template <typename T>
void IterativeSolver<T>::reconstructSolution(const Vector_Dual<T> &l, const Vector_Dual<T> &r, const step::Step &step)
{
    Projector<T> *P = feti.projector;
    DualOperator<T> *F = feti.dualOperator;

    F->toPrimal(l, iKfBtL);
    P->apply_RinvGGtG(r, Ra);
    print(step);
    #pragma omp parallel for
    for (size_t d = 0; d < feti.K.size(); ++d) {
        math::copy(feti.x[d], iKfBtL[d]);
        math::add(feti.x[d], T{1}, Ra[d]);
        math::set(feti.BtL[d], T{0});
        math::spblas::applyT(feti.BtL[d], T{1}, feti.B1[d], feti.D2C[d].data(), l);
    }
}

template <typename T>
void IterativeSolver<T>::reconstructSolution(const Vector_Dual<T> &l, const Vector_Kernel<T> &rbm, const step::Step &step)
{
    Projector<T> *P = feti.projector;
    DualOperator<T> *F = feti.dualOperator;

    F->toPrimal(l, iKfBtL);
    P->apply_R(rbm, Ra);
    print(step);
    #pragma omp parallel for
    for (size_t d = 0; d < feti.K.size(); ++d) {
        math::copy(feti.x[d], iKfBtL[d]);
        math::add(feti.x[d], T{1}, Ra[d]);
        math::set(feti.BtL[d], T{0});
        math::spblas::applyT(feti.BtL[d], T{1}, feti.B1[d], feti.D2C[d].data(), l);
    }
}

template <>
void IterativeSolver<double>::setInfo(IterativeSolverInfo &info, const FETIConfiguration &configuration, const double &ww)
{
    eslog::info("       = ----------------------------------------------------------------------------- = \n");
//    eslog::info("       - ITERATION     RELATIVE NORM      ABSOLUTE NORM      ARIOLI NORM      TIME [s] - \n");
    eslog::info("       - ITERATION     RELATIVE NORM      ABSOLUTE NORM                       TIME [s] - \n");
    eslog::info("       - ----------------------------------------------------------------------------- - \n");

    info.iterations = 1;
    info.norm.dual.absolute = info.norm.dual.initial = std::sqrt(ww);
    info.norm.dual.relative = 1;
    info.norm.dual.arioli = std::numeric_limits<double>::infinity();
    info.norm.dual.ksi = 0;
    info.converged = configuration.stopping_criterion == FETIConfiguration::STOPPING_CRITERION::ABSOLUTE && info.norm.dual.absolute < configuration.precision;
    info.time.total = info.time.current = eslog::time();

    switch (configuration.stopping_criterion) {
    case FETIConfiguration::STOPPING_CRITERION::ABSOLUTE: info.norm.dual.criteria = info.norm.dual.absolute; break;
    case FETIConfiguration::STOPPING_CRITERION::RELATIVE: info.norm.dual.criteria = info.norm.dual.relative; break;
    case FETIConfiguration::STOPPING_CRITERION::ARIOLI:   info.norm.dual.criteria = info.norm.dual.arioli; break;
    }

    if (info.converged) {
        eslog::info("       - %9d        %9.4e         %9.4e                                - \n", info.iterations, info.norm.dual.relative, info.norm.dual.absolute);
    }
    info.stagnation.buffer.resize(configuration.max_stagnation, info.norm.dual.criteria);
}

template <>
void IterativeSolver<std::complex<double> >::setInfo(IterativeSolverInfo &info, const FETIConfiguration &configuration, const std::complex<double> &ww)
{
    eslog::info("       = ----------------------------------------------------------------------------- = \n");
    eslog::info("       - ITERATION     RELATIVE NORM      ABSOLUTE NORM      ARIOLI NORM      TIME [s] - \n");
    eslog::info("       - ----------------------------------------------------------------------------- - \n");
}

template <>
void IterativeSolver<double>::updateInfo(IterativeSolverInfo &info, const FETIConfiguration &configuration, const double &ww, const double &psi, const double &ry)
{
//    info.norm.dual.ksi += psi;
    info.norm.dual.absolute = std::sqrt(std::max(0., ww));
    info.norm.dual.relative = info.norm.dual.absolute / info.norm.dual.initial;
//    info.norm.dual.arioli = std::sqrt(info.norm.dual.ksi / ry);

    switch (configuration.stopping_criterion) {
    case FETIConfiguration::STOPPING_CRITERION::ABSOLUTE: info.norm.dual.criteria = info.norm.dual.absolute; break;
    case FETIConfiguration::STOPPING_CRITERION::RELATIVE: info.norm.dual.criteria = info.norm.dual.relative; break;
    case FETIConfiguration::STOPPING_CRITERION::ARIOLI:   info.norm.dual.criteria = info.norm.dual.arioli; break;
    }
    info.converged = info.norm.dual.criteria < configuration.precision;

    info.stagnation.buffer[info.stagnation.p] = info.norm.dual.criteria;
    info.stagnation.p = (info.stagnation.p + 1) % configuration.max_stagnation;
    if (configuration.max_stagnation <= info.iterations && info.stagnation.buffer[info.stagnation.p] < info.norm.dual.criteria) {
        info.error = IterativeSolverInfo::ERROR::STAGNATION;
        info.converged = true;
    }

    if (info.iterations == feti.configuration.max_iterations && !info.converged) {
        info.error = IterativeSolverInfo::ERROR::MAX_ITERATIONS_REACHED;
        info.converged = true;
    }

    if (std::isnan(info.norm.dual.criteria)) {
        info.error = IterativeSolverInfo::ERROR::INVALID_DATA;
        info.converged = true;
    }

    if (info.converged || (info.iterations - 1) % feti.configuration.print_iteration == 0) {
        eslog::info("       - %9d        %9.4e         %9.4e                       %7.2e - \n", info.iterations, info.norm.dual.relative, info.norm.dual.absolute, eslog::time() - info.time.current);
    }
    if (info.converged) {
        info.time.total = eslog::time() - info.time.total;
    }
    if (!info.converged) {
        info.iterations++;
    }
    info.time.current = eslog::time();
}

template <>
void IterativeSolver<std::complex<double> >::updateInfo(IterativeSolverInfo &info, const FETIConfiguration &configuration, const std::complex<double> &ww, const std::complex<double> &psi, const std::complex<double> &ry)
{

}

template <typename T>
void IterativeSolver<T>::print(const step::Step &step)
{
    if (info::ecf->output.print_matrices > 1) {
        eslog::storedata(" STORE: feti/iterativesolver/{iKfBtL, Ra}\n");
        for (size_t d = 0; d < iKfBtL.size(); ++d) {
            math::store(iKfBtL[d], utils::filename(utils::debugDirectory(step) + "/feti/iterativesolver", "iKfBtL").c_str());
            math::store(Ra[d], utils::filename(utils::debugDirectory(step) + "/feti/iterativesolver", "Ra").c_str());
        }
    }
}


template class IterativeSolver<double>;
template class IterativeSolver<std::complex<double> >;

}
