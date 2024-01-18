
#include "iterativesolver.h"

#include "cpg.h"
#include "pcpg.h"
#include "orthocpg.h"
#include "feti/projector/projector.h"
#include "feti/dualoperator/dualoperator.h"

#include "feti/common/vector_kernel.h"

#include "esinfo/eslog.hpp"

#include <limits>

namespace espreso {

template <typename T>
IterativeSolver<T>* IterativeSolver<T>::set(FETI<T> &feti, const step::Step &step)
{
	switch (feti.configuration.iterative_solver) {
	case FETIConfiguration::ITERATIVE_SOLVER::PCG:
		if (feti.configuration.preconditioner == FETIConfiguration::PRECONDITIONER::NONE) {
			eslog::info(" = ITERATIVE SOLVER                                             CONJUGATE PROJECTED GRADIENT = \n");
			return new CPG<T>(feti);
		} else {
			eslog::info(" = ITERATIVE SOLVER                              PRECONDITIONED CONJUGATE PROJECTED GRADIENT = \n");
			return new PCPG<T>(feti);
		}
	case FETIConfiguration::ITERATIVE_SOLVER::pipePCG:
	case FETIConfiguration::ITERATIVE_SOLVER::orthogonalPCG:
		if (feti.configuration.preconditioner == FETIConfiguration::PRECONDITIONER::NONE) {
			eslog::info(" = ITERATIVE SOLVER                      CONJUGATE PROJECTED GRADIENT WITH ORTHOGONALIZATION = \n");
			return new OrthogonalizedCPG<T>(feti);
		} else {
			eslog::info(" = ITERATIVE SOLVER       PRECONDITIONED CONJUGATE PROJECTED GRADIENT WITH ORTHOGONALIZATION = \n");
//			return new OrthogonalizedPCPG<T>(feti);
			return new OrthogonalizedCPG<T>(feti);
		}
	case FETIConfiguration::ITERATIVE_SOLVER::GMRES:
	case FETIConfiguration::ITERATIVE_SOLVER::BICGSTAB:
	case FETIConfiguration::ITERATIVE_SOLVER::QPCE:
	case FETIConfiguration::ITERATIVE_SOLVER::orthogonalPCG_CP:
	case FETIConfiguration::ITERATIVE_SOLVER::PCG_CP:
	default: return nullptr;
	}
}

template <typename T>
void IterativeSolver<T>::reconstructSolution(const Vector_Dual<T> &l, const Vector_Dual<T> &r)
{
	Projector<T> *P = feti.projector;
	DualOperator<T> *F = feti.dualOperator;

	F->toPrimal(l, iKfBtL);
	P->applyRInvGGtG(r, Ra);
	#pragma omp parallel for
	for (size_t d = 0; d < feti.K.size(); ++d) {
		math::copy(feti.x[d], iKfBtL.domains[d]);
		math::add(feti.x[d], T{1}, Ra.domains[d]);
	}
}

template <>
void IterativeSolver<double>::setInfo(IterativeSolverInfo &info, const FETIConfiguration &configuration, const double &ww)
{
	eslog::info("       = ----------------------------------------------------------------------------- = \n");
	eslog::info("       - ITERATION     RELATIVE NORM      ABSOLUTE NORM      ARIOLI NORM      TIME [s] - \n");
	eslog::info("       - ----------------------------------------------------------------------------- - \n");

	info.iterations = 1;
	info.norm.dual.absolute = info.norm.dual.initial = std::sqrt(ww);
	info.norm.dual.relative = 1;
	info.norm.dual.arioli = std::numeric_limits<double>::infinity();
	info.norm.dual.ksi = 0;
	info.converged = info.norm.dual.absolute < configuration.precision;
	info.time.total = info.time.current = eslog::time();

	switch (configuration.stopping_criterion) {
	case FETIConfiguration::STOPPING_CRITERION::ABSOLUTE: info.norm.dual.criteria = info.norm.dual.absolute; break;
	case FETIConfiguration::STOPPING_CRITERION::RELATIVE: info.norm.dual.criteria = info.norm.dual.relative; break;
	case FETIConfiguration::STOPPING_CRITERION::ARIOLI:   info.norm.dual.criteria = info.norm.dual.arioli; break;
	}

	if (info.converged) {
		eslog::info("       - %9d        %9.4e         %9.4e        %9.4e               - \n", info.iterations, info.norm.dual.relative, info.norm.dual.absolute, info.norm.dual.arioli);
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
//	info.norm.dual.ksi += psi;
	info.norm.dual.absolute = std::sqrt(ww);
	info.norm.dual.relative = info.norm.dual.absolute / info.norm.dual.initial;
//	info.norm.dual.arioli = std::sqrt(info.norm.dual.ksi / ry);

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
		eslog::info("       - %9d        %9.4e         %9.4e        %9.4e      %7.2e - \n", info.iterations, info.norm.dual.relative, info.norm.dual.absolute, info.norm.dual.arioli, eslog::time() - info.time.current);
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

template class IterativeSolver<double>;
template class IterativeSolver<std::complex<double> >;

}
