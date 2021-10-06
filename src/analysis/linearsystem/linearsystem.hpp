
#ifndef SRC_ANALYSIS_LINEARSYSTEM_LINEARSYSTEM_HPP_
#define SRC_ANALYSIS_LINEARSYSTEM_LINEARSYSTEM_HPP_

#include "mklpdsssystem.h"
#include "fetisystem.h"

namespace espreso {

template <typename Assembler, typename Solver, typename Analysis>
void initSystem(AX_LinearSystem<Assembler, Solver>* &system, Analysis *analysis)
{
	switch (analysis->configuration.solver) {
	case LoadStepSolverConfiguration::SOLVER::FETI:    system = new AX_FETISystem<Analysis>(analysis); break;
	case LoadStepSolverConfiguration::SOLVER::HYPRE:   break;
	case LoadStepSolverConfiguration::SOLVER::MKLPDSS: system = new AX_MKLPDSSSystem<Analysis>(analysis); break;
	case LoadStepSolverConfiguration::SOLVER::PARDISO: break;
	case LoadStepSolverConfiguration::SOLVER::SUPERLU: break;
	case LoadStepSolverConfiguration::SOLVER::WSMP:    break;
	}
}


}

#endif /* SRC_ANALYSIS_LINEARSYSTEM_LINEARSYSTEM_HPP_ */
