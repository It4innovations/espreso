
#ifndef SRC_ANALYSIS_ANALYSIS_ANALYSIS_H_
#define SRC_ANALYSIS_ANALYSIS_ANALYSIS_H_

#include "analysis/linearsystem/empty.h"
#include "analysis/linearsystem/fetisolver.h"
#include "analysis/linearsystem/mklpdsssolver.h"
#include "analysis/linearsystem/suitesparse.h"

#include "analysis/pattern/pattern.h"

namespace espreso {

namespace step { struct Step; }

struct Physics {
    virtual ~Physics() {}

    virtual bool analyze(step::Step &step) =0;
    virtual bool run(step::Step &step, Physics *prev) =0;
};

template <typename T, typename Configuration>
LinearSystemSolver<T>* setSolver(Configuration &configuration)
{
    switch (configuration.solver) {
    case LoadStepSolverConfiguration::SOLVER::FETI:        return new FETILinearSystemSolver<T>(configuration.feti);
    case LoadStepSolverConfiguration::SOLVER::HYPRE:       break;
    case LoadStepSolverConfiguration::SOLVER::MKLPDSS:     return new MKLPDSSLinearSystemSolver<T>(configuration.mklpdss);
    case LoadStepSolverConfiguration::SOLVER::PARDISO:     break;
    case LoadStepSolverConfiguration::SOLVER::SUPERLU:     break;
    case LoadStepSolverConfiguration::SOLVER::SUITESPARSE: return new SuiteSparseLinearSystemSolver<T>(configuration.suitesparse);
    case LoadStepSolverConfiguration::SOLVER::WSMP:        break;
    case LoadStepSolverConfiguration::SOLVER::NONE:        break;
    }
    return new EmptySystemSolver<T>();
}

}

#endif /* SRC_ANALYSIS_ANALYSIS_ANALYSIS_H_ */
