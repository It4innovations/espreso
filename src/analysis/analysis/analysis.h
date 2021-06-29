
#ifndef SRC_ANALYSIS_ANALYSIS_ANALYSIS_H_
#define SRC_ANALYSIS_ANALYSIS_ANALYSIS_H_

#include "analysis/linearsolver/linearsolver.h"

namespace espreso {

struct Analysis {

template<class Solver, class Analysis, class Configuration>
static LinearSolver* init(Analysis *analysis, Configuration &configuration)
{
	Solver *solver = new Solver(configuration);
	solver->init(analysis);
	analysis->scheme.init(solver);
	analysis->mode.init(solver);
	analysis->assembler.init();
	return solver;
}

};
}

#endif /* SRC_ANALYSIS_ANALYSIS_ANALYSIS_H_ */
