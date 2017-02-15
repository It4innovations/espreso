
#ifndef SRC_ASSEMBLER_SOLVER_NEWTONRHAPSON_H_
#define SRC_ASSEMBLER_SOLVER_NEWTONRHAPSON_H_

#include "solver.h"

namespace espreso {

struct NonLinearSolverBase;

class NewtonRhapson: public Solver
{
public:
	NewtonRhapson(
			Mesh *mesh,
			std::vector<Physics*> &physics,
			std::vector<Instance*> &instances,
			std::vector<LinearSolver*> &linearSolvers,
			store::ResultStore* store,
			const NonLinearSolverBase &configuration);

	virtual void run(Step &step);

protected:
	const NonLinearSolverBase &_configuration;
};

}



#endif /* SRC_ASSEMBLER_SOLVER_NEWTONRHAPSON_H_ */
