
#ifndef SRC_ASSEMBLER_SOLVER_LINEAR_H_
#define SRC_ASSEMBLER_SOLVER_LINEAR_H_

#include "solver.h"

namespace espreso {

class Linear: public Solver
{
public:
	Linear(
			Mesh *mesh,
			std::vector<NewPhysics*> &physics,
			std::vector<NewInstance*> &instances,
			std::vector<LinearSolver*> &linearSolvers,
			store::ResultStore* store);

	virtual void init();
	virtual void solve(std::vector<std::vector<double> > &solution);
	virtual void finalize();
};

}


#endif /* SRC_ASSEMBLER_SOLVER_LINEAR_H_ */
