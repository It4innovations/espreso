
#ifndef SRC_ASSEMBLER_SOLVER_LINEAR_H_
#define SRC_ASSEMBLER_SOLVER_LINEAR_H_

#include "solver.h"

namespace espreso {

class Linear: public Solver
{
public:
	Linear(
			Mesh *mesh,
			Physics* physics,
			LinearSolver* linearSolver,
			store::ResultStore* store);

	virtual void run(Step &step);
};

}


#endif /* SRC_ASSEMBLER_SOLVER_LINEAR_H_ */
