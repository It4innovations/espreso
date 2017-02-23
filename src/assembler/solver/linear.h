
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
			store::ResultStore* store,
			Matrices restriction = Matrices::NONE);

	virtual void run(Step &step);

	virtual void init(Step &step);
	virtual void preprocess(Step &step);
	virtual void solve(Step &step);
	virtual void postprocess(Step &step);
	virtual void finalize(Step &step);
};

}


#endif /* SRC_ASSEMBLER_SOLVER_LINEAR_H_ */
