
#ifndef SRC_ASSEMBLER_SOLVER_TRANSIENTFIRSTORDERIMPLICIT_H_
#define SRC_ASSEMBLER_SOLVER_TRANSIENTFIRSTORDERIMPLICIT_H_

#include "solver.h"

namespace espreso {

class TransientFirstOrderImplicit: public Solver
{
public:
	TransientFirstOrderImplicit(
			Mesh *mesh,
			Physics* physics,
			LinearSolver* linearSolver,
			store::ResultStore* store);

	virtual void run(Step &step);

	virtual void init(Step &step);
	virtual void preprocess(Step &step);
	virtual void solve(Step &step);
	virtual void postprocess(Step &step);
	virtual void finalize(Step &step);
};

}



#endif /* SRC_ASSEMBLER_SOLVER_TRANSIENTFIRSTORDERIMPLICIT_H_ */
