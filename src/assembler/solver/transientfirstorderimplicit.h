
#ifndef SRC_ASSEMBLER_SOLVER_TRANSIENTFIRSTORDERIMPLICIT_H_
#define SRC_ASSEMBLER_SOLVER_TRANSIENTFIRSTORDERIMPLICIT_H_

#include "solver.h"

namespace espreso {

struct TransientSolver;

class TransientFirstOrderImplicit: public Solver
{
public:
	TransientFirstOrderImplicit(
			Mesh *mesh,
			Physics* physics,
			LinearSolver* linearSolver,
			store::ResultStore* store,
			const TransientSolver &configuration);

	virtual void run(Step &step);

	virtual void init(Step &step);
	virtual void preprocess(Step &step);
	virtual void solve(Step &step);
	virtual void postprocess(Step &step);
	virtual void finalize(Step &step);

protected:
	const TransientSolver &_configuration;
};

}



#endif /* SRC_ASSEMBLER_SOLVER_TRANSIENTFIRSTORDERIMPLICIT_H_ */
