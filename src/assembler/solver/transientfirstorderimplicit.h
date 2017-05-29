
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
			output::Store* store,
			const TransientSolver &configuration,
			double duration);

	virtual void run(Step &step);

	virtual void init(Step &step);
	virtual void preprocess(Step &step);
	virtual void solve(Step &step);
	virtual void postprocess(Step &step);
	virtual void finalize(Step &step);

protected:
	static size_t offset;
	static size_t lastStep;

	enum SolutionIndex: size_t {
		SOLUTION   = 0,
		DELTA      = 1,
		VELOCITY   = 2,
		X          = 3,
		Y          = 4,

		SIZE       = 3
	};

	const TransientSolver &_configuration;
};

}



#endif /* SRC_ASSEMBLER_SOLVER_TRANSIENTFIRSTORDERIMPLICIT_H_ */
