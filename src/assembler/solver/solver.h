
#ifndef SRC_ASSEMBLER_SOLVER_SOLVER_H_
#define SRC_ASSEMBLER_SOLVER_SOLVER_H_

#include <vector>

namespace espreso {

struct Step;
class Mesh;
class Physics;
class Instance;
class LinearSolver;
namespace store { class ResultStore; }
class TimeEval;


class Solver
{
public:
	Solver(
			Mesh *mesh,
			std::vector<Physics*> &physics,
			std::vector<Instance*> &instances,
			std::vector<LinearSolver*> &linearSolvers,
			store::ResultStore* store);

	virtual void run(const Step &step) =0;

	virtual ~Solver();

	std::vector<Physics*> physics;
	std::vector<Instance*> instances;
	std::vector<LinearSolver*> linearSolvers;

protected:
	void assembleStiffnessMatrices(const Step &step);
	void assembleB1(const Step &step);
	void makeStiffnessMatricesRegular();
	void assembleB0(const Step &step);

	void initLinearSolver();
	void startLinearSolver(const Step &step);

	void finalizeLinearSolver();

	Mesh *_mesh;
	store::ResultStore* _store;

	TimeEval *_timeStatistics;
};

}



#endif /* SRC_ASSEMBLER_SOLVER_SOLVER_H_ */
