
#ifndef SRC_ASSEMBLER_SOLVER_SOLVER_H_
#define SRC_ASSEMBLER_SOLVER_SOLVER_H_

#include <vector>

namespace espreso {

class Mesh;
class NewPhysics;
class NewInstance;
class LinearSolver;
namespace store { class ResultStore; }
class TimeEval;


class Solver
{
public:
	Solver(
			Mesh *mesh,
			std::vector<NewPhysics*> &physics,
			std::vector<NewInstance*> &instances,
			std::vector<LinearSolver*> &linearSolvers,
			store::ResultStore* store);

	virtual void init() =0;
	virtual void solve(std::vector<std::vector<double> > &solution) =0;
	virtual void finalize() =0;

	virtual ~Solver();

protected:
	void meshPreprocessing();
	void assembleStiffnessMatrices();
	void assembleB1();
	void makeStiffnessMatricesRegular();
	void assembleB0();

	void initLinearSolver();
	void startLinearSolver(std::vector<std::vector<double> > &solution);

	void finalizeLinearSolver();

	Mesh *_mesh;
	std::vector<NewPhysics*> &_physics;
	std::vector<NewInstance*> &_instances;
	std::vector<LinearSolver*> &_linearSolvers;
	store::ResultStore* _store;

	TimeEval *_timeStatistics;
};

}



#endif /* SRC_ASSEMBLER_SOLVER_SOLVER_H_ */
