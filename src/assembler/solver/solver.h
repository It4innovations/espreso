
#ifndef SRC_ASSEMBLER_SOLVER_SOLVER_H_
#define SRC_ASSEMBLER_SOLVER_SOLVER_H_

#include <cstddef>
#include <vector>
#include <string>

namespace espreso {

struct Step;
class Mesh;
class Physics;
class Instance;
class LinearSolver;
namespace store { class ResultStore; }
class TimeEval;
class SparseMatrix;


class Solver
{
public:
	Solver(
			Mesh *mesh,
			std::vector<Physics*> &physics,
			std::vector<Instance*> &instances,
			std::vector<LinearSolver*> &linearSolvers,
			store::ResultStore* store);

	virtual void run(Step &step) =0;

	virtual ~Solver();

	std::vector<Physics*> physics;
	std::vector<Instance*> instances;
	std::vector<LinearSolver*> linearSolvers;

protected:
	void assembleStiffnessMatrices(const Step &step);
	void subtractResidualForces(const Step &step);
	void assembleB1(const Step &step);
	void subtractSolutionFromB1c(const Step &step);
	void makeStiffnessMatricesRegular(const Step &step);
	void assembleB0(const Step &step);

	void postProcessDelta(Physics *physics, const std::vector<std::vector<double> > &previous);
	void addToSolution(Physics *physics, const std::vector<std::vector<double> > &previous);

	void storeSolution(const Step &step);

	void storeData(const Step &step, std::vector<SparseMatrix> &matrices, const std::string &name, const std::string &description);
	void storeData(const Step &step, std::vector<std::vector<double> > &vectors, const std::string &name, const std::string &description);

	void initLinearSolver();
	void startLinearSolver();
	void finalizeLinearSolver();

	Mesh *_mesh;
	store::ResultStore* _store;

	TimeEval *_timeStatistics;
};

}



#endif /* SRC_ASSEMBLER_SOLVER_SOLVER_H_ */
