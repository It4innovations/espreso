
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

	virtual void run(const Step &step) =0;

	virtual ~Solver();

	std::vector<Physics*> physics;
	std::vector<Instance*> instances;
	std::vector<LinearSolver*> linearSolvers;

protected:
	double norm(const std::vector<std::vector<double> > &v) const;

	void assembleStiffnessMatrices(const Step &step);
	void subtractResidualForces(const Step &step);
	void assembleB1(const Step &step);
	void subtractSolutionFromB1c();
	void makeStiffnessMatricesRegular();
	void assembleB0(const Step &step);
	void addToPrimar(size_t instance, const std::vector<std::vector<double> > &values);
	void storeSolution(const Step &step);

	void storeInput(std::vector<SparseMatrix> &matrices, const std::string &name, const std::string &description);
	void storeInput(std::vector<std::vector<double> > &vectors, const std::string &name, const std::string &description);

	void initLinearSolver();
	void startLinearSolver();
	void finalizeLinearSolver();

	Mesh *_mesh;
	store::ResultStore* _store;

	TimeEval *_timeStatistics;
};

}



#endif /* SRC_ASSEMBLER_SOLVER_SOLVER_H_ */
