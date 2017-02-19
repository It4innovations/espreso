
#ifndef SRC_ASSEMBLER_SOLVER_SOLVER_H_
#define SRC_ASSEMBLER_SOLVER_SOLVER_H_

#include <cstddef>
#include <vector>
#include <string>

#include "../instance.h"

namespace espreso {

struct Step;
class Mesh;
class Physics;
class LinearSolver;
namespace store { class ResultStore; }
class TimeEval;
class SparseMatrix;


class Solver
{
public:
	Solver(
			const std::string &name,
			Mesh *mesh,
			Physics* physics,
			LinearSolver* linearSolver,
			store::ResultStore* store);

	virtual void run(Step &step) =0;
	const std::string& name() const { return _name; }

	virtual ~Solver();

	Physics* physics;
	LinearSolver* linearSolver;

protected:
	void assembleMatrices(const Step &step, Matrices matrices);
	void updateMatrices(const Step &step, Matrices matrices, const std::vector<Solution*> &solution);

	void composeGluing(const Step &step, Matrices matrices);
	void regularizeMatrices(const Step &step, Matrices matrices);
	void processSolution(const Step &step);

	void initLinearSolver();
	void updateLinearSolver(Matrices matrices);
	void runLinearSolver();
	void finalizeLinearSolver();

	void subtractSolutionFromB1c(const Step &step);


	void lineSearch(const std::vector<std::vector<double> > &U, std::vector<std::vector<double> > &deltaU, std::vector<std::vector<double> > &F_ext, Physics *physics, const Step &step);
	void sumVectors(std::vector<std::vector<double> > &result, const std::vector<std::vector<double> > &a, const std::vector<std::vector<double> > &b, double alpha = 1, double beta = 1);

	void storeData(const Step &step, std::vector<SparseMatrix> &matrices, const std::string &name, const std::string &description);
	void storeData(const Step &step, std::vector<std::vector<double> > &vectors, const std::string &name, const std::string &description);
	void storeSolution(const Step &step);

	std::string _name;
	Mesh *_mesh;
	store::ResultStore* _store;

	TimeEval *_timeStatistics;
};

}



#endif /* SRC_ASSEMBLER_SOLVER_SOLVER_H_ */
