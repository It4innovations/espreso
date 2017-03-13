
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
namespace output { class Store; }
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
			output::Store* store,
			Matrices restriction);

	virtual void run(Step &step) =0;

	virtual void init(Step &step) =0;
	virtual void preprocess(Step &step) =0;
	virtual void solve(Step &step) =0;
	virtual void postprocess(Step &step) =0;
	virtual void finalize(Step &step) =0;

	const std::string& name() const { return _name; }

	virtual ~Solver();

	Physics* physics;
	Instance* instance;
	LinearSolver* linearSolver;

protected:
	void assembleMatrices(const Step &step, Matrices matrices);
	void updateMatrices(const Step &step, Matrices matrices, const std::vector<Solution*> &solution);

	/// z = a * x + b + y
	void sum(std::vector<std::vector<double> > &z, double a, const std::vector<std::vector<double> > &x, double b, const std::vector<std::vector<double> > &y, const std::string &description);
	/// z = a * x + b + y (prefix variant)
	void sum(std::vector<std::vector<double> > &z, double a, const std::vector<std::vector<double> > &x, double b, const std::vector<std::vector<double> > &y, const std::vector<size_t> &prefix, const std::string &description);
	/// A += beta * B
	void sum(std::vector<SparseMatrix> &A, double beta, std::vector<SparseMatrix> &B, const std::string &description);

	void subtractDirichlet();

	void multiply(const Step &step, Matrices v1, std::vector<std::vector<double> > &v2, std::vector<std::vector<double> > &solution, double beta = 1, const std::string v2name = "{?}", const std::string solutionName = "{?}");
	void multiply(const Step &step, Matrices v, double beta);

	void composeGluing(const Step &step, Matrices matrices);
	void regularizeMatrices(const Step &step, Matrices matrices);
	void processSolution(const Step &step);

	void initLinearSolver(const Step &step);
	void updateLinearSolver(const Step &step, Matrices matrices);
	void runLinearSolver(const Step &step);
	void finalizeLinearSolver(const Step &step);

	double lineSearch(const std::vector<std::vector<double> > &U, std::vector<std::vector<double> > &deltaU, std::vector<std::vector<double> > &F_ext, Physics *physics, const Step &step);

	double maxAbsValue(const std::vector<std::vector<double> > &v) const;

	void storeData(const Step &step, std::vector<SparseMatrix> &matrices, const std::string &name, const std::string &description);
	void storeData(const Step &step, std::vector<std::vector<double> > &vectors, const std::string &name, const std::string &description);
	void storeSolution(const Step &step);
	void storeSubSolution(const Step &step);

	std::string _name;
	Mesh *_mesh;
	output::Store* _store;
	Matrices _restriction;

	TimeEval *_timeStatistics;
};

}



#endif /* SRC_ASSEMBLER_SOLVER_SOLVER_H_ */
