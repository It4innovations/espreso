
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
			store::ResultStore* store,
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

	void sum(const Step &step, Matrices v1, Matrices v2, double alpha = 1, double beta = 1);
	void sum(const Step &step, Matrices v1, const std::vector<std::vector<double> > &v2, double alpha = 1, double beta = 1, const std::string v2name = "{?}");
	void sum(std::vector<std::vector<double> > &result, const std::vector<std::vector<double> > &a, const std::vector<std::vector<double> > &b, double alpha = 1, double beta = 1, const std::string resultName = "{?}", const std::string aName = "{?}", const std::string bName = "{?}");
	void multiply(const Step &step, Matrices v1, std::vector<std::vector<double> > &v2, std::vector<std::vector<double> > &solution, double beta = 1, const std::string v2name = "{?}", const std::string solutionName = "{?}");

	void composeGluing(const Step &step, Matrices matrices);
	void regularizeMatrices(const Step &step, Matrices matrices);
	void processSolution(const Step &step);

	void initLinearSolver();
	void updateLinearSolver(Matrices matrices);
	void runLinearSolver();
	void finalizeLinearSolver();

	void lineSearch(const std::vector<std::vector<double> > &U, std::vector<std::vector<double> > &deltaU, std::vector<std::vector<double> > &F_ext, Physics *physics, const Step &step);


	void storeData(const Step &step, std::vector<SparseMatrix> &matrices, const std::string &name, const std::string &description);
	void storeData(const Step &step, std::vector<std::vector<double> > &vectors, const std::string &name, const std::string &description);
	void storeSolution(const Step &step);
	void storeSubSolution(const Step &step);

	std::string _name;
	Mesh *_mesh;
	store::ResultStore* _store;
	Matrices _restriction;

	TimeEval *_timeStatistics;
};

}



#endif /* SRC_ASSEMBLER_SOLVER_SOLVER_H_ */
