
#ifndef SRC_ASSEMBLER_SOLVER_SOLVER_H_
#define SRC_ASSEMBLER_SOLVER_SOLVER_H_

#include <cstddef>
#include <vector>
#include <string>

#include "solverbase.h"
#include "../instance.h"

namespace espreso {

class Physics;
class FETISolver;
namespace output { class Store; }
class SparseMatrix;


class Solver: public SolverBase
{
public:
	Solver(
			const std::string &name,
			Mesh *mesh,
			Physics* physics,
			FETISolver* linearSolver,
			output::Store* store,
			double duration,
			Matrices restriction);

	virtual void run(Step &step) =0;

	virtual void init(Step &step) =0;
	virtual void preprocess(Step &step) =0;
	virtual void solve(Step &step) =0;
	virtual void postprocess(Step &step) =0;
	virtual void finalize(Step &step) =0;

	virtual ~Solver();

	Physics* physics;
	Instance* instance;
	FETISolver* linearSolver;

protected:
	void preprocessData(const Step &step);
	void updateMatrices(const Step &step, Matrices matrices);
	void updateMatrices(const Step &step, Matrices matrices, const std::vector<Solution*> &solution);

	/// z = a * x + b + y
	void sum(std::vector<std::vector<double> > &z, double a, const std::vector<std::vector<double> > &x, double b, const std::vector<std::vector<double> > &y, const std::string &description);
	/// z = a * x + b + y (prefix variant)
	void sum(std::vector<std::vector<double> > &z, double a, const std::vector<std::vector<double> > &x, double b, const std::vector<std::vector<double> > &y, const std::vector<size_t> &prefix, const std::string &description);
	/// A += beta * B
	void sum(std::vector<SparseMatrix> &A, double beta, std::vector<SparseMatrix> &B, const std::string &description);

	void subtractDirichlet();

	/// y = A * x
	void multiply(std::vector<std::vector<double> > &y, std::vector<SparseMatrix> &A, std::vector<std::vector<double> > &x, const std::string &description);
	/// x = a * x
	void multiply(std::vector<std::vector<double> > &x, double a, const std::string &description);

	void composeGluing(const Step &step, Matrices matrices);
	void regularizeMatrices(const Step &step, Matrices matrices);
	void setEmptyRegularization(const Step &step, Matrices matrices);
	void processSolution(const Step &step);

	void initLinearSolver(const Step &step);
	void updateLinearSolver(const Step &step, Matrices matrices);
	void runLinearSolver(const Step &step);
	void finalizeLinearSolver(const Step &step);

	double lineSearch(const std::vector<std::vector<double> > &U, std::vector<std::vector<double> > &deltaU, std::vector<std::vector<double> > &F_ext, Physics *physics, const Step &step);

	double maxAbsValue(const std::vector<std::vector<double> > &v) const;

	void storeData(const Step &step, SparseMatrix &matrix, size_t domain, const std::string &name, const std::string &description);
	void storeData(const Step &step, std::vector<SparseMatrix> &matrices, const std::string &name, const std::string &description);
	void storeData(const Step &step, std::vector<std::vector<double> > &vectors, const std::string &name, const std::string &description);
	void storeSolution(const Step &step);
	void storeSubSolution(const Step &step);

	output::Store* _store;
	Matrices _restriction;
};

}



#endif /* SRC_ASSEMBLER_SOLVER_SOLVER_H_ */
