
#ifndef SRC_ASSEMBLER_PHYSICSSOLVER_ASSEMBLER_H_
#define SRC_ASSEMBLER_PHYSICSSOLVER_ASSEMBLER_H_

#include <functional>
#include <vector>

namespace espreso {

struct Step;
struct Instance;
struct Physics;
class Mesh;
struct LinearSolver;
class TimeEval;
class SparseMatrix;
class Solution;
class Store;
enum Matrices: int;
enum class ElementType;
enum class SumOperation;
enum class SumRestriction;

class Assembler {

public:
	Assembler(Instance &instance, Physics &physics, Mesh &mesh, Store &store, LinearSolver &linearSolver);
	~Assembler();

	void preprocessData(const Step &step);
	void updateMatrices(const Step &step, Matrices matrices);
	void processSolution(const Step &step);

	void setRegularizationCallback();
	void setEmptyRegularizationCallback();
	void setB0Callback();

	void solve(const Step &step, Matrices updatedMatrices);

	void storeSolution(const Step &step);
	void storeSubSolution(const Step &step);

	void finalize();

	Solution* addSolution(const std::string &name, ElementType eType);

	/// z = a * x + b + y
	void sum(std::vector<std::vector<double> > &z, double a, const std::vector<std::vector<double> > &x, double b, const std::vector<std::vector<double> > &y, const std::string &description);
	/// z = a * x + b + y (prefix variant)
	void sum(std::vector<std::vector<double> > &z, double a, const std::vector<std::vector<double> > &x, double b, const std::vector<std::vector<double> > &y, const std::vector<size_t> &prefix, const std::string &description);
	/// A += beta * B
	void sum(std::vector<SparseMatrix> &A, double beta, std::vector<SparseMatrix> &B, const std::string &description);

	/// y = A * x
	void multiply(std::vector<std::vector<double> > &y, std::vector<SparseMatrix> &A, std::vector<std::vector<double> > &x, const std::string &description);

	double sumSquares(const Step &step, const std::vector<std::vector<double> > &data, SumOperation operation, SumRestriction restriction, const std::string &description);
	void addToDirichletInB1(double a, const std::vector<std::vector<double> > &x);
	double maxAbsValue(const std::vector<std::vector<double> > &v, const std::string &description);
	double lineSearch(const Step &step, const std::vector<std::vector<double> > &U, std::vector<std::vector<double> > &deltaU, std::vector<std::vector<double> > &F_ext);

	Instance &instance;
	Physics &physics;
	Mesh &mesh;
	Store &store;
	LinearSolver &linearSolver;

protected:
	void timeWrapper(const std::string &action, std::function<void(void)> operations);

	bool checkForStore(const std::string &name);
	void storeMatrices(Matrices matrices, size_t domain);
	void storeWrapper(const std::string &name, Matrices matrices);
	void storeWrapper(const std::string &name, Matrices matrices, size_t domain);
	void storeWrapper(const std::string &name, std::vector<SparseMatrix> &matrices);
	void storeWrapper(const std::string &name, std::vector<std::vector<double> > &data);

	TimeEval *_timeStatistics;
};

}



#endif /* SRC_ASSEMBLER_PHYSICSSOLVER_ASSEMBLER_H_ */
