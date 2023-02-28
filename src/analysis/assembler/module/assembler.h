
#ifndef SRC_ANALYSIS_ASSEMBLER_MODULE_ASSEMBLER_H_
#define SRC_ANALYSIS_ASSEMBLER_MODULE_ASSEMBLER_H_

#include "basis/evaluator/evaluator.h"
#include "analysis/assembler/operator.h"
#include "math/primitives/vector_sparse.h"

#include "analysis/scheme/steadystate.h"

#include <vector>
#include <map>
#include <functional>

namespace espreso {

struct ECFExpression;
struct ECFExpressionVector;
struct ECFExpressionOptionalVector;
class ConvectionConfiguration;
struct ImpedanceConfiguration;
struct PointSourceConfiguration;
struct PhysicsConfiguration;

class Assembler
{
public:
	Assembler(PhysicsConfiguration &settings);
	virtual ~Assembler();

	static Evaluator* getEvaluator(size_t interval, std::map<std::string, ECFExpression> &settings);
	static Evaluator* getEvaluator(size_t interval, std::map<std::string, ECFExpressionVector> &settings, int dim);

	PhysicsConfiguration &settings;
	std::vector<int> etype, bfilter;
	std::vector<std::vector<int> > btype;
	std::vector<std::vector<ActionOperator*> > elementOps;
	std::vector<std::vector<std::vector<ActionOperator*> > > boundaryOps;

protected:
	double assemble(ActionOperator::Action action);
	virtual double instantiate(ActionOperator::Action action, int code, int etype, const std::vector<ActionOperator*> &ops, size_t interval, esint elements) { return 0; }

	template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype>
	double loop(ActionOperator::Action action, const std::vector<ActionOperator*> &ops, esint elements);

	bool checkExpression(const std::string &name, ECFExpression &expression);
	bool checkElementParameter(const std::string &name, std::map<std::string, ECFExpression> &settings);
	bool checkElementParameter(const std::string &name, std::map<std::string, ECFExpressionVector> &settings);
	bool checkElementParameter(const std::string &name, std::map<std::string, ECFExpressionVector> &settings, int dim);
	bool checkBoundaryParameter(const std::string &name, std::map<std::string, ECFExpression> &settings);
	bool checkBoundaryParameter(const std::string &name, std::map<std::string, ECFExpressionVector> &settings);
	bool checkBoundaryParameter(const std::string &name, std::map<std::string, ECFExpressionOptionalVector> &settings);

	void printElementVolume(std::vector<double> &volume);
	void printBoundarySurface(std::vector<double> &surface);

	void printMaterials(const std::map<std::string, std::string> &settings);
	template<typename Ttype>
	void validateRegionSettings(const std::string &name, const std::map<std::string, Ttype> &settings);

	// TODO: remove
	Matrix_Base<double> *K;
};

template <typename T>
static void reset(T *t)
{
	if (t) {
		t->set(0);
	}
}

template <typename T, typename ...Other>
static void reset(T *t, Other... other)
{
	reset(t);
	reset(other...);
}

template <typename T>
static void update(T *t)
{
	if (t) {
		t->synchronize();
		t->touched = true;
	}
}

template <typename T, typename ...Other>
static void update(T *t, Other... other)
{
	update(t);
	update(other...);
}

}

#endif /* SRC_ANALYSIS_ASSEMBLER_MODULE_ASSEMBLER_H_ */
