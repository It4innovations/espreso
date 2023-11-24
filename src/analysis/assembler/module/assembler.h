
#ifndef SRC_ANALYSIS_ASSEMBLER_MODULE_ASSEMBLER_H_
#define SRC_ANALYSIS_ASSEMBLER_MODULE_ASSEMBLER_H_

#include "basis/evaluator/evaluator.h"
#include "math/primitives/vector_sparse.h"

#include <vector>
#include <map>
#include <functional>

namespace espreso {

struct ECFExpression;
struct ECFExpressionVector;
struct ECFExpressionOptionalVector;
struct ConvectionConfiguration;
struct ImpedanceConfiguration;
struct PointSourceConfiguration;
struct PhysicsConfiguration;

class Assembler
{
public:
	enum Action: int {
		VOID       = 1 << 0,
		ASSEMBLE   = 1 << 1,
		REASSEMBLE = 1 << 2,
		FILL       = 1 << 3,
		SOLUTION   = 1 << 4,
		ITERATION  = 1 << 5,
		PREPROCESS = 1 << 6
	};

	Assembler(PhysicsConfiguration &settings);
	virtual ~Assembler();

	static ECFExpression* getExpression(size_t interval, std::map<std::string, ECFExpression> &settings);
	static ECFExpressionVector* getExpression(size_t interval, std::map<std::string, ECFExpressionVector> &settings);
	static ECFExpression* getExpression(const std::string &name, std::map<std::string, ECFExpression> &settings);
	static ECFExpressionVector* getExpression(const std::string &name, std::map<std::string, ECFExpressionVector> &settings);

	static Evaluator* getEvaluator(size_t interval, std::map<std::string, ECFExpression> &settings);
	static Evaluator* getEvaluator(size_t interval, std::map<std::string, ECFExpressionVector> &settings, int dim);

	PhysicsConfiguration &settings;

protected:
	void assemble(Action action);
	virtual void run(Action action, size_t interval) =0;
	virtual void run(Action action, size_t region, size_t interval) =0;

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
};

inline Assembler::Action  operator| (Assembler::Action  a1, Assembler::Action a2) { return static_cast<Assembler::Action>(static_cast<int>(a1) | static_cast<int>(a2)); }
inline Assembler::Action  operator& (Assembler::Action  a1, Assembler::Action a2) { return static_cast<Assembler::Action>(static_cast<int>(a1) & static_cast<int>(a2)); }
inline Assembler::Action& operator|=(Assembler::Action &a1, Assembler::Action a2) { a1 = a1 | a2; return a1; }
inline Assembler::Action& operator&=(Assembler::Action &a1, Assembler::Action a2) { a1 = a1 & a2; return a1; }

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
