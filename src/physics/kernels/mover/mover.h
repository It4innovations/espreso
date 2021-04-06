
#ifndef SRC_PHYSICS_KERNELS_MOVER_MOVER_H_
#define SRC_PHYSICS_KERNELS_MOVER_MOVER_H_

#include "basis/containers/point.h"
#include "basis/evaluator/evaluator.h"
#include "config/holders/expression.h"
#include <vector>
#include <map>

namespace espreso {

template <typename TEBoundaries, typename TEData> class serializededata;
struct Element;
struct NodeData;
struct ElementData;
class MaterialConfiguration;
class MoverParameter;
template <class ...TArgs> struct Move;

struct MoverIncrement {
	int nodes, elements;

	MoverIncrement(): nodes(0), elements(0) {}
	MoverIncrement(int nodes, int elements)
	: nodes(nodes), elements(elements) {}

	int operator*(const MoverIncrement &other)
	{
		return nodes * other.nodes + elements * other.elements;
	}
};

struct InputCoordinates {
	int dimension;
	serializededata<esint, Point> *values;

	InputCoordinates(): dimension(0), values(NULL) {}

	bool isConst() const { return false; }
};

struct InputExpression {
	ECFExpression *ecf;
	Evaluator::Params *params;

	InputExpression(): ecf(NULL), params(NULL) {}

	bool isConst() const { return ecf->evaluator->isConstant(); }
};

struct InputHarmonicMagnitudeExpression {
	ECFHarmonicExpression *ecf;
	Evaluator::Params *params;

	InputHarmonicMagnitudeExpression(): ecf(NULL), params(NULL) {}

	bool isConst() const { return ecf->magnitude.evaluator->isConstant(); }
};

struct InputHarmonicPhaseExpression {
	ECFHarmonicExpression *ecf;
	Evaluator::Params *params;

	InputHarmonicPhaseExpression(): ecf(NULL), params(NULL) {}

	bool isConst() const { return ecf->phase.evaluator->isConstant(); }
};

struct InputExpressionVector {
	ECFExpressionVector *ecf;
	Evaluator::Params *params;

	InputExpressionVector(): ecf(NULL), params(NULL) {}

	bool isConst() const
	{
		return
				ecf->x.evaluator->isConstant() &&
				ecf->y.evaluator->isConstant() &&
				ecf->z.evaluator->isConstant();
	}
};

struct InputHarmonicMagnitudeExpressionVector {
	ECFHarmonicExpressionVector *ecf;
	Evaluator::Params *params;

	InputHarmonicMagnitudeExpressionVector(): ecf(NULL), params(NULL) {}

	bool isConst() const
	{
		return
				ecf->magnitude.x.evaluator->isConstant() &&
				ecf->magnitude.y.evaluator->isConstant() &&
				ecf->magnitude.z.evaluator->isConstant();
	}
};

struct InputHarmonicPhaseExpressionVector {
	ECFHarmonicExpressionVector *ecf;
	Evaluator::Params *params;

	InputHarmonicPhaseExpressionVector(): ecf(NULL), params(NULL) {}

	bool isConst() const
	{
		return
				ecf->phase.x.evaluator->isConstant() &&
				ecf->phase.y.evaluator->isConstant() &&
				ecf->phase.z.evaluator->isConstant();
	}
};

struct InputExpressionOptionalVector {
	ECFExpressionOptionalVector *ecf;
	Evaluator::Params *params;

	InputExpressionOptionalVector(): ecf(NULL), params(NULL) {}

	bool isConst() const
	{
		return
				ecf->all.isSet() && ecf->all.evaluator->isConstant() &&
				ecf->x.isSet() && ecf->x.evaluator->isConstant() &&
				ecf->y.isSet() && ecf->y.evaluator->isConstant() &&
				ecf->z.isSet() && ecf->z.evaluator->isConstant();
	}
};

struct InputExpressionMap {
	std::map<std::string, ECFExpression> *ecf;
	Evaluator::Params *params;

	InputExpressionMap(): ecf(NULL), params(NULL) {}

	bool isConst() const
	{
		return ecf->size() <= 1 && ECFExpression::forall(*ecf, [] (const ECFExpression &expr) {
			return expr.evaluator->isConstant();
		});
	}
};

struct InputExpressionVectorMap {
	std::map<std::string, ECFExpressionVector> *ecf;
	Evaluator::Params *params;

	InputExpressionVectorMap(): ecf(NULL), params(NULL) {}

	bool isConst() const
	{
		return ecf->size() <= 1 && ECFExpressionVector::forall(*ecf, [] (const ECFExpression &expr) {
			return expr.evaluator->isConstant();
		});
	}
};

struct InputExpressionOptionalVectorMap {
	std::map<std::string, ECFExpressionOptionalVector> *ecf;
	Evaluator::Params *params;

	InputExpressionOptionalVectorMap(): ecf(NULL), params(NULL) {}

	bool isConst() const
	{
		return ecf->size() <= 1 && ECFExpressionOptionalVector::forall(*ecf, [] (const ECFExpression &expr) {
			return expr.evaluator->isConstant();
		});
	}
};

struct OutputNodes {
	NodeData *data;

	OutputNodes(): data(NULL) {}
};

struct OutputElements {
	ElementData *data;

	OutputElements(): data(NULL) {}
};

struct ElementNodeValues {
	int dimension;
	int isConst;
	serializededata<esint, esint> *nodes;
	serializededata<esint, double> *values;

	ElementNodeValues()
	: dimension(0), isConst(0), nodes(NULL), values(NULL) {}

	void clear();
};

struct Mover {
	virtual ~Mover() {}
	virtual void operator()() =0;
};

struct MoverInstance {
	Element* element;

	Evaluator::Params nodeparams;
	Evaluator::Params kernelparams;

	std::vector<Mover*> stepMovers;
	std::vector<Mover*> solutionMovers;
	std::vector<MoverParameter*> inputs;
	std::vector<MoverParameter*> outputs;

	void clear();
	void initKernelIterator(esint offset);
	void initOutputIterator(esint offset);
	void next(const MoverIncrement &counters);

	void nextSubstep();
	void solutionChanged();

	MoverInstance();
	MoverInstance(const MoverInstance &other);
	~MoverInstance();

	void registerSolution(MoverParameter &parameter, const MoverIncrement &increment);
	void registerInput(MoverParameter &parameter, const MoverIncrement &increment);
	void registerOutput(MoverParameter &parameter, const MoverIncrement &increment);

	template <class ...TArgs>
	void afterStepChange(TArgs ...args) { solutionMovers.push_back(new Move<TArgs...>(args...)); }

	template <class ...TArgs>
	void afterSolutionChange(TArgs ...args) { solutionMovers.push_back(new Move<TArgs...>(args...)); }

	template <class ...TArgs>
	void now(TArgs ...args) { Move<TArgs...>(args...)(); }
};

struct ElementIterator: public MoverInstance {
	const MaterialConfiguration* material;
	esint offset;

	ElementIterator(): material(NULL), offset(0) {}
};

struct BoundaryIterator: public MoverInstance {
	BoundaryIterator() {}
};

#define MOVE(FROM, TO)                    \
	template <>                           \
	struct Move<FROM, TO>: public Mover { \
	FROM from; TO to; bool moved;         \
	Move(const FROM &from, const TO &to); \
	void operator()();                    \
	void now() { operator()(); }          \
	};                                    \

MOVE(InputCoordinates                      , ElementNodeValues)
MOVE(InputExpression                       , ElementNodeValues)
MOVE(InputExpressionVector                 , ElementNodeValues)
MOVE(InputExpressionOptionalVector         , ElementNodeValues)
MOVE(InputExpressionMap                    , OutputNodes)
MOVE(InputExpressionMap                    , ElementNodeValues)
MOVE(InputExpressionVectorMap              , ElementNodeValues)
MOVE(InputExpressionOptionalVectorMap      , ElementNodeValues)
MOVE(InputHarmonicMagnitudeExpression      , ElementNodeValues)
MOVE(InputHarmonicPhaseExpression          , ElementNodeValues)
MOVE(InputHarmonicMagnitudeExpressionVector, ElementNodeValues)
MOVE(InputHarmonicPhaseExpressionVector    , ElementNodeValues)
MOVE(OutputNodes                           , ElementNodeValues)
MOVE(OutputElements                        , ElementNodeValues)
MOVE(ElementNodeValues                     , OutputNodes)
MOVE(ElementNodeValues                     , OutputElements)

}

#endif /* SRC_PHYSICS_KERNELS_MOVER_MOVER_H_ */
