
#ifndef SRC_CONFIG_HOLDERS_EXPRESSION_H_
#define SRC_CONFIG_HOLDERS_EXPRESSION_H_

#include "config/metadata.h"
#include "config/description.h"
#include <vector>
#include <map>
#include <functional>

namespace espreso {

class Evaluator;
struct BoundaryRegionStore;
struct ElementsRegionStore;

struct ECFExpression {
	enum class Scope {
		GLOBAL,
		ELEMENT,
		ENODES,
		EGPS,
		BNODES,
		BGPS,
		NODE
	};

	static std::vector<ECFExpression*> parametrized; // parametrized expressions

	std::string value;
	Scope scope;
	std::vector<std::string> parameters;
	Evaluator *evaluator;
	bool isset;

	ECFExpression(Scope scope);
	ECFExpression(const std::string &initialValue, Scope scope);
	ECFExpression(const ECFExpression &other);
	ECFExpression& operator=(const ECFExpression &other);
	~ECFExpression();

	void createEvaluator();

	static bool forall(const std::map<std::string, ECFExpression> &parameter, std::function<bool(const ECFExpression &expr)> fnc);
};

struct ECFHarmonicExpression: public ECFDescription {
	enum class Type {
		COMPONENTS
	};

	Type type;
	ECFExpression magnitude, phase;

	ECFHarmonicExpression();
	ECFHarmonicExpression(const std::string &initialValue);

	static bool forall(const std::map<std::string, ECFHarmonicExpression> &parameter, std::function<bool(const ECFExpression &expr)> fnc);

protected:
	void init();
};

struct ECFExpressionOptionalVector;

struct ECFExpressionVector: public ECFDescription {
	ECFExpression data[3];
	ECFExpression &x = data[0], &y = data[1], &z = data[2];
	DIMENSION *dimension;

	ECFExpressionVector(const ECFExpressionVector &other);
	ECFExpressionVector& operator=(const ECFExpressionVector &other);
	ECFExpressionVector(DIMENSION *dimension, ECFExpression::Scope scope);
	ECFExpressionVector(DIMENSION *dimension, const std::string &initialValue, ECFExpression::Scope scope);

	static bool forall(const std::map<std::string, ECFExpressionVector> &parameter, std::function<bool(const ECFExpression &expr)> fnc);

protected:
	void init();
};

struct ECFHarmonicExpressionVector: public ECFDescription {
	enum class Type {
		COMPONENTS
	};

	Type type;
	ECFExpressionVector magnitude, phase;

	ECFHarmonicExpressionVector(DIMENSION *dimension, ECFExpression::Scope scope);
	ECFHarmonicExpressionVector(DIMENSION *dimension, const std::string &initialValue, ECFExpression::Scope scope);

	static bool forall(const std::map<std::string, ECFHarmonicExpressionVector> &parameter, std::function<bool(const ECFExpression &expr)> fnc);

protected:
	void init();
};

struct ECFExpressionOptionalVector: public ECFExpressionVector {
	ECFExpression all;

	ECFExpressionOptionalVector(DIMENSION *dimension, ECFExpression::Scope scope);

	static bool forall(const std::map<std::string, ECFExpressionOptionalVector> &parameter, std::function<bool(const ECFExpression &expr)> fnc);
};

}

#endif /* SRC_CONFIG_HOLDERS_EXPRESSION_H_ */
