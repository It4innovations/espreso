
#ifndef SRC_BASIS_EVALUATOR_EXPRESSIONEVALUATOR_H_
#define SRC_BASIS_EVALUATOR_EXPRESSIONEVALUATOR_H_

#include "evaluator.h"

#include <string>
#include <vector>

namespace espreso {

class Expression;

/**
 * Evaluator can be called from various threads.
 * Create one instance for each worker in order to avoid race conditions.
 */
class ExpressionEvaluator: public Evaluator {

public:
//	static std::vector<std::string> variables() { return { "X", "Y", "Z", "INITIAL_TEMPERATURE", "TEMPERATURE", "TIME", "FREQUENCY", "R", "DISPLACEMENT" }; }

	ExpressionEvaluator(const std::string &expression, std::vector<std::string> &variables);
	ExpressionEvaluator(const ExpressionEvaluator &other);
	~ExpressionEvaluator();

	virtual Evaluator* copy() const { return new ExpressionEvaluator(*this); }

	void evalVectorInc(esint size, esint increment, const Params &params, double *results) const;
	void evalVectorSimdInc(esint size, esint increment, const Params &params, double *results) const;
	void evalFilteredInc(esint size, esint increment, const esint *elements, const esint *distribution, const Params &params, double *results) const;
	void evalSelectedSparseInc(esint size, esint increment, const esint *selection, const Params &params, double *results) const;
	void evalSelectedDenseInc(esint size, esint increment, const esint *selection, const Params &params, double *results) const;

	double evaluate(double r) const;

	std::string getEXPRTKForm() const;
	std::string toString() const { return getEXPRTKForm(); }

protected:
	std::vector<Expression*> _expressions;
};

}


#endif /* SRC_BASIS_EVALUATOR_EXPRESSIONEVALUATOR_H_ */
