
#ifndef SRC_BASIS_EVALUATOR_CONSTEVALUATOR_H_
#define SRC_BASIS_EVALUATOR_CONSTEVALUATOR_H_

#include "evaluator.h"

namespace espreso {

class ConstEvaluator: public Evaluator {

public:
	ConstEvaluator(double value): _value(value) {}

	Type type() { return Type::CONST; }
	virtual Evaluator* copy() const { return new ConstEvaluator(*this); }

	void evalVector(esint size, esint increment, const Params &params, double *results) const;
	void evalFiltered(esint size, esint increment, const esint *elements, const esint *distribution, const Params &params, double *results) const;

	void evalSelectedSparse(esint size, esint increment, const esint *selection, const Params &params, double *results) const
	{
		evalVector(size, increment, params, results);
	}

	void evalSelectedDense(esint size, esint increment, const esint *selection, const Params &params, double *results) const;

	double evaluate(double r) const { return _value; }

	std::string getEXPRTKForm() const { return std::to_string(_value); }

protected:
	double _value;
};

}


#endif /* SRC_BASIS_EVALUATOR_CONSTEVALUATOR_H_ */
