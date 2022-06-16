
#ifndef SRC_BASIS_EVALUATOR_CONSTEVALUATOR_H_
#define SRC_BASIS_EVALUATOR_CONSTEVALUATOR_H_

#include "evaluator.h"

namespace espreso {

class ConstEvaluator: public Evaluator {

public:
	ConstEvaluator(double value): _value(value) {}

	virtual Evaluator* copy() const { return new ConstEvaluator(*this); }

	void evalVectorInc(esint size, esint increment, const Params &params, double *results) const;
	void evalVectorSimdInc(esint size, esint increment, const Params &params, double *results) const;
	void evalFilteredInc(esint size, esint increment, const esint *elements, const esint *distribution, const Params &params, double *results) const;

	void evalSelectedSparseInc(esint size, esint increment, const esint *selection, const Params &params, double *results) const
	{
		evalVectorInc(size, increment, params, results);
	}

	void evalSelectedDenseInc(esint size, esint increment, const esint *selection, const Params &params, double *results) const;

	double evaluate(double r) const { return _value; }

	std::string getEXPRTKForm() const { return std::to_string(_value); }

protected:
	double _value;
};

}


#endif /* SRC_BASIS_EVALUATOR_CONSTEVALUATOR_H_ */
