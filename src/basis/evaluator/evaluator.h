
#ifndef SRC_BASIS_EVALUATOR_EVALUATOR_H_
#define SRC_BASIS_EVALUATOR_EVALUATOR_H_

#include "esinfo/stepinfo.h"
#include "math/simd/simd.h"

#include <string>
#include <vector>

namespace espreso {

struct Variable;

class Evaluator {

public:
	struct Params {
		struct General {
			double *val;
			int offset, increment;
			Variable *variable;
		};

		std::vector<General> general;
	};

	virtual Evaluator* copy() const { return new Evaluator(*this); }

	virtual ~Evaluator() {};

	double eval() const
	{
		return eval(params);
	}

	double eval(const Params &params) const
	{
		double value;
		evalVector(1, params, &value);
		return value;
	}

	virtual void evalVector(esint size, const Params &params, double *results) const
	{
		evalVector(size, 1, params, results);
	}

	virtual void evalVector(esint size, esint increment, const Params &params, double *results) const
	{
		for (esint i = 0; i < size; ++i) {
			results[i * increment] = 0;
		}
	}

	virtual void evalVectorSimd(esint size, esint increment, const Params &params, double *results) const
	{
		evalVector(size, 1, params, results);
	}

	virtual void evalFiltered(esint size, const esint *elements, const esint *distribution, const Params &params, double *results) const
	{
		evalFiltered(size, 1, elements, distribution, params, results);
	}

	virtual void evalFiltered(esint size, esint increment, const esint *elements, const esint *distribution, const Params &params, double *results) const
	{
		for (esint i = 0; i < size; ++i) {
			for (esint e = distribution[elements[i]]; e < distribution[elements[i] + 1]; ++e) {
				results[e * increment] = 0;
			}
		}
	}

	virtual void evalSelectedSparse(esint size, const esint *selection, const Params &params, double *results) const
	{
		evalSelectedSparse(size, 1, selection, params, results);
	}

	virtual void evalSelectedSparse(esint size, esint increment, const esint *selection, const Params &params, double *results) const
	{
		for (esint i = 0; i < size; ++i) {
			results[i * increment] = 0;
		}
	}

	virtual void evalSelectedDense(esint size, const esint *selection, const Params &params, double *results) const
	{
		evalSelectedDense(size, 1, selection, params, results);
	}

	virtual void evalSelectedDense(esint size, esint increment, const esint *selection, const Params &params, double *results) const
	{
		for (esint i = 0; i < size; ++i) {
			results[selection[i] * increment] = 0;
		}
	}

	virtual double evaluate(double r) const { return 0; }

	virtual std::string getEXPRTKForm() const { return ""; }

	Params params;
	std::vector<std::string> variables;
	bool isset;
};

}



#endif /* SRC_BASIS_EVALUATOR_EVALUATOR_H_ */
