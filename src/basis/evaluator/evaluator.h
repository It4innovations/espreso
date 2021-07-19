
#ifndef SRC_BASIS_EVALUATOR_EVALUATOR_H_
#define SRC_BASIS_EVALUATOR_EVALUATOR_H_

#include "esinfo/stepinfo.h"
#include "math/simd/simd.h"

#include <string>
#include <vector>

namespace espreso {

class Evaluator {

protected:
	enum class Type: int {
		DEFAULT,
		CONST,
		EXPRESSION,
		TABLE,
		TABLE_INTERPOLATION,
		ARRAY
	};

public:
	struct Params {
		struct General {
			double *val;
			int offset, increment;
		};

		int _ncoors;
		const double* _coors;
		const double* _inittemp;
		const double* _temp;
		const double* _disp;
		double _time;
		double _frequency;

		std::vector<General> general;

		Params(): _ncoors(0), _coors(NULL), _inittemp(NULL), _temp(NULL), _disp(NULL), _time(step::time.current), _frequency(step::frequency.current) {}

		int ncoords() { return _ncoors; }
		const double* coords() { return _coors; }
		const double* inittemp() { return _inittemp; }
		const double* temp() { return _temp; }
		double time() { return _time; }
		double freq() { return _frequency; }
		const double* disp() { return _disp; }

		Params& coords(int ncoords, const double *coords) { _ncoors = ncoords; _coors = coords; return *this; }
		Params& inittemp(const double *inittemp) { _inittemp = inittemp; return *this; }
		Params& temp(const double *temp) { _temp = temp; return *this; }
		Params& time(double  time) { _time = time; return *this; }
		Params& freq(double  freq) { _frequency = freq; return *this; }
		Params& disp(const double *disp) { _disp = disp; return *this; }
	};

	virtual Type type() { return Type::DEFAULT; }
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

	bool isConstant() const;
	bool isCoordinateDependent() const;
	bool isTimeDependent() const;
	bool isFrequencyDependent() const;
	bool isTemperatureDependent() const;

	Params params;
	std::vector<std::string> variables;
	bool isset;
};

}



#endif /* SRC_BASIS_EVALUATOR_EVALUATOR_H_ */
