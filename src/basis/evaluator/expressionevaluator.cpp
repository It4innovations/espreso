
#include "esinfo/envinfo.h"
#include "expressionevaluator.h"

#include "basis/containers/point.h"
#include "basis/utilities/parser.h"
#include "basis/expression/expression.h"
#include "omp.h"

using namespace espreso;

ExpressionEvaluator::ExpressionEvaluator(const std::string &expression)
{
	size_t threads = info::env::OMP_NUM_THREADS;

	_expressions.resize(info::env::OMP_NUM_THREADS);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		_expressions[t] = new Expression(expression, ExpressionEvaluator::variables());
	}

	_coordinateDependency = StringCompare::contains(expression, { "X", "Y", "Z" });
	_timeDependency = StringCompare::contains(expression, { "TIME" });
	_frequencyDependency = StringCompare::contains(expression, { "FREQUENCY" });
	_temperatureDependency = StringCompare::contains(expression, { "INITIAL_TEMPERATURE", "TEMPERATURE" });
}

ExpressionEvaluator::~ExpressionEvaluator()
{
	size_t threads = info::env::OMP_NUM_THREADS;

	for (size_t t = 0; t < threads; t++) {
		delete _expressions[t];
	}
}

ExpressionEvaluator::ExpressionEvaluator(const ExpressionEvaluator &other)
{
	size_t threads = info::env::OMP_NUM_THREADS;

	_expressions.resize(info::env::OMP_NUM_THREADS);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		_expressions[t] = new Expression(*other._expressions[t]);
	}

	_coordinateDependency = other._coordinateDependency;
	_timeDependency = other._timeDependency;
	_frequencyDependency = other._frequencyDependency;
	_temperatureDependency = other._temperatureDependency;
}

void ExpressionEvaluator::evalVector(esint size, esint increment, const Params &params, double *results) const
{
	int thread = omp_get_thread_num();
	for (esint i = 0; i < size; ++i) {
		if (params._coors != NULL) {
			_expressions[thread]->values[0] = params._coors[i * params._ncoors + 0];
			_expressions[thread]->values[1] = params._coors[i * params._ncoors + 1];
			_expressions[thread]->values[2] = params._ncoors == 3 ? params._coors[i * params._ncoors + 2] : 0;
		}
		if (params._inittemp != NULL) {
			_expressions[thread]->values[3] = params._inittemp[i];
		}
		if (params._temp != NULL) {
			_expressions[thread]->values[4] = params._temp[i];
		}
		_expressions[thread]->values[5] = params._time;
		_expressions[thread]->values[6] = params._frequency;
		if (params._disp != NULL) {
			_expressions[thread]->values[8] = params._disp[i];
		}
		results[i * increment] = _expressions[thread]->evaluate();
	}
}

void ExpressionEvaluator::evalFiltered(esint size, esint increment, const esint *elements, const esint *distribution, const Params &params, double *results) const
{
	int thread = omp_get_thread_num();
	for (esint i = 0; i < size; ++i) {
		for (esint e = distribution[elements[i]]; e < distribution[elements[i] + 1]; ++e) {
			if (params._coors != NULL) {
				_expressions[thread]->values[0] = params._coors[e * params._ncoors + 0];
				_expressions[thread]->values[1] = params._coors[e * params._ncoors + 1];
				_expressions[thread]->values[2] = params._ncoors == 3 ? params._coors[e * params._ncoors + 2] : 0;
			}
			if (params._inittemp != NULL) {
				_expressions[thread]->values[3] = params._inittemp[e];
			}
			if (params._temp != NULL) {
				_expressions[thread]->values[4] = params._temp[e];
			}
			_expressions[thread]->values[5] = params._time;
			_expressions[thread]->values[6] = params._frequency;
			if (params._disp != NULL) {
				_expressions[thread]->values[8] = params._disp[i];
			}
			results[e * increment] = _expressions[thread]->evaluate();
		}
	}
}

void ExpressionEvaluator::evalSelectedSparse(esint size, esint increment, const esint *selection, const Params &params, double *results) const
{
	int thread = omp_get_thread_num();
	for (esint i = 0; i < size; ++i) {
		if (params._coors != NULL) {
			_expressions[thread]->values[0] = params._coors[selection[i] * params._ncoors + 0];
			_expressions[thread]->values[1] = params._coors[selection[i] * params._ncoors + 1];
			_expressions[thread]->values[2] = params._ncoors == 3 ? params._coors[selection[i] * params._ncoors + 2] : 0;
		}
		if (params._inittemp != NULL) {
			_expressions[thread]->values[3] = params._inittemp[selection[i]];
		}
		if (params._temp != NULL) {
			_expressions[thread]->values[4] = params._temp[selection[i]];
		}
		_expressions[thread]->values[5] = params._time;
		_expressions[thread]->values[6] = params._frequency;
		if (params._disp != NULL) {
			_expressions[thread]->values[8] = params._disp[i];
		}
		results[i * increment] = _expressions[thread]->evaluate();
	}
}

void ExpressionEvaluator::evalSelectedDense(esint size, esint increment, const esint *selection, const Params &params, double *results) const
{
	int thread = omp_get_thread_num();
	for (esint i = 0; i < size; ++i) {
		if (params._coors != NULL) {
			_expressions[thread]->values[0] = params._coors[selection[i] * params._ncoors + 0];
			_expressions[thread]->values[1] = params._coors[selection[i] * params._ncoors + 1];
			_expressions[thread]->values[2] = params._ncoors == 3 ? params._coors[selection[i] * params._ncoors + 2] : 0;
		}
		if (params._inittemp != NULL) {
			_expressions[thread]->values[3] = params._inittemp[selection[i]];
		}
		if (params._temp != NULL) {
			_expressions[thread]->values[4] = params._temp[selection[i]];
		}
		_expressions[thread]->values[5] = params._time;
		_expressions[thread]->values[6] = params._frequency;
		if (params._disp != NULL) {
			_expressions[thread]->values[8] = params._disp[i];
		}
		results[selection[i] * increment] = _expressions[thread]->evaluate();
	}
}

double ExpressionEvaluator::evaluate(double r) const
{
	int thread = omp_get_thread_num();
	_expressions[thread]->values[7] = r;
	return _expressions[thread]->evaluate();
}

std::string ExpressionEvaluator::getEXPRTKForm() const
{
	return _expressions.front()->expression();
}




