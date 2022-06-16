
#include "esinfo/envinfo.h"
#include "expressionevaluator.h"

#include "basis/containers/point.h"
#include "basis/utilities/parser.h"
#include "basis/expression/expression.h"
#include "omp.h"

using namespace espreso;

ExpressionEvaluator::ExpressionEvaluator(const std::string &expression, std::vector<std::string> &variables)
{
	this->variables = variables;
	_expressions.resize(info::env::OMP_NUM_THREADS);

	#pragma omp parallel for
	for (int t = 0; t < info::env::OMP_NUM_THREADS; t++) {
		_expressions[t] = new Expression(expression, variables);
	}
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
	variables = other.variables;
}

void ExpressionEvaluator::evalVectorInc(esint size, esint increment, const Params &params, double *results) const
{
	int thread = omp_get_thread_num();
	for (esint i = 0; i < size; ++i) {
		for (size_t p = 0; p < params.general.size(); ++p) {
			_expressions[thread]->values[p] = *(params.general[p].val + i * params.general[p].increment + params.general[p].offset);
		}
		results[i * increment] = _expressions[thread]->evaluate();
	}
}

void ExpressionEvaluator::evalVectorSimdInc(esint size, esint increment, const Params &params, double *results) const
{
	int thread = omp_get_thread_num();
	for (esint i = 0; i < size; i+=SIMD::size) {
		for (esint simdLane = 0; simdLane < SIMD::size; ++simdLane) {
			for (size_t p = 0; p < params.general.size(); ++p) {
				_expressions[thread]->values[p] = *(params.general[p].val + i * params.general[p].increment + params.general[p].offset + simdLane);
			}
			results[i * increment + simdLane] = _expressions[thread]->evaluate();
		}
	}
}

void ExpressionEvaluator::evalFilteredInc(esint size, esint increment, const esint *elements, const esint *distribution, const Params &params, double *results) const
{
	int thread = omp_get_thread_num();
	for (esint i = 0; i < size; ++i) {
		for (esint e = distribution[elements[i]]; e < distribution[elements[i] + 1]; ++e) {
			for (size_t p = 0; p < params.general.size(); ++p) {
				_expressions[thread]->values[p] = *(params.general[p].val + e * params.general[p].increment + params.general[p].offset);
			}
			results[e * increment] = _expressions[thread]->evaluate();
		}
	}
}

void ExpressionEvaluator::evalSelectedSparseInc(esint size, esint increment, const esint *selection, const Params &params, double *results) const
{
	int thread = omp_get_thread_num();
	for (esint i = 0; i < size; ++i) {
		for (size_t p = 0; p < params.general.size(); ++p) {
			_expressions[thread]->values[p] = *(params.general[p].val + selection[i] * params.general[p].increment + params.general[p].offset);
		}
		results[i * increment] = _expressions[thread]->evaluate();
	}
}

void ExpressionEvaluator::evalSelectedDenseInc(esint size, esint increment, const esint *selection, const Params &params, double *results) const
{
	int thread = omp_get_thread_num();
	for (esint i = 0; i < size; ++i) {
		for (size_t p = 0; p < params.general.size(); ++p) {
					_expressions[thread]->values[p] = *(params.general[p].val + selection[i] * params.general[p].increment + params.general[p].offset);
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




