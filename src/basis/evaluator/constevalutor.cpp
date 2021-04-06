
#include "constevaluator.h"

using namespace espreso;

void ConstEvaluator::evalVector(esint size, esint increment, const Params &params, double *results) const
{
	for (esint i = 0; i < size; ++i) {
		results[i * increment] = _value;
	}
}


void ConstEvaluator::evalFiltered(esint size, esint increment, const esint *elements, const esint *distribution, const Params &params, double *results) const
{
	for (esint i = 0; i < size; ++i) {
		for (esint e = distribution[elements[i]]; e < distribution[elements[i] + 1]; ++e) {
			results[e * increment] = _value;
		}
	}
}

void ConstEvaluator::evalSelectedDense(esint size, esint increment, const esint *selection, const Params &params, double *results) const
{
	for (esint i = 0; i < size; ++i) {
		results[selection[i] * increment] = _value;
	}
}


