
#include "constevaluator.h"

using namespace espreso;

void ConstEvaluator::evalVectorInc(esint size, esint increment, const Params &params, double *results) const
{
	for (esint i = 0; i < size; ++i) {
		results[i * increment] = _value;
	}
}

void ConstEvaluator::evalVectorSimdInc(esint size, esint increment, const Params &params, double *results) const
{
	for (esint i = 0; i < size; i += SIMD::size) {
		for (esint simdLane = 0; simdLane < SIMD::size; ++simdLane) {
			results[i * increment + simdLane] = _value;
		}
	}
}

void ConstEvaluator::evalFilteredInc(esint size, esint increment, const esint *elements, const esint *distribution, const Params &params, double *results) const
{
	for (esint i = 0; i < size; ++i) {
		for (esint e = distribution[elements[i]]; e < distribution[elements[i] + 1]; ++e) {
			results[e * increment] = _value;
		}
	}
}

void ConstEvaluator::evalSelectedDenseInc(esint size, esint increment, const esint *selection, const Params &params, double *results) const
{
	for (esint i = 0; i < size; ++i) {
		results[selection[i] * increment] = _value;
	}
}


