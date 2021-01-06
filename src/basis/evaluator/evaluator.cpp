
#include "evaluator.h"
#include "basis/utilities/parser.h"

using namespace espreso;

bool Evaluator::isConstant() const
{
	return parameters.empty();
}

bool Evaluator::isCoordinateDependent() const
{
	for (auto it = parameters.begin(); it != parameters.end(); ++it) {
		if (StringCompare::caseInsensitiveEq(*it, "X") || StringCompare::caseInsensitiveEq(*it, "Y") || StringCompare::caseInsensitiveEq(*it, "Z")) {
			return true;
		}
	}
	return false;
}

bool Evaluator::isTimeDependent() const
{
	for (auto it = parameters.begin(); it != parameters.end(); ++it) {
		if (StringCompare::caseInsensitiveEq(*it, "TIME")) {
			return true;
		}
	}
	return false;
}

bool Evaluator::isFrequencyDependent() const
{
	for (auto it = parameters.begin(); it != parameters.end(); ++it) {
		if (StringCompare::caseInsensitiveEq(*it, "FREQUENCY")) {
			return true;
		}
	}
	return false;
}

bool Evaluator::isTemperatureDependent() const
{
	for (auto it = parameters.begin(); it != parameters.end(); ++it) {
		if (StringCompare::caseInsensitiveEq(*it, "TEMPERATURE")) {
			return true;
		}
	}
	return false;
}
