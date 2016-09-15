
#include "evaluator.h"

using namespace espreso;

Evaluator* Evaluator::create(std::ifstream &is, const Coordinates &coordinates)
{
	Evaluator::Type type;
	is.read(reinterpret_cast<char *>(&type), sizeof(Evaluator::Type));

	switch (type) {
	case Type::DEFAULT: return new Evaluator();
	case Type::CONST: return new ConstEvaluator(is);
	case Type::COORDINATE: return new CoordinatesEvaluator(is, coordinates);
	default: ESINFO(GLOBAL_ERROR) << "Unknown evaluator type"; return NULL;
	}
}


