
#include "material.h"
#include "coordinates.h"
#include "../settings/evaluator.h"
#include "../../config/configuration.h"

using namespace espreso;

Material::Material(const Coordinates &coordinates, const Configuration &configuration)
: _coordinates(coordinates)
{
	_values.resize(static_cast<int>(MATERIAL_PARAMETER::SIZE), NULL);
	if (configuration.parameters.find("MODEL") != configuration.parameters.end()) {
		_model = (MATERIAL_MODEL)configuration.parameters.find("MODEL")->second->option();
	} else {
		_model = MATERIAL_MODEL::SIZE;
	}

	for (size_t p = 0; p < configuration.orderedParameters.size(); p++) {
		if (!StringCompare::caseInsensitiveEq(configuration.orderedParameters[p]->name, "MODEL")) {
			const std::string &value = configuration.orderedParameters[p]->get();
			if (StringCompare::contains(value, "xyzt")) {
				_values.push_back(new CoordinatesEvaluator(value, _coordinates));
			} else {
				espreso::Expression expr(value, {});
				_values.push_back(new ConstEvaluator(expr.evaluate({})));
			}
		}
	}
}

Material::~Material()
{
	for (size_t i = 0; i < _values.size(); i++) {
		if (_values[i] != NULL) {
			delete _values[i];
		}
	}
}

void Material::set(size_t index, const std::string &value)
{
	if (_values[index] != NULL) {
		delete _values[index];
	}

	if (StringCompare::contains(value, "xyzt")) {
		_values[index] = new CoordinatesEvaluator(value, _coordinates);
	} else {
		espreso::Expression expr(value, {});
		_values[index] = new ConstEvaluator(expr.evaluate({}));
	}
}



