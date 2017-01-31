
#include "material.h"
#include "coordinates.h"
#include "../settings/evaluator.h"
#include "../../config/configuration.h"

using namespace espreso;

Material::Material(const Coordinates &coordinates): _coordinates(coordinates), _models((size_t)PHYSICS::SIZE, MATERIAL_MODEL::SIZE)
{
	for (size_t i = 0; i < (size_t)MATERIAL_PARAMETER::SIZE; i++) {
		_values.push_back(new ConstEvaluator(0));
	}
}

Material::Material(const Coordinates &coordinates, const Configuration &configuration)
: _coordinates(coordinates)
{
	for (size_t i = 0; i < (size_t)MATERIAL_PARAMETER::SIZE; i++) {
		_values.push_back(new ConstEvaluator(0));
	}
	if (configuration.parameters.find("MODEL") != configuration.parameters.end()) {
		_models.resize((size_t)PHYSICS::SIZE, (MATERIAL_MODEL)configuration.parameters.find("MODEL")->second->index());
	} else {
		_models.resize((size_t)PHYSICS::SIZE, MATERIAL_MODEL::SIZE);
	}

	for (size_t p = 0; p < configuration.orderedParameters.size(); p++) {
		if (!StringCompare::caseInsensitiveEq(configuration.orderedParameters[p]->name, "MODEL")) {
			const std::string &value = configuration.orderedParameters[p]->get();
			delete _values[configuration.orderedParameters[p]->index()];
			if (StringCompare::contains(value, "xyzt")) {
				_values[configuration.orderedParameters[p]->index()] = new CoordinatesEvaluator(value, _coordinates);
			} else {
				espreso::Expression expr(value, {});
				_values[configuration.orderedParameters[p]->index()] = new ConstEvaluator(expr.evaluate({}));
			}
		}
	}
}

Material::~Material()
{
	for (size_t i = 0; i < _values.size(); i++) {
		delete _values[i];
	}
}

void Material::set(MATERIAL_PARAMETER parameter, const std::string &value)
{
	delete _values[(size_t)parameter];

	if (StringCompare::contains(value, "xyzt")) {
		_values[(size_t)parameter] = new CoordinatesEvaluator(value, _coordinates);
	} else {
		espreso::Expression expr(value, {});
		_values[(size_t)parameter] = new ConstEvaluator(expr.evaluate({}));
	}
}

void Material::set(MATERIAL_PARAMETER parameter, Evaluator* value)
{
	delete _values[(size_t)parameter];
	_values[(size_t)parameter] = value;
}

void Material::store(std::ofstream& os)
{
	int size = _models.size();
	os.write(reinterpret_cast<const char *>(&size), sizeof(int));
	for (size_t i = 0; i < _models.size(); i++) {
		int model = (int)_models[i];
		os.write(reinterpret_cast<const char *>(&model), sizeof(int));
	}

	size = _values.size();
	os.write(reinterpret_cast<const char *>(&size), sizeof(int));
	for (size_t i = 0; i < _values.size(); i++) {
		_values[i]->store(os);
	}
}

void Material::load(std::ifstream& is)
{
	int data;
	is.read(reinterpret_cast<char *>(&data), sizeof(int));
	_models.resize(data, MATERIAL_MODEL::SIZE);
	ESTEST(MANDATORY) << (_models.size() == (size_t)MATERIAL_MODEL::SIZE ? ESPRESOTest::TEST_PASSED : ESPRESOTest::TEST_PASSED) << "Cannot read old binary format.";
	for (size_t i = 0; i < _models.size(); i++) {
		is.read(reinterpret_cast<char *>(&data), sizeof(int));
		_models[i] = (MATERIAL_MODEL)data;
	}

	for (size_t i = 0; i < _values.size(); i++) {
		delete _values[i];
	}
	is.read(reinterpret_cast<char *>(&data), sizeof(int));
	_values.resize(data);
	ESTEST(MANDATORY) << (_values.size() == (size_t)MATERIAL_PARAMETER::SIZE ? ESPRESOTest::TEST_PASSED : ESPRESOTest::TEST_PASSED) << "Cannot read old binary format.";
	for (size_t i = 0; i < _values.size(); i++) {
		_values[i] = Evaluator::create(is, _coordinates);
	}
}



