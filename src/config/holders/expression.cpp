
#include "expression.h"
#include "config/configuration.hpp"
#include "esinfo/eslog.hpp"
#include "basis/utilities/parser.h"
#include "basis/utilities/utils.h"
#include "basis/expression/expression.h"
#include "basis/evaluator/constevaluator.h"
#include "basis/evaluator/expressionevaluator.h"
#include "basis/evaluator/tableinterpolationevaluator.h"

#include "mesh/store/elementsregionstore.h"
#include "mesh/store/boundaryregionstore.h"

using namespace espreso;

bool ECFExpression::forall(const std::map<std::string, ECFExpression> &parameter, std::function<bool(const ECFExpression &expr)> fnc)
{
	for (auto it = parameter.begin(); it != parameter.end(); ++it) {
		if (!fnc(it->second)) {
			return false;
		}
	}
	return true;
}

bool ECFHarmonicExpression::forall(const std::map<std::string, ECFHarmonicExpression> &parameter, std::function<bool(const ECFExpression &expr)> fnc)
{
	for (auto it = parameter.begin(); it != parameter.end(); ++it) {
		if (!fnc(it->second.magnitude)) {
			return false;
		}
		if (!fnc(it->second.phase)) {
			return false;
		}
	}
	return true;
}

bool ECFExpressionVector::forall(const std::map<std::string, ECFExpressionVector> &parameter, std::function<bool(const ECFExpression &expr)> fnc)
{
	for (auto it = parameter.begin(); it != parameter.end(); ++it) {
		if (!fnc(it->second.x) || !fnc(it->second.y) || !fnc(it->second.z)) {
			return false;
		}
	}
	return true;
}

bool ECFHarmonicExpressionVector::forall(const std::map<std::string, ECFHarmonicExpressionVector> &parameter, std::function<bool(const ECFExpression &expr)> fnc)
{
	for (auto it = parameter.begin(); it != parameter.end(); ++it) {
		if (!fnc(it->second.magnitude.x) || !fnc(it->second.magnitude.y) || !fnc(it->second.magnitude.z)) {
			return false;
		}
		if (!fnc(it->second.phase.x) || !fnc(it->second.phase.y) || !fnc(it->second.phase.z)) {
			return false;
		}
	}
	return true;
}

bool ECFExpressionOptionalVector::forall(const std::map<std::string, ECFExpressionOptionalVector> &parameter, std::function<bool(const ECFExpression &expr)> fnc)
{
	for (auto it = parameter.begin(); it != parameter.end(); ++it) {
		if (it->second.all.isSet()) {
			if (!fnc(it->second.all)) {
				return false;
			}
		} else {
			if (!fnc(it->second.x) || !fnc(it->second.y) || !fnc(it->second.z)) {
				return false;
			}
		}
	}
	return true;
}

ECFExpression::ECFExpression(const std::vector<std::string> &variables)
: variables(variables), evaluator(NULL)
{
	createEvaluator();
}

ECFExpression::ECFExpression(const std::vector<std::string> &variables, const std::string &initialValue)
: value(initialValue), variables(variables), evaluator(NULL)
{
	createEvaluator();
}

ECFExpression::~ECFExpression()
{
	if (evaluator) {
		delete evaluator;
	}
}

ECFExpression::ECFExpression(const ECFExpression &other)
{
	value = other.value;
	variables = other.variables;
	evaluator = NULL;
	if (other.evaluator != NULL) {
		evaluator = other.evaluator->copy();
	}
}

ECFExpression& ECFExpression::operator=(const ECFExpression &other)
{
	if (this != &other) {
		value = other.value;
		variables = other.variables;
		if (evaluator != NULL) {
			delete evaluator;
			evaluator = NULL;
		}
		if (other.evaluator != NULL) {
			evaluator = other.evaluator->copy();
		}
	}
	return *this;
}

void ECFExpression::createEvaluator()
{
	if (evaluator != NULL) {
		delete evaluator;
		evaluator = NULL;
	}
	if (value == "") {
		evaluator = new Evaluator();
		return;
	}
	if (StringCompare::contains(value, { "TABULAR" })) {
		std::string value = Parser::strip(value.substr(value.find_first_of("[")));
		value = value.substr(1, value.size() - 3);
		std::vector<std::string> lines = Parser::split(value, ";");
		std::vector<std::pair<double, double> > table;

		for (size_t i = 0; i < lines.size(); i++) {
			if (lines[i].size() == 0) {
				continue;
			}
			std::vector<std::string> line = Parser::split(lines[i], ",");
			if (line.size() != 2) {
				eslog::globalerror("Invalid TABULAR data: %s\n", value.c_str());
			}
			table.push_back(std::make_pair(std::stod(line[0]), std::stod(line[1])));
		}
		evaluator = new TableInterpolationEvaluator(table);
		return;
	}
	if (StringCompare::contains(value, ExpressionEvaluator::variables())) {
		if (Expression::isValid(value, variables)) {
			evaluator = new ExpressionEvaluator(value);
			return;
		}
	} else {
		if (Expression::isValid(value, variables)) {
			Expression expr(value, {});
			evaluator = new ConstEvaluator(expr.evaluate());
			return;
		}
	}
}

void ECFHarmonicExpression::init()
{
	type = Type::COMPONENTS;
	REGISTER(type, ECFMetaData()
			.setdescription({ "Description type." })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("COMPONENTS").setdescription("Separated components.")));

	REGISTER(magnitude, ECFMetaData().setdescription({ "Magnitude." }));
	REGISTER(phase, ECFMetaData().setdescription({ "Phase." }));
}


ECFHarmonicExpression::ECFHarmonicExpression(const std::vector<std::string> &variables)
: type(Type::COMPONENTS), magnitude(variables), phase(variables)
{
	init();
}

ECFHarmonicExpression::ECFHarmonicExpression(const std::vector<std::string> &variables, const std::string &initialValue)
: type(Type::COMPONENTS), magnitude(variables, initialValue), phase(variables, initialValue)
{
	init();
}

void ECFHarmonicExpressionVector::init()
{
	type = Type::COMPONENTS;
	REGISTER(type, ECFMetaData()
			.setdescription({ "Description type." })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("COMPONENTS").setdescription("Separated components.")));

	REGISTER(magnitude, ECFMetaData().setdescription({ "Magnitude." }));
	REGISTER(phase, ECFMetaData().setdescription({ "Phase." }));
}

ECFHarmonicExpressionVector::ECFHarmonicExpressionVector(DIMENSION *dimension, const std::vector<std::string> &variables)
: type(Type::COMPONENTS), magnitude(dimension, variables), phase(dimension, variables)
{
	init();
}

ECFHarmonicExpressionVector::ECFHarmonicExpressionVector(DIMENSION *dimension, const std::vector<std::string> &variables, const std::string &initialValue)
: type(Type::COMPONENTS), magnitude(dimension, variables, initialValue), phase(dimension, variables, initialValue)
{
	init();
}

void ECFExpressionVector::init()
{
	REGISTER(x, ECFMetaData()
			.setdescription({ "x-direction." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.allowonly([&] () { return *dimension != DIMENSION::Z; }));
	REGISTER(y, ECFMetaData()
			.setdescription({ "y-direction." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.allowonly([&] () { return *dimension == DIMENSION::D2 || *dimension == DIMENSION::D3; }));
	REGISTER(z, ECFMetaData()
			.setdescription({ "z-direction." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.allowonly([&] () { return *dimension == DIMENSION::Z || *dimension == DIMENSION::D3; }));
}

ECFExpressionVector::ECFExpressionVector(const ECFExpressionVector &other)
: data{other.data[0], other.data[1], other.data[2]}, dimension(other.dimension)
{

}

ECFExpressionVector& ECFExpressionVector::operator=(const ECFExpressionVector &other)
{
	if (this != &other) {
		data[0] = other.data[0];
		data[1] = other.data[1];
		data[2] = other.data[2];
		dimension = other.dimension;
	}
	return *this;
}

ECFExpressionVector::ECFExpressionVector(DIMENSION *dimension, const std::vector<std::string> &variables)
: data{variables, variables, variables}, dimension(dimension)
{
	init();
}

ECFExpressionVector::ECFExpressionVector(DIMENSION *dimension, const std::vector<std::string> &variables, const std::string &initialValue)
: data{ {variables,initialValue}, {variables,initialValue}, {variables,initialValue}}, dimension(dimension)
{
	init();
}

ECFExpressionOptionalVector::ECFExpressionOptionalVector(DIMENSION *dimension, const std::vector<std::string> &variables)
: ECFExpressionVector(dimension, variables), all(variables)
{
	REGISTER(all, ECFMetaData()
			.setdescription({ "all-directions." })
			.setdatatype({ ECFDataType::EXPRESSION }));
}


