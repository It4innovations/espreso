
#include "generator.h"

using namespace espreso::input;

static bool checkInterval(const std::map<std::string, espreso::Interval> &intervals, const std::string &name)
{
	if (intervals.find(name) != intervals.end()) {
		return true;
	}
	return false;
}

static std::string parseValue(const std::string &param, const std::string &value)
{
	std::vector<std::string> values = espreso::Parser::split(espreso::Parser::strip(value), ",");
	if (values.size() == 1 && !espreso::Parser::contains(values[0], ":")) {
		return values[0];
	}
	for (size_t i = 0; i < values.size(); i++) {
		if (espreso::StringCompare::caseInsensitiveEq(espreso::Parser::getParameter(values[i], ":"), param)) {
			return espreso::Parser::getValue(values[i], ":");
		}
	}
	return "0";
}

static void setElements(
		std::vector<espreso::Evaluator*> &evaluators,
		std::vector<espreso::Element*> &elements,
		const espreso::Interval &interval,
		const std::string &values,
		const std::vector<std::string> &parameters,
		const std::vector<espreso::Property> &properties)
{
	if (interval.all()) {
		for (size_t p = 0; p < properties.size(); p++) {
			std::string value = parseValue((p < parameters.size()) ? parameters[p] : "ALL", values);
			evaluators.push_back(new espreso::ExpressionEvaluator(value));
			for (size_t i = 0; i < elements.size(); i++) {
				elements[i]->settings(properties[p]).push_back(evaluators.back());
			}
		}
	} else {
		ESINFO(espreso::GLOBAL_ERROR) << "Only selection of all elements is supported.";
	}
}


static void setFaces(
		std::vector<espreso::Evaluator*> &evaluators,
		std::vector<espreso::Element*> &faces,
		const espreso::Interval &interval,
		const std::string &values,
		const std::vector<std::string> &parameters,
		const std::vector<espreso::Property> &properties)
{
	ESINFO(espreso::GLOBAL_ERROR) << "Faces selection is not implemented.";
}

static void setNodes(
		std::vector<espreso::Evaluator*> &evaluators,
		std::vector<espreso::Element*> &nodes,
		const espreso::Coordinates &coordinates,
		const espreso::Interval &interval,
		const std::string &values,
		const std::vector<std::string> &parameters,
		const std::vector<espreso::Property> &properties)
{
	for (size_t p = 0; p < properties.size(); p++) {
		std::string value = parseValue((p < parameters.size()) ? parameters[p] : "ALL", values);
		evaluators.push_back(new espreso::ExpressionEvaluator(value));
		for (size_t i = 0; i < coordinates.clusterSize(); i++) {
			if (interval.isIn(coordinates[i])) {
				nodes[i]->settings(properties[p]).push_back(evaluators.back());
			}
		}
	}
}

void Generator::loadProperties(
		std::vector<Evaluator*> &evaluators,
		std::vector<Element*> &elements,
		std::vector<Element*> &faces,
		std::vector<Element*> &edges,
		std::vector<Element*> &nodes,
		const std::string &name,
		std::vector<std::string> parameters,
		std::vector<Property> properties)
{
	if (_settings.properties.find(name) == _settings.properties.end()) {
		return;
	}

	const std::map<std::string, std::string> &values = _settings.properties.find(name)->second;

	for (auto it = values.begin(); it != values.end(); ++it) {
		if (checkInterval(_settings.nodes, it->first)) {
			setNodes(evaluators, nodes, mesh.coordinates(), _settings.nodes.find(it->first)->second, it->second, parameters, properties);
		}
		if (checkInterval(_settings.faces, it->first)) {
			setFaces(evaluators, faces, _settings.faces.find(it->first)->second, it->second, parameters, properties);
		}
		if (checkInterval(_settings.elements, it->first)) {
			setElements(evaluators, elements, _settings.elements.find(it->first)->second, it->second, parameters, properties);
		}
	}
}


