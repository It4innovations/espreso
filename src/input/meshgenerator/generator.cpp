
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
	return "";
}

static void setProperty(
		const espreso::Mesh &mesh,
		std::vector<espreso::Evaluator*> &evaluators,
		std::vector<espreso::Element*> &elements,
		const std::string &values,
		const std::vector<std::string> &parameters,
		const std::vector<espreso::Property> &properties)
{
	for (size_t p = 0; p < properties.size(); p++) {
		std::string value = parseValue((p < parameters.size()) ? parameters[p] : "ALL", values);
		if (value.size()) {
			if (value.find("xyz") == std::string::npos) {
				espreso::Expression expr(value, {});
				evaluators.push_back(new espreso::ConstEvaluator(expr.evaluate({})));
			} else {
				evaluators.push_back(new espreso::CoordinatesEvaluator(value, mesh.coordinates()));
			}
			for (size_t i = 0; i < elements.size(); i++) {
				elements[i]->settings(properties[p]).push_back(evaluators.back());
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

	std::vector<Element*> selection;

	for (auto it = values.begin(); it != values.end(); ++it) {
		if (checkInterval(_settings.nodes, it->first)) {
			pickNodesInInterval(nodes, selection, _settings.nodes.find(it->first)->second);
			setProperty(mesh, evaluators, selection, it->second, parameters, properties);
		}
		if (checkInterval(_settings.faces, it->first)) {
			generateFacesInInterval(selection, _settings.faces.find(it->first)->second);
			setProperty(mesh, evaluators, selection, it->second, parameters, properties);
			faces.insert(faces.end(), selection.begin(), selection.end());
		}
		if (checkInterval(_settings.edges, it->first)) {
			generateEdgesInInterval(selection, _settings.edges.find(it->first)->second);
			setProperty(mesh, evaluators, selection, it->second, parameters, properties);
			edges.insert(edges.end(), selection.begin(), selection.end());
		}
		if (checkInterval(_settings.elements, it->first)) {
			if (_settings.elements.find(it->first)->second.all()) {
				setProperty(mesh, evaluators, elements, it->second, parameters, properties);
			} else {
				pickElementsInInterval(elements, selection, _settings.elements.find(it->first)->second);
				setProperty(mesh, evaluators, elements, it->second, parameters, properties);
			}

		}
	}
}


