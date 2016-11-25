
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


void Generator::materials(std::vector<Material> &materials)
{
	materials.resize(2, Material(mesh.coordinates()));

	for (auto it = _settings.material1.begin(); it != _settings.material1.end(); ++it) {
		if (!materials[0].setParameter(it->first, it->second)) {
			ESINFO(ALWAYS) << Info::TextColor::YELLOW << "Unknown material parameter " << it->first;
		}
	}
	for (auto it = _settings.material2.begin(); it != _settings.material2.end(); ++it) {
		if (!materials[1].setParameter(it->first, it->second)) {
			ESINFO(ALWAYS) << Info::TextColor::YELLOW << "Unknown material parameter " << it->first;
		}
	}
}

void Generator::regions(
			std::vector<Evaluator*> &evaluators,
			std::vector<Region> &regions,
			std::vector<Element*> &elements,
			std::vector<Element*> &faces,
			std::vector<Element*> &edges,
			std::vector<Element*> &nodes)
{
	this->loadProperties(evaluators, elements, faces, edges, nodes, "DIRICHLET", { "T", "x", "y", "z" }, { Property::TEMPERATURE, Property::DISPLACEMENT_X, Property::DISPLACEMENT_Y, Property::DISPLACEMENT_Z });
	this->loadProperties(evaluators, elements, faces, edges, nodes, "NEUMAN", { "P"}, { Property::PRESSURE });
	this->loadProperties(evaluators, elements, faces, edges, nodes, "HEAT_SOURCES", { "T" }, { Property::HEAT_SOURCE });
	this->loadProperties(evaluators, elements, faces, edges, nodes, "TRANSLATION_MOTIONS", { "x", "y", "z" }, { Property::TRANSLATION_MOTION_X, Property::TRANSLATION_MOTION_Y, Property::TRANSLATION_MOTION_Z });
	this->loadProperties(evaluators, elements, faces, edges, nodes, "ACCELERATION", { "x", "y", "z" }, { Property::ACCELERATION_X, Property::ACCELERATION_Y, Property::ACCELERATION_Z });
	this->loadProperties(evaluators, elements, faces, edges, nodes, "THICKNESS", { }, { Property::THICKNESS });
	this->loadProperties(evaluators, elements, faces, edges, nodes, "INITIAL_TEMPERATURE", { }, { Property::INITIAL_TEMPERATURE });
	this->loadProperties(evaluators, elements, faces, edges, nodes, "TEMPERATURE", { }, { Property::TEMPERATURE });
	this->loadProperties(evaluators, elements, faces, edges, nodes, "OBSTACLE", { }, { Property::OBSTACLE });
	this->loadProperties(evaluators, elements, faces, edges, nodes, "NORMAL_DIRECTION", { }, { Property::NORMAL_DIRECTION });
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
			if (value.find("x") == std::string::npos && value.find("y") == std::string::npos && value.find("z") == std::string::npos && value.find("t") == std::string::npos) {
				espreso::Expression expr(value, {});
				evaluators.push_back(new espreso::ConstEvaluator(expr.evaluate({}), properties[p]));
			} else {
				evaluators.push_back(new espreso::CoordinatesEvaluator(value, mesh.coordinates(), properties[p]));
			}
			for (size_t i = 0; i < elements.size(); i++) {
				elements[i]->addSettings(properties[p], evaluators.back());
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
		std::vector<Element*> selection;
		if (checkInterval(_settings.nodes, it->first)) {
			if (_settings.nodes.find(it->first)->second.all()) {
				setProperty(mesh, evaluators, nodes, it->second, parameters, properties);
			} else {
				pickNodesInInterval(nodes, selection, _settings.nodes.find(it->first)->second);
				setProperty(mesh, evaluators, selection, it->second, parameters, properties);
			}
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


