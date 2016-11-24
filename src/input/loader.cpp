
#include "loader.h"

#include "../mesh/settings/property.h"
#include "../config/description.h"

using namespace espreso::input;

void Loader::boundaryConditions()
{
	auto getValueIndex = [] (const std::vector<std::string> &values, const std::string &parameter) -> size_t {
		if (values.size() == 1 && !Parser::contains(values[0], ":=")) {
			return 0;
		}
		for (size_t i = 0; i < values.size(); i++) {
			if (StringCompare::caseInsensitiveEq(parameter, Parser::strip(Parser::split(values[i], ":=")[0]))) {
				return i;
			}
		}
		return values.size();
	};

	auto loadProperty = [&] (const std::map<std::string, std::string> &regions, const std::vector<std::string> &parameters, const std::vector<Property> &properties) {
		for (auto it = regions.begin(); it != regions.end(); ++it) {
			Region &region = mesh.region(it->first);
			std::vector<std::string> values = Parser::split(it->second, ",;");

			for (size_t p = 0; p < properties.size(); p++) {
				size_t index = getValueIndex(values, parameters[p]);
				if (index < values.size()) {
					std::string value = Parser::contains(values[index], ":=") ? Parser::split(values[index], ":=")[1] : values[index];
					if (value.find("xyzt") == std::string::npos) {
						espreso::Expression expr(value, {});
						mesh._evaluators.push_back(new espreso::ConstEvaluator(expr.evaluate({}), properties[p]));
					} else {
						mesh._evaluators.push_back(new espreso::CoordinatesEvaluator(value, mesh.coordinates(), properties[p]));
					}
					for (size_t i = 0; i < region.elements.size(); i++) {
						region.elements[i]->addSettings(properties[p], mesh._evaluators.back());
					}
				}
			}
		}
	};

	loadProperty(configuration.displacement.values, { "x", "y", "z" }, { Property::DISPLACEMENT_X, Property::DISPLACEMENT_Y, Property::DISPLACEMENT_Z });
	loadProperty(configuration.normal_presure.values, { "p" }, { Property::PRESSURE });

}



