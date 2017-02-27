
#include "configuration.h"

#include "../basis/logging/logging.h"

#include <algorithm>

using namespace espreso;

ParameterBase::ParameterBase(const std::string &name, const std::string &description, const std::string &allowedValue)
: name(name), description(description), allowedValue(allowedValue)
{

}

size_t ParameterBase::index() const
{
	ESINFO(ERROR) << "Call index of non-indexable parameter";
	return 0;
}

ParameterBase::~ParameterBase()
{

}

Configuration::Configuration(): copy(false)
{

}

Configuration::Configuration(const Configuration &other)
: copy(true),
  parameters(other.parameters),
  subconfigurations(other.subconfigurations),
  orderedParameters(other.orderedParameters),
  orderedSubconfiguration(other.orderedSubconfiguration),
  name(other.name),
  description(other.description)
{

}

Configuration& Configuration::operator=(const Configuration &other)
{
	if (this != &other) {
		copy = true;
		parameters = other.parameters;
		subconfigurations = other.subconfigurations;
		orderedParameters = other.orderedParameters;
		orderedSubconfiguration = other.orderedSubconfiguration;
		name = other.name;
		description = other.description;
	}
	return *this;
}

Configuration& Configuration::operator[](const std::string &subconfiguration)
{
	if (subconfigurations.find(subconfiguration) != subconfigurations.end()) {
		return *subconfigurations.find(subconfiguration)->second;
	} else {
		ESINFO(GLOBAL_ERROR) << "Unknown sub-configuration '" << subconfiguration << "'";
		exit(EXIT_FAILURE);
	}
}

bool Configuration::set(const std::string &parameter, const std::string &value)
{
	if (parameters.find(parameter) != parameters.end()) {
		return parameters.find(parameter)->second->set(value);
	} else {
		ESINFO(GLOBAL_ERROR) << "Unknown parameter '" << parameter << "'";
		return false;
	}
}

Configuration::~Configuration()
{
	if (!copy) {
		for (auto it = parameters.begin(); it != parameters.end(); ++it) {
			delete it->second;
		}
		std::for_each(toDelete.begin(), toDelete.end(), [] (Configuration *c) { delete c; });
	}
}




