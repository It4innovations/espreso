
#ifndef SRC_CONFIGURATION_CONFIGURATION_H_
#define SRC_CONFIGURATION_CONFIGURATION_H_

#include <string>
#include <vector>
#include <map>

#include "../basis/utilities/parser.h"

namespace espreso {

struct ParameterBase {
	std::string name;
	std::string description;

	ParameterBase() {};
	ParameterBase(const std::string &name, const std::string &description);
};

struct Parameter: public ParameterBase {
	std::string allowedValue;

	Parameter(const std::string &name, const std::string &description, const std::string &allowedValue);

	virtual bool set(const std::string &value) =0;
	virtual std::string get() const =0;
	virtual size_t index() const;

	virtual std::string XMLAttributeType() const =0;
	virtual void XMLChildsElements(std::ostream &os, size_t indent) const {};

	virtual ~Parameter();
};

struct Configuration: public ParameterBase {
	bool copy;
	std::map<std::string, Parameter*, StringCompare> parameters;
	std::map<std::string, Configuration*, StringCompare> subconfigurations;
	std::vector<Parameter*> orderedParameters;
	std::vector<Configuration*> orderedSubconfiguration;
	std::vector<Configuration*> toDelete;

	Configuration();
	Configuration(const Configuration &other);

	Configuration& operator=(const Configuration &other);

	template <typename Ttype>
	static Ttype create(const std::string &name, const std::string &description, Configuration* conf)
	{
		Ttype configuration;
		conf->subconfigurations[name] = &configuration;
		conf->orderedSubconfiguration.push_back(&configuration);
		configuration.name = name;
		configuration.description = description;
		return configuration;
	}

	virtual Configuration& operator[](const std::string &subconfiguration);
	virtual bool set(const std::string &parameter, const std::string &value);

	virtual const std::vector<Configuration*>& storeConfigurations() const { return orderedSubconfiguration; }
	virtual const std::vector<Parameter*>& storeParameters() const { return orderedParameters; }
	virtual const Configuration* configurationPattern() const { return NULL; }
	virtual const Parameter* parameterPattern() const { return NULL; }

	virtual ~Configuration();
};

}



#endif /* SRC_CONFIGURATION_CONFIGURATION_H_ */
